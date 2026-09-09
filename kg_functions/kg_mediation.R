# ==============================================================================
# kg_functions/kg_mediation.R  --  kgMediation: does a MEDIATOR carry an EXPOSURE's
# effect on an OUTCOME? (op#13)
#
#   exposure    a metadata COLUMN holding the exposure levels (e.g. "final_cluster",
#               "diagnosis"). The caller resolves the surface form -> column, exactly as
#               F35/kgPhenoGroup does; this function never guesses a column.
#   levels      "A,B" -- the TWO levels compared (e.g. "C2,C0"). Two levels only: a
#               proportion mediated is defined for ONE contrast, not for an omnibus F.
#   mediators   named set ("MAP1LC3A" | "MAP1LC3A,PRSS1,..."). NEVER "all" -- see maxMediators.
#   phenotypes  a leaf metadata column | a PARENT concept ("gsis") expanded to its child
#               columns by .kgExpandPhenotypes | a comma list. Continuous outcomes only.
#   layers      "all" (each mediator's own layers, resolver-driven) | comma list of keys.
#   covariates  "default" (the 5) | "none" | list. METADATA COLUMNS ONLY, as op1 -- the
#               MEDIATOR is not a covariate here, it is the tested term.
#   minN        donor floor, default 30. RAISED FROM 10 (user ruling 2026-08-10) because
#               10 was a floor for a ONE-model op and this is a THREE-model design: measured
#               on real data, a pseudobulk pair reached n=12 against ~8 parameters, leaving
#               df=4, and the PM it produced carried a bootstrap CI of [-16, +34]. A number
#               that wide is not a result.
#   minDF       residual-df floor on the DIRECT model, default 10, applied ON TOP of minN.
#               minN alone cannot protect the fit, because the parameter count moves with the
#               covariate set: n=30 with 8 parameters is not the same design as n=30 with 3.
#               THE FLOORS SKIP A PAIR; they never quietly weaken it.
#   boot        0 (STAGE 1: point estimates only) | B (STAGE 2: bootstrap CI on PM).
#   maxPairs    refuse above this many (mediator x layer x outcome) pairs. NOT a silent
#               truncation -- returns RES-TOO-MANY with the arithmetic so the caller can
#               tell the user which factor to cut.
#
# THE MODEL -- three fits per pair, on ONE complete-case donor set:
#     a-path   mediator ~ exposure + covariates          -> a_coef, a_p
#     total    outcome  ~ exposure + covariates          -> b_total
#     direct   outcome  ~ exposure + mediator + covars   -> b_direct
#     indirect = b_total - b_direct ;  PM = indirect / b_total
#
# WHY THE THREE FITS LIVE IN ONE FUNCTION AND NOT IN THREE op1 CALLS. The complete-case
# set MOVES with the mediator: the protein layer covers ~234 donors where RNA-seq covers
# ~371, so a `total` model fitted without the mediator would be fitted on donors the
# `direct` model never sees, and a difference of coefficients across two donor sets is not
# a proportion of anything. All three fits here are restricted to the SAME rows, and the
# row count is reported per pair so the reader can check it.
#
# WHY NOT .kgLmCell. It fits ONE model and derives its OWN complete cases, which is the
# exact thing that must not happen three times independently here. The design matrix is
# built once per pair and solved with .lm.fit -- measured 219us per pair at n=371 versus
# ~4.5ms through lm() formulas, which is what makes the stage-1 sweep sub-second.
#
# THE LIMIT THIS FUNCTION CANNOT ESCAPE, and the caller must carry it into the answer:
# exposure, mediator and outcome are measured in the SAME donors at the SAME time. A
# non-zero PM is CONSISTENT WITH mediation; it does not establish it, and it cannot be
# distinguished from confounding by a shared cause or from reverse direction. Nothing in
# this output should be phrased as "the mediator explains X% of the effect".
#
# Writes kg_mediation.csv (one row per mediator x layer x phenotype). Returns
#   "RES-OK;<n_rows>" | "RES-TOO-MANY;<pairs>;<mediators>;<layers>;<outcomes>;<cap>"
#   | "RES-NO" | "RES-NO-DB" | "RES-NO-META" | "RES-NO-PHENO" | "RES-NO-FEATURE"
#   | "RES-NO-RESOLVE" | "RES-NO-RESOLVER" | "RES-NO-EXPOSURE" | "RES-NO-LEVELS".
# ==============================================================================

.kgMedDefaultCovs <- c("donorage", "donorsex", "bodymassindex",
                       "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# Solve one design matrix; return the coefficient and p for column `which`.
# .lm.fit gives coefficients and residuals but no inference, so the SE is formed from the
# residual variance and (X'X)^-1 -- the textbook OLS quantities, not an approximation.
.kgMedFit <- function(X, y, which){
  f <- try(.lm.fit(X, y), silent = TRUE)
  if(inherits(f, "try-error") || any(is.na(f$coefficients))) return(NULL)
  n <- nrow(X); p <- ncol(X); df <- n - p
  if(df < 1) return(NULL)
  s2 <- sum(f$residuals^2) / df
  xtxi <- try(chol2inv(chol(crossprod(X))), silent = TRUE)
  if(inherits(xtxi, "try-error")) return(NULL)
  se <- sqrt(s2 * diag(xtxi))[which]
  b  <- f$coefficients[which]
  if(!is.finite(b) || !is.finite(se) || se <= 0) return(NULL)
  list(b = b, se = se, p = 2 * stats::pt(abs(b / se), df, lower.tail = FALSE), df = df)
}

.kgMediationImpl <- function(exposure, levels, mediators, phenotypes,
                        layers = "all", covariates = "default", use_raw = FALSE,
                        minN = 30, boot = 0, seed = 1,
                        maxPairs = 1000, maxMediators = 150, minDF = 10,
                        subset = "all"){   # `subset` LAST: positional callers unaffected
  .kgSetPaths()
  suppressPackageStartupMessages({ library(DBI); library(RSQLite) })

  # ---- EVERY NUMERIC ARRIVES AS A STRING OVER Rserve ---------------------------
  # KgRCenter assigns character values (kgpg_minn = "10"), so an uncoerced `minN`
  # would make `length(don) < minN` a STRING comparison -- "9" < "30" is FALSE, and the
  # floor would silently invert. Coerced once, here, with the documented default on any
  # unparseable value rather than a guess.
  .int <- function(x, d){ v <- suppressWarnings(as.integer(x)); if(is.na(v)) d else v }
  minN         <- .int(minN, 30L);   minDF        <- .int(minDF, 10L)
  boot         <- .int(boot, 0L);    seed         <- .int(seed, 1L)
  maxPairs     <- .int(maxPairs, 1000L)
  maxMediators <- .int(maxMediators, 150L)
  use_raw      <- isTRUE(use_raw) || identical(tolower(as.character(use_raw)), "true")

  # db.path is LOCAL to each kg_* entry point (kg_assoc.R:172, kg_correlation.R:30) --
  # never a global; defining it here keeps this file independent in the same way.
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")
  if(!file.exists(db.path)) return("RES-NO-DB")

  m <- .kgLoadMetaFrame(use_raw); if(is.null(m) || !nrow(m)) return("RES-NO-META")
  # the donor row filter -- ONE implementation, `.kgDonorSubset` in kg_common.R. Applied right after
  # the frame loads, so every vector, every N and every model below describes the SAME donor set.
  # NULL means the filter could not be applied, and that REFUSES: returning the unrestricted
  # analysis as the answer to a restricted question is a wider answer wearing the asked one's face.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")
  meta_names <- names(m); rec <- as.character(m[["record_id"]])

  # ---- exposure: a real column, reduced to TWO named levels -------------------
  if(!nzchar(exposure) || !(exposure %in% meta_names)) return("RES-NO-EXPOSURE")
  lv <- trimws(strsplit(levels, "[,;]")[[1]]); lv <- lv[nzchar(lv)]
  if(length(lv) != 2) return("RES-NO-LEVELS")
  expo_raw <- as.character(m[[exposure]]); names(expo_raw) <- rec
  keep_don <- rec[!is.na(expo_raw) & expo_raw %in% lv]
  if(length(keep_don) < minN) return("RES-NO-LEVELS")
  # 0 = lv[2] (reference), 1 = lv[1]; the sign of every coefficient follows this.
  expo <- stats::setNames(as.numeric(expo_raw[keep_don] == lv[1]), keep_don)

  # ---- outcomes: parent concepts expand to their children ---------------------
  ptoks  <- trimws(strsplit(phenotypes, "[,;]")[[1]]); ptoks <- ptoks[nzchar(ptoks)]
  if(length(ptoks) == 0) return("RES-NO-PHENO")
  phenos <- .kgExpandPhenotypes(ptoks, meta_names)
  phenos <- phenos[phenos %in% meta_names]
  phenos <- phenos[vapply(phenos, function(p) identical(.kgPhenoType(m, p), "cont"), logical(1))]
  if(length(phenos) == 0) return("RES-NO-PHENO")

  # ---- mediators: named only, and capped ---------------------------------------
  mtoks <- trimws(strsplit(mediators, "[,;]")[[1]]); mtoks <- mtoks[nzchar(mtoks)]
  mtoks <- unique(mtoks[!(tolower(mtoks) %in% c("all", ""))])
  if(length(mtoks) == 0) return("RES-NO-FEATURE")
  if(length(mtoks) > maxMediators)
    return(sprintf("RES-TOO-MANY;%d;%d;%d;%d;%d",
                   length(mtoks), length(mtoks), 1L, length(phenos), as.integer(maxMediators)))

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  feats <- .kgResolveFeatures(res, mtoks, "auto")
  if(length(feats) == 0) return("RES-NO-RESOLVE")
  lay_filter <- if(identical(layers, "all")) NULL
                else tolower(trimws(strsplit(layers, "[,;]")[[1]]))

  # ---- the pair budget, CHECKED BEFORE ANY WORK -------------------------------
  # pairs = mediator x layer x outcome. Over the cap this REFUSES and reports every
  # factor, so the caller can say which one to cut rather than silently trimming.
  lay_of <- list(); n_pairs <- 0L
  for(tk in names(feats)){
    ls <- feats[[tk]]$layers
    if(!is.null(lay_filter)) ls <- Filter(function(l) tolower(l$display) %in% lay_filter, ls)
    if(length(ls)){ lay_of[[tk]] <- ls; n_pairs <- n_pairs + length(ls) * length(phenos) }
  }
  if(length(lay_of) == 0) return("RES-NO-RESOLVE")
  n_lay <- sum(vapply(lay_of, length, integer(1)))
  if(n_pairs > maxPairs)
    return(sprintf("RES-TOO-MANY;%d;%d;%d;%d;%d",
                   n_pairs, length(lay_of), n_lay, length(phenos), as.integer(maxPairs)))

  # ---- covariates: METADATA COLUMNS ONLY, as op1 -------------------------------
  covuniv <- if(covariates %in% c("none", "NA", "")) character(0)
             else if(identical(covariates, "default")) .kgMedDefaultCovs
             else trimws(strsplit(covariates, "[,;]")[[1]])
  covuniv <- unique(covuniv[nzchar(covuniv) & covuniv %in% meta_names])
  covbase <- if(length(covuniv)){ cb <- .kgCoerce(m, covuniv); rownames(cb) <- rec; cb } else NULL

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  phframe <- .kgCoerce(m, phenos); rownames(phframe) <- rec
  if(boot > 0) set.seed(seed)
  rows <- list()

  # ---- ONE READ PER LAYER, NOT ONE PER FEATURE (the N+1 already fixed elsewhere) ------
  # MEASURED 2026-08-10 before this change: 31 mediators x 137 feature-layers took 53.8s,
  # and the FITS account for 0.03s of it -- every second was a separate `SELECT *`, ~0.4s
  # each. `.kgFeatureRows` reads all of a layer's requested features in ONE IN-query and
  # returns the feature x donor matrix, so the query count falls from one-per-feature-layer
  # to one-per-layer. Identical rows, identical order; only the number of round trips moves.
  groups <- list()
  for(tk in names(lay_of)) for(la in lay_of[[tk]]){
    key <- paste(la$table, la$read_col, sep = "\r")
    if(is.null(groups[[key]]))
      groups[[key]] <- list(table = la$table, read_col = la$read_col,
                            toks = character(0), vals = character(0), disp = character(0))
    groups[[key]]$toks <- c(groups[[key]]$toks, tk)
    groups[[key]]$vals <- c(groups[[key]]$vals, la$read_val)
    groups[[key]]$disp <- c(groups[[key]]$disp, la$display)
  }

  for(g in groups){
    fr <- .kgFeatureRows(con, g$table, g$read_col, g$vals)
    if(is.null(fr)) next
    for(gi in seq_along(g$toks)){
      tk <- g$toks[gi]; la <- list(display = g$disp[gi])
      ri <- match(g$vals[gi], rownames(fr$mat))
      if(is.na(ri)) next
      med_all <- fr$mat[ri, ]; names(med_all) <- colnames(fr$mat)
      if(all(is.na(med_all))) next
      for(ph in phenos){
        # ---- THE OUTCOME IS NEVER ALSO A COVARIATE (P5, outcome-covariate leak) --
        # kg_assoc does the same with `covcols_for`: adjusting an outcome for itself
        # regresses away the very effect under test. Real here, not hypothetical --
        # .kgExpandPhenotypes("gsis") returns culturetime2 among its columns.
        covcols <- setdiff(covuniv, ph)
        cb_all  <- if(length(covcols)) covbase[, covcols, drop = FALSE] else NULL

        # ---- ONE complete-case set for all three fits --------------------------
        don <- intersect(keep_don, names(med_all))
        don <- don[!is.na(med_all[don])]
        yv  <- suppressWarnings(as.numeric(phframe[don, ph]))
        don <- don[!is.na(yv)]
        if(!is.null(cb_all)) don <- don[stats::complete.cases(cb_all[don, , drop = FALSE])]
        # the exposure must still carry BOTH levels after the drops
        if(length(don) < minN) next
        e <- expo[don]; if(length(unique(e)) < 2) next
        y   <- suppressWarnings(as.numeric(phframe[don, ph]))
        med <- as.numeric(med_all[don])

        # ---- CATEGORICAL COVARIATES MUST BE DUMMY-CODED -------------------------
        # .kgCoerce returns a FACTOR for a categorical covariate (donorsex), and
        # as.matrix() on a frame holding one yields a CHARACTER matrix that .lm.fit
        # cannot solve. A formula-based lm gets this coding for free; .lm.fit does not,
        # so model.matrix does it here -- same contrasts, intercept dropped because the
        # design already carries one. A factor left with a single level after the
        # complete-case drops is removed rather than producing a rank-deficient fit.
        C <- NULL
        if(!is.null(cb_all)){
          cdf  <- droplevels(cb_all[don, , drop = FALSE])
          keep <- vapply(cdf, function(z) !is.factor(z) || nlevels(z) > 1, logical(1))
          cdf  <- cdf[, keep, drop = FALSE]
          if(ncol(cdf)){
            C <- try(stats::model.matrix(~ ., data = cdf)[, -1, drop = FALSE], silent = TRUE)
            if(inherits(C, "try-error") || nrow(C) != length(don)) next
          }
        }
        if(!is.null(C) && any(!is.finite(C))) next

        X_tot <- if(is.null(C)) cbind(1, e) else cbind(1, e, C)
        X_dir <- cbind(X_tot, med)
        # a-path shares the TOTAL model's design; only the response changes.
        fa <- .kgMedFit(X_tot, med, 2)
        ft <- .kgMedFit(X_tot, y,   2)
        fd <- .kgMedFit(X_dir, y,   2)
        if(is.null(fa) || is.null(ft) || is.null(fd)) next
        # THE DESIGN-AWARE FLOOR. minN counts donors; this counts what is left after the
        # design spends them. Skipped, not down-weighted -- a pair that cannot support the
        # direct model cannot support a proportion derived from it.
        if(fd$df < minDF) next

        indirect <- ft$b - fd$b
        # PM IS UNDEFINED WHEN THE TOTAL EFFECT IS ~0 -- a ratio to a near-zero
        # denominator is noise, not a proportion. Reported as NA with a reason,
        # never as a large number and never silently dropped.
        pm_ok <- is.finite(ft$b) && abs(ft$b) > 1e-8 && abs(ft$b / ft$se) > 1
        pm    <- if(pm_ok) indirect / ft$b else NA_real_
        note  <- if(pm_ok) "" else "PM undefined: total effect indistinguishable from zero"

        # exposure x mediator interaction -- the difference method assumes it away,
        # so it is TESTED and reported rather than assumed.
        X_int <- cbind(X_dir, e * med)
        fi <- .kgMedFit(X_int, y, ncol(X_int))
        int_p <- if(is.null(fi)) NA_real_ else fi$p

        # ---- STAGE 2: bootstrap CI on the INDIRECT EFFECT and on PM -----------
        # !! THE INTERVAL THAT GETS PLOTTED IS THE INDIRECT ONE. PM is a RATIO, so its bootstrap
        # distribution is heavy-tailed whenever the total effect is small -- measured on the
        # published MAP1LC3A example, PM 0.92 came back as [-4.97, 7.48] while the same replicates
        # give a perfectly ordinary interval for (c - c'). The indirect effect is on the OUTCOME's
        # scale, which is what a forest plot needs and what the literature reports.
        # !! AND IT IS COMPUTED EVEN WHEN PM IS NOT. The loop was gated on `pm_ok`, i.e. skipped
        # entirely when the total effect is indistinguishable from zero -- but c - c' is perfectly
        # well defined there, and that is exactly the case where the ratio must not be quoted and
        # the absolute effect is the only honest thing to show. 15 of 40 pairs on that example.
        # The replicate ALREADY computes bt$b - bd$b to form the ratio; keeping it costs nothing.
        pm_lo <- NA_real_; pm_hi <- NA_real_; n_boot_ok <- 0L
        ind_lo <- NA_real_; ind_hi <- NA_real_
        if(boot > 0){
          n <- length(don); pms <- rep(NA_real_, boot); inds <- rep(NA_real_, boot)
          for(b in seq_len(boot)){
            i  <- sample.int(n, n, TRUE)
            bt <- .kgMedFit(X_tot[i, , drop = FALSE], y[i], 2)
            bd <- .kgMedFit(X_dir[i, , drop = FALSE], y[i], 2)
            if(is.null(bt) || is.null(bd)) next
            inds[b] <- bt$b - bd$b
            if(abs(bt$b) >= 1e-8) pms[b] <- (bt$b - bd$b) / bt$b
          }
          inds <- inds[is.finite(inds)]; n_boot_ok <- length(inds)
          if(n_boot_ok >= 100){
            qi <- stats::quantile(inds, c(0.025, 0.975), names = FALSE)
            ind_lo <- qi[1]; ind_hi <- qi[2]
          }
          pms <- pms[is.finite(pms)]
          if(pm_ok && length(pms) >= 100){
            qq <- stats::quantile(pms, c(0.025, 0.975), names = FALSE)
            pm_lo <- qq[1]; pm_hi <- qq[2]
          }
        }

        # ---- THE REMINDER, NOT A FILTER (user ruling 2026-08-10) --------------
        # A pair can clear both floors and still rest on very little. This states the
        # residual df in words so the reader is told rather than left to derive it from N.
        # It is a DISPLAY threshold: nothing is dropped, and the caller must surface it.
        warn <- if(fd$df < 30)
                  sprintf("few donors: N=%d, residual df=%d on the direct model - PM is unstable at this size",
                          length(don), fd$df)
                else ""

        rows[[length(rows) + 1L]] <- data.frame(
          exposure = exposure, level_1 = lv[1], level_0 = lv[2],
          mediator = tk, mediator_type = feats[[tk]]$type, layer = la$display,
          phenotype = ph, N = length(don),
          a_coef = fa$b, a_p = fa$p,
          b_total = ft$b, b_total_p = ft$p,
          b_direct = fd$b, b_direct_p = fd$p,
          indirect = indirect, indirect_lo = ind_lo, indirect_hi = ind_hi,
          PM = pm, PM_lo = pm_lo, PM_hi = pm_hi,
          boot_n = n_boot_ok, interaction_p = int_p, df_direct = fd$df,
          note = note, warn = warn, stringsAsFactors = FALSE)
      }
    }
  }

  if(length(rows) == 0) return("RES-NO")
  out <- do.call(rbind, rows)
  # BH across exactly the tests asked -- the a-path and the direct model each get their
  # own family, because they answer different questions (R6: never one pooled FDR).
  out$a_adj_p        <- stats::p.adjust(out$a_p,        method = "BH")
  out$b_direct_adj_p <- stats::p.adjust(out$b_direct_p, method = "BH")
  out <- out[order(-abs(ifelse(is.na(out$PM), -Inf, out$PM))), , drop = FALSE]
  # SAME OUTPUT CONVENTION AS EVERY OTHER kg_* op (kg_assoc.R:395-397): a timestamped
  # copy plus the stable name the caller reads, both in the working directory the
  # servlet sets to the user's folder -- never an absolute path.
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 30))
  utils::write.csv(out, paste0("kg_mediation_", safe(mediators), "_", safe(phenotypes),
                               "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_mediation.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgMediation <- function(...) .kgGuard("kgMediation", .kgMediationImpl, ...)
