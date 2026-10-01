# ==============================================================================
# kg_functions/kg_integration.R  --  the true MULTI-OMICS INTEGRATION (op#7, "recommend":
# heavier than the concordance proxy). Three methods, self-contained via the underlying
# CRAN/Bioc packages (never web_tool_functions):
#   kgDiablo  supervised   (mixOmics::block.splsda) -- multi-omics signature separating a group
#   kgWgcna   co-expression (WGCNA)                 -- modules within a layer + module-trait corr
#   kgMofa    unsupervised (MOFA2)                  -- latent factors + factor-phenotype assoc
#
# Shared helpers (below) load an omics layer as a donors x features matrix, top-variance
# filtered + mean-imputed, keyed by symbol; and align donors across blocks.
# ==============================================================================

.kgIntLayer <- function(layer){
  switch(tolower(layer),
    "rnaseq"      = list(tbl="proc_rnaseq",           id="accession", disp="RNA-seq"),
    "nanostring"  = list(tbl="proc_nanostring_merge", id="gene_id",   disp="Nanostring"),
    "protein"     = list(tbl="proc_prot_combat",      id="id",        disp="Protein"),
    "pbrna_alpha" = list(tbl="proc_pbrna_Alpha",      id="symbol",    disp="Pseudobulk (alpha)"),
    "pbrna_beta"  = list(tbl="proc_pbrna_Beta",       id="symbol",    disp="Pseudobulk (beta)"),
    "metabolite"  = list(tbl="proc_metabolite_combat_HG", id="compound", disp="Metabolite (HG)"),
    NULL)
}

# load a layer as donors x features (rownames = R### donors, colnames = symbol/compound),
# top-`topVar`-variance features, mean-imputed. Returns a matrix or NULL.
.kgIntBlock <- function(con, lp, topVar){
  ft <- try(DBI::dbReadTable(con, lp$tbl), silent = TRUE); if(inherits(ft, "try-error")) return(NULL)
  key.col <- if("symbol" %in% names(ft) && lp$id != "compound") "symbol" else lp$id
  keys <- as.character(ft[[key.col]])
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  ft <- ft[keep, ]; keys <- keys[keep]
  donor_cols <- grep("^R[0-9]+$", names(ft), value = TRUE)
  if(length(donor_cols) == 0) return(NULL)
  mat <- suppressWarnings(as.matrix(sapply(ft[, donor_cols, drop = FALSE], function(z) as.numeric(as.character(z)))))
  rownames(mat) <- keys                                   # features x donors
  # top-variance features
  v <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat <- mat[order(-v)[seq_len(min(topVar, sum(!is.na(v))))], , drop = FALSE]
  # mean-impute remaining NA per feature
  for(i in seq_len(nrow(mat))){ na <- is.na(mat[i, ]); if(any(na)) mat[i, na] <- mean(mat[i, ], na.rm = TRUE) }
  t(mat)                                                  # donors x features
}

# ==============================================================================
# kgDiablo -- supervised multi-omics integration (mixOmics DIABLO / block.splsda).
# Finds a sparse multi-omics signature that separates the groups AND correlates across
# blocks. phenotype must be categorical; `contrast`="A-B" restricts to two groups.
# layers (>=2), ncomp, keepX (selected features per block per component), topVar
# (variance pre-filter per block). Writes kg_diablo.csv (Block/Component/Feature/Loading).
# ==============================================================================
kgDiablo <- function(phenotype, contrast = "", layers = "rnaseq,protein",
                     ncomp = "2", keepX = "20", topVar = "2000", mode = "tool"){
  suppressMessages(library(mixOmics)); library(RSQLite); library(DBI)
  .kgSetPaths()
  ncomp  <- max(1L, suppressWarnings(as.integer(ncomp)));  if(is.na(ncomp))  ncomp  <- 2L
  keepXn <- max(2L, suppressWarnings(as.integer(keepX)));  if(is.na(keepXn)) keepXn <- 20L
  topVar <- max(50L, suppressWarnings(as.integer(topVar))); if(is.na(topVar)) topVar <- 2000L

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  if(!(phenotype %in% names(m))){
    j <- match(.kgNorm(phenotype), tolower(names(m))); if(is.na(j)) return("RES-NO-PHENO"); phenotype <- names(m)[j]
  }
  if(.kgInferType(m[[phenotype]]) != "disc") return("RES-NO-DISC")   # DIABLO is supervised

  lids <- trimws(strsplit(layers, "[,;]")[[1]]); lids <- lids[nzchar(lids)]
  lps  <- lapply(lids, .kgIntLayer); names(lps) <- lids; lps <- lps[!sapply(lps, is.null)]
  if(length(lps) < 2) return("RES-NO-OMICS")

  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite")); on.exit(dbDisconnect(con), add = TRUE)
  blocks <- list(); disp <- character(0)
  for(nm in names(lps)){ b <- .kgIntBlock(con, lps[[nm]], topVar); if(!is.null(b)){ blocks[[nm]] <- b; disp[nm] <- lps[[nm]]$disp } }
  if(length(blocks) < 2) return("RES-NO")

  # outcome + donor alignment
  yv <- setNames(as.character(m[[phenotype]]), as.character(m[["record_id"]]))
  grp <- if(nzchar(contrast)){ p <- trimws(strsplit(contrast,"-",fixed=TRUE)[[1]]); if(length(p)==2) p else NULL } else NULL
  keep.donor <- Reduce(intersect, lapply(blocks, rownames))
  keep.donor <- keep.donor[!is.na(yv[keep.donor]) & yv[keep.donor] != "" & yv[keep.donor] != "NA"]
  if(!is.null(grp)) keep.donor <- keep.donor[yv[keep.donor] %in% grp]
  if(length(keep.donor) < 20) return("RES-NO")

  X <- lapply(blocks, function(b) b[keep.donor, , drop = FALSE])
  Y <- factor(yv[keep.donor])
  if(length(levels(Y)) < 2) return("RES-NO")
  ncomp <- min(ncomp, length(levels(Y)) - 1L + 1L)   # splsda supports up to (#classes-1)+ comps; keep modest

  design <- matrix(0.1, nrow = length(X), ncol = length(X)); diag(design) <- 0
  keepX.list <- lapply(X, function(z) rep(min(keepXn, ncol(z)), ncomp)); names(keepX.list) <- names(X)

  fit <- try(mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX.list, design = design),
             silent = TRUE)
  if(inherits(fit, "try-error")) return(paste0("RES-NO-FIT;", conditionMessage(attr(fit, "condition"))))

  rows <- list()
  for(nm in names(X)){
    for(cc in seq_len(ncomp)){
      sv <- try(mixOmics::selectVar(fit, block = nm, comp = cc), silent = TRUE)
      if(inherits(sv, "try-error") || is.null(sv[[nm]])) next
      nmv <- sv[[nm]]$name; ld <- sv[[nm]]$value$value.var
      if(length(nmv) == 0) next
      rows[[length(rows)+1]] <- data.frame(
        Block = unname(disp[nm]), Component = cc, Feature = as.character(nmv),
        Loading = signif(unname(ld), 3), stringsAsFactors = FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows); out <- out[order(out$Component, -abs(out$Loading)), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_diablo_", safe(phenotype), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_diablo.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out), ";n_donors=", length(keep.donor))
}


# ==============================================================================
# kgWgcna -- co-expression network (WGCNA) within one omics layer + module-trait
# correlation. Builds signed co-expression modules over the top-variance features,
# then correlates each module's eigengene with the requested donor traits.
# omics_layer (one), traits (comma list of numeric metadata columns), topVar, power
# (soft-threshold), minModuleSize. Writes kg_wgcna.csv (Module/N_genes/Trait/Correlation/P_value).
# ==============================================================================
kgWgcna <- function(omics_layer = "rnaseq", traits = "hba1c,donorage,bodymassindex",
                    topVar = "3000", power = "6", minModuleSize = "30", mode = "tool"){
  suppressMessages(library(WGCNA)); library(RSQLite); library(DBI)
  .kgSetPaths()
  topVar <- max(100L, suppressWarnings(as.integer(topVar))); if(is.na(topVar)) topVar <- 3000L
  power  <- suppressWarnings(as.integer(power)); if(is.na(power)) power <- 6L
  minMod <- max(5L, suppressWarnings(as.integer(minModuleSize))); if(is.na(minMod)) minMod <- 30L

  lp <- .kgIntLayer(omics_layer); if(is.null(lp)) return("RES-NO-OMICS")
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  datExpr <- .kgIntBlock(con, lp, topVar); dbDisconnect(con)
  if(is.null(datExpr) || nrow(datExpr) < 15) return("RES-NO")

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  tnames <- trimws(strsplit(traits, "[,;]")[[1]]); tnames <- tnames[tnames %in% names(m)]
  if(length(tnames) == 0) return("RES-NO-TRAIT")
  tmat <- .kgCoerce(m, tnames); rownames(tmat) <- as.character(m[["record_id"]])
  tmat <- data.frame(lapply(tmat, function(z) suppressWarnings(as.numeric(as.character(z)))), row.names = rownames(tmat))
  names(tmat) <- tnames

  common <- intersect(rownames(datExpr), rownames(tmat))
  if(length(common) < 15) return("RES-NO")
  datExpr <- datExpr[common, , drop = FALSE]; tmat <- tmat[common, , drop = FALSE]

  net <- try(WGCNA::blockwiseModules(datExpr, power = power, networkType = "signed",
              minModuleSize = minMod, numericLabels = TRUE, maxBlockSize = ncol(datExpr),
              saveTOMs = FALSE, verbose = 0), silent = TRUE)
  if(inherits(net, "try-error")) return(paste0("RES-NO-FIT;", conditionMessage(attr(net, "condition"))))

  MEs <- net$MEs                                    # donors x modules (ME0 = unassigned)
  MEs <- MEs[, colnames(MEs) != "ME0", drop = FALSE]
  if(ncol(MEs) == 0) return("RES-NO")
  nGenes <- table(net$colors)

  rows <- list()
  for(mod in colnames(MEs)){
    modnum <- sub("^ME", "", mod)
    ng <- as.integer(nGenes[modnum]); if(is.na(ng)) ng <- NA_integer_
    for(tr in tnames){
      ok <- is.finite(MEs[[mod]]) & is.finite(tmat[[tr]])
      if(sum(ok) < 10) next
      r  <- stats::cor(MEs[ok, mod], tmat[ok, tr])
      p  <- WGCNA::corPvalueStudent(r, sum(ok))
      rows[[length(rows)+1]] <- data.frame(
        Module = paste0("M", modnum), N_genes = ng, Trait = tr,
        Correlation = signif(r, 3), P_value = signif(p, 3), N = sum(ok),
        stringsAsFactors = FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows); out <- out[order(out$P_value), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_wgcna_", safe(omics_layer), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_wgcna.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out), ";n_modules=", ncol(MEs))
}


# ==============================================================================
# kgMofa -- unsupervised multi-omics integration (MOFA2). Learns latent factors shared
# across the omics views, then correlates each factor with the requested phenotypes.
# layers (>=2), num_factors, phenotypes (numeric metadata columns), topVar. Writes
# kg_mofa.csv (Factor/Phenotype/Correlation/P_value). NOTE: MOFA2 trains via a Python
# backend (basilisk); the first run may install it.
# ==============================================================================
kgMofa <- function(layers = "rnaseq,protein", num_factors = "5",
                   phenotypes = "hba1c,donorage,bodymassindex", topVar = "2000", mode = "tool"){
  suppressMessages(library(MOFA2)); library(RSQLite); library(DBI)
  .kgSetPaths()
  nfac   <- max(2L, suppressWarnings(as.integer(num_factors))); if(is.na(nfac)) nfac <- 5L
  topVar <- max(100L, suppressWarnings(as.integer(topVar)));    if(is.na(topVar)) topVar <- 2000L

  lids <- trimws(strsplit(layers, "[,;]")[[1]]); lids <- lids[nzchar(lids)]
  lps  <- lapply(lids, .kgIntLayer); names(lps) <- lids; lps <- lps[!sapply(lps, is.null)]
  if(length(lps) < 2) return("RES-NO-OMICS")

  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  views <- list()
  for(nm in names(lps)){ b <- .kgIntBlock(con, lps[[nm]], topVar); if(!is.null(b)) views[[lps[[nm]]$disp]] <- t(b) }  # features x donors
  dbDisconnect(con)
  if(length(views) < 2) return("RES-NO")

  common <- Reduce(intersect, lapply(views, colnames))
  if(length(common) < 15) return("RES-NO")
  views <- lapply(views, function(v) v[, common, drop = FALSE])

  obj <- try(MOFA2::create_mofa(views), silent = TRUE)
  if(inherits(obj, "try-error")) return(paste0("RES-NO-MOFA;", conditionMessage(attr(obj, "condition"))))
  dopts <- MOFA2::get_default_data_options(obj)
  mopts <- MOFA2::get_default_model_options(obj); mopts$num_factors <- nfac
  topts <- MOFA2::get_default_training_options(obj); topts$verbose <- FALSE; topts$seed <- 1L
  obj <- MOFA2::prepare_mofa(obj, data_options = dopts, model_options = mopts, training_options = topts)

  outf <- tempfile(fileext = ".hdf5")
  trained <- try(suppressMessages(MOFA2::run_mofa(obj, outfile = outf, use_basilisk = TRUE, save_data = FALSE)), silent = TRUE)
  if(inherits(trained, "try-error")) return(paste0("RES-NO-MOFA-BACKEND;", conditionMessage(attr(trained, "condition"))))

  Z <- MOFA2::get_factors(trained, factors = "all")[[1]]   # donors x factors
  if(is.null(Z) || nrow(Z) == 0) return("RES-NO")

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  pnames <- trimws(strsplit(phenotypes, "[,;]")[[1]]); pnames <- pnames[pnames %in% names(m)]
  if(length(pnames) == 0) return("RES-NO-TRAIT")
  pmat <- data.frame(lapply(pnames, function(c) suppressWarnings(as.numeric(m[[c]]))))
  names(pmat) <- pnames; rownames(pmat) <- as.character(m[["record_id"]])

  rows <- list()
  for(f in colnames(Z)){
    for(p in pnames){
      d <- intersect(rownames(Z), rownames(pmat))
      x <- Z[d, f]; y <- pmat[d, p]; ok <- is.finite(x) & is.finite(y)
      if(sum(ok) < 10) next
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok]))
      rows[[length(rows)+1]] <- data.frame(Factor = f, Phenotype = p,
        Correlation = signif(unname(ct$estimate), 3), P_value = signif(ct$p.value, 3), N = sum(ok),
        stringsAsFactors = FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows); out <- out[order(out$P_value), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_mofa_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_mofa.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out), ";n_factors=", ncol(Z), ";n_donors=", nrow(Z))
}
