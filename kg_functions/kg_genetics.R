# ==============================================================================
# kg_functions/kg_genetics.R  --  the GENETIC validation family (op#6): cis-QTL x cohort
# GWAS COLOCALIZATION (coloc.abf) and cis Wald-ratio MENDELIAN RANDOMIZATION, computed live
# from the local GWAS/QTL parquet warehouse. coloc + MR share the warehouse access, gene
# resolution, phenotype labels and allele harmonisation, so they live in ONE file.
#
# Self-contained: reads the parquet warehouse directly via arrow (predicate pushdown), reuses
# only .kgSetPaths()/.kgNorm() (kg_common.R). Nothing here calls another folder's function.
# The method REPRODUCES api/services/coloc_service.py exactly (train==serve): Wakefield
# approximate Bayes factor per variant -> the five coloc.abf posteriors (Giambartolomei 2014),
# priors p1=p2=1e-4/p12=1e-5, prior sd 0.15; MR = lead-cis-variant Wald ratio, alleles
# harmonised to the QTL effective allele, delta-method SE.
#
# Warehouse (other.tables.path/gwasqtl/, relative -- deploys unchanged):
#   gene_catalog.parquet          gene_symbol / ensg -> which qtl_type exist
#   gwasqtl_phenotype_labels.csv  raw_name -> display + group
#   qtl/{eQTL,pQTL,sQTL,caQTL,eQTLexon}.parquet   cis slope/slope_se per variant
#   gwas/<domain>/<phenotype>.parquet             genome-wide beta/se/p (a1 = effect allele)
#
# kgColoc(gene, qtl_type, domain, phenotype, ...) : SCAN a gene's cis-QTL vs every islet GWAS
#   (ranked by PP.H4), or TARGETED (pass domain+phenotype) to colocalize ONE pair -- what
#   validates a significant MR. Writes kg_coloc.csv. "RES-OK;<n>".
# kgMR(gene, domain, phenotype, qtl_type, ...) : two-sample cis Wald-ratio MR of the gene on
#   ONE phenotype. Writes kg_mr.csv. "RES-OK;1". Weak-instrument gate is reported, not hidden.
#
# Returns "RES-OK;<n>" | "RES-NO-WAREHOUSE" | "RES-NO-GENE" | "RES-NO-QTL"
#   | "RES-NO-CIS" | "RES-NO-GWAS" | "RES-NO-HARMONISE" | "RES-NO".
# ==============================================================================

.kgQtlTypes <- c("eQTL", "pQTL", "sQTL", "caQTL", "eQTLexon")

# warehouse dir (relative to the shared other.tables.path; never hardcoded).
.kgGwqDir <- function(){
  .kgSetPaths()
  d <- paste0(other.tables.path, "gwasqtl/")
  if(!file.exists(paste0(d, "gene_catalog.parquet"))) return(NULL)
  d
}

# gene -> list(ensg, symbol, types[]) via gene_catalog.parquet (symbol, else ENSG prefix).
.kgGwqGene <- function(dir, gene){
  gc <- tryCatch(arrow::read_parquet(paste0(dir, "gene_catalog.parquet")), error = function(e) NULL)
  if(is.null(gc)) return(NULL)
  g   <- trimws(as.character(gene))
  # which() (not logical indexing) so NA gene_symbol/ensg cells never create phantom NA rows.
  sub <- gc[which(toupper(as.character(gc$gene_symbol)) == toupper(g)), , drop = FALSE]
  if(nrow(sub) == 0 && grepl("^ENSG", toupper(g))) sub <- gc[which(as.character(gc$ensg) == g), , drop = FALSE]
  if(nrow(sub) == 0) return(NULL)
  ensg_v <- as.character(sub$ensg[!is.na(sub$ensg)])
  if(length(ensg_v) == 0) return(NULL)
  list(ensg = ensg_v[1], symbol = as.character(sub$gene_symbol[1]),
       types = unique(as.character(sub$qtl_type[!is.na(sub$qtl_type)])))
}

# pick the qtl_type: requested if the gene has it, else eQTL, else the first available.
.kgGwqPickQtl <- function(qtl_type, types){
  if(qtl_type %in% types) return(qtl_type)
  if("eQTL" %in% types) return("eQTL")
  sort(types)[1]
}

# phenotype labels (raw_name -> c(display, group)), cached by mtime.
.kgGwqLabels <- function(dir){
  f <- paste0(dir, "gwasqtl_phenotype_labels.csv")
  if(!file.exists(f)) return(list())
  mt <- file.info(f)$mtime
  if(is.null(.GlobalEnv$.kg_gwqlab_cache) || !identical(.GlobalEnv$.kg_gwqlab_mtime, mt)){
    lab <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    L <- list()
    for(i in seq_len(nrow(lab))){
      disp <- if(!is.null(lab$display) && nzchar(as.character(lab$display[i]))) as.character(lab$display[i]) else as.character(lab$raw_name[i])
      grp  <- if(!is.null(lab$group)) as.character(lab$group[i]) else ""
      L[[as.character(lab$raw_name[i])]] <- c(disp, grp)
    }
    .GlobalEnv$.kg_gwqlab_cache <- L; .GlobalEnv$.kg_gwqlab_mtime <- mt
  }
  .GlobalEnv$.kg_gwqlab_cache
}

# a gwas file stem is "<model>__<raw_name>"; label is keyed by raw_name. -> c(display,group,model)
.kgGwqLabel <- function(labels, raw){
  if(grepl("__", raw, fixed = TRUE)){ p <- strsplit(raw, "__", fixed = TRUE)[[1]]; model <- p[1]; rest <- paste(p[-1], collapse = "__") }
  else { model <- ""; rest <- raw }
  lv <- labels[[rest]]
  c(if(!is.null(lv)) lv[1] else rest, if(!is.null(lv)) lv[2] else "", model)
}

# ---- coloc.abf math (reproduces coloc_service.py exactly) ---------------------
.kgLabf   <- function(beta, se, W){ v <- se * se; r <- W / (v + W); 0.5 * log(1 - r) + 0.5 * r * (beta * beta / v) }
.kgLogsum <- function(x){ m <- max(x); m + log(sum(exp(x - m))) }
.kgColocAbf <- function(b1, s1, b2, s2, W = 0.15^2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5){
  l1 <- .kgLabf(b1, s1, W); l2 <- .kgLabf(b2, s2, W)
  ls1 <- .kgLogsum(l1); ls2 <- .kgLogsum(l2); ls12 <- .kgLogsum(l1 + l2)
  lH1 <- log(p1) + ls1; lH2 <- log(p2) + ls2
  base <- ls1 + ls2; diff_arg <- ls12 - base
  lH3 <- if(diff_arg >= 0) -Inf else log(p1) + log(p2) + base + log1p(-exp(diff_arg))
  lH4 <- log(p12) + ls12
  arr <- c(0, lH1, lH2, lH3, lH4); fin <- is.finite(arr)
  d <- .kgLogsum(arr[fin]); pp <- ifelse(fin, exp(arr - d), 0)
  list(pp = pp, shared = l1 + l2)          # pp = c(H0,H1,H2,H3,H4)
}

# ---- warehouse readers (arrow predicate pushdown) ----------------------------
# a gene's cis QTL (superset of columns), cleaned: drop NA slope/se/rsid, keep slope_se>0.
.kgGwqCis <- function(dir, qtl_type, ensg){
  qf <- paste0(dir, "qtl/", qtl_type, ".parquet"); if(!file.exists(qf)) return(NULL)
  q <- tryCatch(
    dplyr::collect(dplyr::select(dplyr::filter(arrow::open_dataset(qf), startsWith(phe_id, ensg)),
      rsid, chr, pos, slope, slope_se, nom_pval, effective_allele, non_effective_allele)),
    error = function(e) NULL)
  if(is.null(q) || nrow(q) == 0) return(NULL)
  q <- q[!is.na(q$slope) & !is.na(q$slope_se) & !is.na(q$rsid) & q$slope_se > 0, , drop = FALSE]
  if(nrow(q) == 0) return(NULL)
  q
}

# a gwas cis window (chr==chrom & lo<=pos<=hi), columns as requested, cleaned se>0.
.kgGwqGwasWindow <- function(f, chrom, lo, hi, want_a1){
  d <- tryCatch(
    if(want_a1) dplyr::collect(dplyr::select(dplyr::filter(arrow::open_dataset(f), chr == chrom & pos >= lo & pos <= hi), rsid, a1, beta, se, p))
    else        dplyr::collect(dplyr::select(dplyr::filter(arrow::open_dataset(f), chr == chrom & pos >= lo & pos <= hi), rsid, beta, se, p)),
    error = function(e) NULL)
  if(is.null(d) || nrow(d) == 0) return(NULL)
  d <- d[!is.na(d$beta) & !is.na(d$se) & !is.na(d$rsid) & d$se > 0, , drop = FALSE]
  if(nrow(d) == 0) return(NULL)
  d
}

# ==============================================================================
# kgColoc: SCAN a gene's cis-QTL vs cohort GWAS (rank by PP.H4), or TARGETED one pair.
# ==============================================================================
kgColoc <- function(gene, qtl_type = "eQTL", domain = "", phenotype = "",
                    min_snps = "50", top = "12", strong_p = "5e-6", mode = "tool"){
  library(arrow); library(dplyr)
  dir <- .kgGwqDir(); if(is.null(dir)) return("RES-NO-WAREHOUSE")
  min_snps <- suppressWarnings(as.integer(min_snps)); if(is.na(min_snps)) min_snps <- 50L
  topn     <- suppressWarnings(as.integer(top));      if(is.na(topn)) topn <- 12L
  sp       <- suppressWarnings(as.numeric(strong_p)); if(is.na(sp)) sp <- 5e-6

  info <- .kgGwqGene(dir, gene); if(is.null(info)) return("RES-NO-GENE")
  qt <- .kgGwqPickQtl(qtl_type, info$types)
  q  <- .kgGwqCis(dir, qt, info$ensg); if(is.null(q)) return("RES-NO-QTL")
  if(nrow(q) < min_snps) return("RES-NO-CIS")

  chrom <- as.character(q$chr[1]); lo <- min(q$pos); hi <- max(q$pos)
  qlead <- q[which.min(q$nom_pval), ]
  qsm   <- q[, c("rsid", "slope", "slope_se", "nom_pval")]
  labels <- .kgGwqLabels(dir)

  targeted <- nzchar(domain) && nzchar(phenotype)
  if(targeted){
    tf <- paste0(dir, "gwas/", domain, "/", phenotype, ".parquet")
    if(!file.exists(tf)) return("RES-NO-GWAS")
    files <- tf
  } else {
    files <- list.files(paste0(dir, "gwas"), pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  }
  if(length(files) == 0) return("RES-NO-GWAS")

  rows <- list()
  for(f in files){
    dom <- basename(dirname(f)); raw <- sub("\\.parquet$", "", basename(f))
    g <- .kgGwqGwasWindow(f, chrom, lo, hi, want_a1 = FALSE); if(is.null(g)) next
    mm <- merge(qsm, g, by = "rsid"); if(nrow(mm) < min_snps) next
    ca <- .kgColocAbf(mm$slope, mm$slope_se, mm$beta, mm$se)
    bi <- which.max(ca$shared)
    lab <- .kgGwqLabel(labels, raw)
    rows[[length(rows) + 1]] <- data.frame(
      Gene = info$symbol, Ensembl = info$ensg, QTL_type = qt, Region = paste0(chrom, ":", lo, "-", hi),
      Domain = dom, Phenotype = raw, Display = lab[1], Group = if(nzchar(lab[2])) lab[2] else dom, Model = lab[3],
      PP_H4 = round(ca$pp[5], 3), PP_H3 = round(ca$pp[4], 3), N_snps = nrow(mm),
      Lead_variant = as.character(mm$rsid[bi]), QTL_lead_p = signif(mm$nom_pval[bi], 4), GWAS_lead_p = signif(mm$p[bi], 4),
      N_cis = nrow(q), Instrument_strong = as.character(qlead$nom_pval < sp), stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")
  out <- do.call(rbind, rows); out <- out[order(-out$PP_H4), , drop = FALSE]
  if(nrow(out) > topn) out <- out[seq_len(topn), , drop = FALSE]

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 24))
  utils::write.csv(out, paste0("kg_coloc_", safe(info$symbol), "_", safe(if(targeted) phenotype else "scan"), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_coloc.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ==============================================================================
# kgMR: two-sample cis Wald-ratio MR of a gene on ONE phenotype (lead cis variant).
# ==============================================================================
kgMR <- function(gene, domain, phenotype, qtl_type = "eQTL", strong_p = "5e-6", mode = "tool"){
  library(arrow); library(dplyr)
  dir <- .kgGwqDir(); if(is.null(dir)) return("RES-NO-WAREHOUSE")
  sp <- suppressWarnings(as.numeric(strong_p)); if(is.na(sp)) sp <- 5e-6

  info <- .kgGwqGene(dir, gene); if(is.null(info)) return("RES-NO-GENE")
  qt <- .kgGwqPickQtl(qtl_type, info$types)
  gf <- paste0(dir, "gwas/", domain, "/", phenotype, ".parquet"); if(!file.exists(gf)) return("RES-NO-GWAS")
  q  <- .kgGwqCis(dir, qt, info$ensg); if(is.null(q)) return("RES-NO-QTL")

  lead   <- q[which.min(q$nom_pval), ]
  strong <- as.numeric(lead$nom_pval) < sp
  g <- tryCatch(
    dplyr::collect(dplyr::select(dplyr::filter(arrow::open_dataset(gf), rsid == as.character(lead$rsid)), rsid, a1, beta, se, p)),
    error = function(e) NULL)
  if(is.null(g) || nrow(g) == 0) return("RES-NO-GWAS")
  gr <- g[1, ]

  eff <- toupper(as.character(lead$effective_allele)); nef <- toupper(as.character(lead$non_effective_allele))
  a1  <- toupper(as.character(gr$a1))
  beta_out <- if(a1 == eff) as.numeric(gr$beta) else if(a1 == nef) -as.numeric(gr$beta) else return("RES-NO-HARMONISE")

  beta_exp <- as.numeric(lead$slope); se_exp <- as.numeric(lead$slope_se); se_out <- as.numeric(gr$se)
  theta    <- beta_out / beta_exp
  se_theta <- abs(theta) * sqrt((se_out / beta_out)^2 + (se_exp / beta_exp)^2)
  z <- if(se_theta > 0) theta / se_theta else 0
  p <- 2 * stats::pnorm(-abs(z))                       # == erfc(|z|/sqrt2), two-sided
  labels <- .kgGwqLabels(dir); lab <- .kgGwqLabel(labels, phenotype)

  out <- data.frame(
    Gene = info$symbol, Ensembl = info$ensg, QTL_type = qt, Domain = domain, Phenotype = phenotype,
    Display = lab[1], Group = if(nzchar(lab[2])) lab[2] else domain, Model = lab[3],
    Instrument = as.character(lead$rsid), Instrument_p = signif(as.numeric(lead$nom_pval), 4),
    Instrument_strong = as.character(strong), Beta = round(theta, 4), SE = round(se_theta, 4), P = signif(p, 4),
    Note = if(strong) "" else "instrument is weak (cis-QTL not genome-wide); MR estimate is suggestive only",
    stringsAsFactors = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 24))
  utils::write.csv(out, paste0("kg_mr_", safe(info$symbol), "_", safe(phenotype), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_mr.csv", row.names = FALSE)
  "RES-OK;1"
}

# ==============================================================================
# kgMRTwoSample: TWO-SAMPLE cis Wald-ratio MR -- cohort cis-QTL (exposure) vs a PUBLIC
# hugeamp (KP) trait (outcome). Reproduces mr_service.py / GENETIC_DATA_AND_ANALYSIS_
# REFERENCE.md section 3.2 EXACTLY, but FILE-BACKED (no neo4j):
#   exposure = the q<0.05 cis lead from display_kg/kg_qtl_leads.csv (built from the same
#              local_genetic_leads JSON neo4j's LOCAL_EQTL/LOCAL_PQTL were injected from;
#              q is a global BH-FDR, so it is read precomputed, never recomputed here).
#   outcome  = sqlite/HI_gene_gwas.sqlite gene_trait_variants, source='hugeamp' (the KP
#              layer; gwas_catalog has NO std_err and is excluded). effect_allele = alt,
#              other_allele = ref (KP beta-per-ALT; the exact rule patch_variant_assoc_
#              alleles.py applied to neo4j).
# MR exposures are eQTL / pQTL ONLY (exon/sQTL/caQTL deferred). Wald ratio with the
# first-order (outcome-only) SE the reference specifies -- deliberately different from the
# one-sample kgMR delta-method SE. p via 2*pnorm(-|z|) (== erfc(|z|/sqrt2), tail-stable).
# BH-FDR across exactly the traits this call returns; palindromic variants flagged
# sign-uncertain; a disease NAME expands to its hugeamp trait family.
#
# kgMRTwoSample(gene, trait, qtl_type, max_pairs, mode). trait "" = every trait the
# instrument has a hugeamp beta for. Writes kg_mr.csv. "RES-OK;<n>" | RES-NO-QTL
#   | RES-NO-LEADS | RES-NO-INSTRUMENT | RES-NO-DB | RES-NO-OUTCOME | RES-NO.
# ==============================================================================

# cis-QTL lead table (display_kg/kg_qtl_leads.csv), cached by mtime.
.kgQtlLeads <- function(){
  .kgSetPaths()
  f <- paste0(other.tables.path, "display_kg/kg_qtl_leads.csv")
  if(!file.exists(f)) return(NULL)
  mt <- file.info(f)$mtime
  if(is.null(.GlobalEnv$.kg_leads_cache) || !identical(.GlobalEnv$.kg_leads_mtime, mt)){
    .GlobalEnv$.kg_leads_cache <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    .GlobalEnv$.kg_leads_mtime <- mt
  }
  .GlobalEnv$.kg_leads_cache
}

.kgPalindromic <- function(a, b){
  paste0(sort(toupper(c(a, b))), collapse = "") %in% c("AT", "CG")   # A/T or C/G
}

# normalise a trait string for matching: lowercase, drop everything but [a-z0-9]
# ("Type 2 diabetes" == "type 2 diabetes" == "T2D "-> no; the NAME and the CODE are distinct
# keys, both indexed below -- this only unifies spacing/punctuation/case within one key).
.kgNormTrait <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))

# DATA-DRIVEN trait name->code map from HI_gene_gwas.sqlite `traits` (149 rows; id=code,
# name=human, hugeamp_phenotype=KP code). Every code present in gene_trait_variants is a
# traits.id, so this table IS the authority. Each of id / name / hugeamp_phenotype is indexed
# (normalised) -> the canonical id. Cached by the sqlite mtime. Replaces the old 3-entry
# hardcoded T2D substring map; resolution is now EXACT (a name -> its one code), so a T2D
# query no longer bleeds into complication traits ("CADinT2D") the old CONTAINS matched.
.kgTraitMap <- function(con, db){
  mt <- file.info(db)$mtime
  if(is.null(.GlobalEnv$.kg_traitmap) || !identical(.GlobalEnv$.kg_traitmap_mtime, mt)){
    tr <- try(DBI::dbGetQuery(con, "SELECT id, name, hugeamp_phenotype FROM traits"), silent = TRUE)
    m <- list()
    if(!inherits(tr, "try-error") && nrow(tr) > 0){
      for(i in seq_len(nrow(tr))){
        id <- as.character(tr$id[i]); if(is.na(id) || !nzchar(id)) next
        for(k in c(id, tr$name[i], tr$hugeamp_phenotype[i])){
          nk <- .kgNormTrait(k); if(nzchar(nk) && is.null(m[[nk]])) m[[nk]] <- id
        }
      }
    }
    .GlobalEnv$.kg_traitmap <- m; .GlobalEnv$.kg_traitmap_mtime <- mt
  }
  .GlobalEnv$.kg_traitmap
}

# resolve a user trait string -> the canonical trait CODE (a name, a code, or a KP phenotype
# all map to the code); an unknown string is returned as-is so it exact-matches to nothing.
.kgResolveTrait <- function(con, db, trait){
  m <- .kgTraitMap(con, db); code <- m[[.kgNormTrait(trait)]]
  if(!is.null(code)) code else trimws(trait)
}

kgMRTwoSample <- function(gene, trait = "", qtl_type = "eQTL", max_pairs = "200", mode = "tool"){
  library(RSQLite); library(DBI)
  .kgSetPaths()
  if(!(qtl_type %in% c("eQTL", "pQTL"))) return("RES-NO-QTL")   # MR exposures eQTL/pQTL only
  mp <- suppressWarnings(as.integer(max_pairs)); if(is.na(mp)) mp <- 200L

  leads <- .kgQtlLeads(); if(is.null(leads)) return("RES-NO-LEADS")
  sub <- leads[leads$qtl_type == qtl_type & toupper(leads$symbol) == toupper(gene), , drop = FALSE]
  if(nrow(sub) == 0 && grepl("^ENSG", toupper(gene)))
    sub <- leads[leads$qtl_type == qtl_type & toupper(leads$ensg) == toupper(gene), , drop = FALSE]
  sub <- sub[!is.na(sub$q) & sub$q < 0.05, , drop = FALSE]      # q<0.05 gate (= presence in neo4j)
  if(nrow(sub) == 0) return("RES-NO-INSTRUMENT")
  lead <- sub[which.min(sub$nom_p), ]                          # one lead per gene per modality

  db <- paste0(sqlite.path, "HI_gene_gwas.sqlite"); if(!file.exists(db)) return("RES-NO-DB")
  con <- dbConnect(SQLite(), db); on.exit(dbDisconnect(con), add = TRUE)

  where <- "source='hugeamp' AND rsid=? AND beta IS NOT NULL AND std_err IS NOT NULL"
  params <- list(as.character(lead$rsid))
  if(nzchar(trait)){                              # resolve name/code/KP-phenotype -> exact code
    where <- paste(where, "AND trait = ?"); params <- c(params, .kgResolveTrait(con, db, trait))
  }
  o <- try(dbGetQuery(con, sprintf("SELECT trait, ref, alt, beta, std_err, pvalue, n FROM gene_trait_variants WHERE %s", where), params = params), silent = TRUE)
  if(inherits(o, "try-error") || nrow(o) == 0) return("RES-NO-OUTCOME")
  o <- o[!duplicated(o$trait), , drop = FALSE]                 # one row per (variant,trait)
  if(nrow(o) > mp) o <- o[seq_len(mp), , drop = FALSE]

  exp_eff <- toupper(lead$effect_allele); exp_oth <- toupper(lead$other_allele); slope <- as.numeric(lead$slope)
  pal <- .kgPalindromic(exp_eff, exp_oth)
  rows <- list()
  for(i in seq_len(nrow(o))){
    out_eff <- toupper(o$alt[i]); out_oth <- toupper(o$ref[i])  # KP rule: effect=alt, other=ref
    if(exp_eff == out_eff){ slope_al <- slope;  harm <- "matched" }
    else if(exp_eff == out_oth){ slope_al <- -slope; harm <- "flipped" }
    else next
    if(slope_al == 0) next
    beta_out <- as.numeric(o$beta[i]); se_out <- as.numeric(o$std_err[i])
    theta <- beta_out / slope_al; se_theta <- se_out / abs(slope_al)
    if(se_theta <= 0) next
    z <- theta / se_theta; p <- 2 * stats::pnorm(-abs(z))       # == erfc(|z|/sqrt2)
    rows[[length(rows)+1]] <- data.frame(
      Gene = as.character(lead$symbol), QTL_type = qtl_type, Trait = o$trait[i], Design = "two-sample",
      Instrument = as.character(lead$rsid), Instrument_p = signif(as.numeric(lead$nom_p), 4),
      Instrument_q = signif(as.numeric(lead$q), 4),
      Beta = round(theta, 6), SE = round(se_theta, 6),
      CI_low = round(theta - 1.96 * se_theta, 6), CI_high = round(theta + 1.96 * se_theta, 6),
      P = signif(p, 4), Z = round(z, 4), N_snps = 1L, Harmonization = harm,
      Palindromic = as.character(pal), Outcome_n = o$n[i], Outcome_p = signif(as.numeric(o$pvalue[i]), 4),
      Note = if(pal) "palindromic variant (A/T or C/G): SIGN of the estimate is uncertain" else "",
      stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")
  out <- do.call(rbind, rows)
  out$P_adj <- signif(stats::p.adjust(out$P, method = "BH"), 4)   # FDR across THIS call's traits
  out <- out[order(out$P), , drop = FALSE]

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 24))
  utils::write.csv(out, paste0("kg_mr_2s_", safe(as.character(lead$symbol)), "_", safe(if(nzchar(trait)) trait else "all"), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_mr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}
