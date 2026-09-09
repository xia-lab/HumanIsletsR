# ==============================================================================
# kg_functions/kg_multiomics.R  --  multi-omics CONCORDANCE PROXY (op#7 fast proxy).
# ------------------------------------------------------------------------------
# The immediate answer to a multi-omics question ("do the omics layers agree about
# phenotype P?"): run the omics differential for P in each requested layer, then measure
# how well the per-gene differential EFFECTS agree between every pair of layers
# (Pearson + Spearman of the effect sizes over shared genes, plus the fraction with the
# same direction). This is the fast proxy the true integration (DIABLO/MOFA/WGCNA,
# kg_integration.R) is recommended after.
#
# Self-contained: reuses .kgOmicsDE + the correlation helpers (kg_common.R). Effects are
# keyed by SYMBOL (the common identifier across gene layers).
#
# phenotype   : a metadata_sum_norm column
# contrast    : group phenotype -> "A-B" or "anova"; ""=continuous
# covariates  : "default" (the 4) | "none" | comma/semicolon list
# layers      : comma list of gene layers (rnaseq | nanostring | protein | pbrna_alpha | pbrna_beta)
# minShared   : min shared genes to report a layer pair (default 20)
# method      : "both" | "pearson" | "spearman"
# Writes kg_multiomics_concordance.csv (one row per layer pair). "RES-OK;<n_pairs>".
# ==============================================================================

.kgMoDefaultCovs <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

.kgMoLayer <- function(layer){
  switch(tolower(layer),
    "rnaseq"      = list(tbl="proc_rnaseq",           id="accession", disp="RNA-seq"),
    "nanostring"  = list(tbl="proc_nanostring_merge", id="gene_id",   disp="Nanostring"),
    "protein"     = list(tbl="proc_prot_combat",      id="id",        disp="Protein"),
    "pbrna_alpha" = list(tbl="proc_pbrna_Alpha",      id="symbol",    disp="Pseudobulk (alpha)"),
    "pbrna_beta"  = list(tbl="proc_pbrna_Beta",       id="symbol",    disp="Pseudobulk (beta)"),
    NULL)
}

kgMultiOmicsConcordance <- function(phenotype, contrast = "", ref = "",
                                    covariates = "default",
                                    layers = "rnaseq,protein,nanostring",
                                    minShared = "20", method = "both", mode = "tool",
                                    subset = "all"){   # `subset` LAST: positional callers unaffected
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minShared <- suppressWarnings(as.integer(minShared)); if(is.na(minShared)) minShared <- 20L
  method <- tolower(method)

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  # the donor row filter -- ONE implementation, `.kgDonorSubset` in kg_common.R. Applied right after
  # the frame loads, so every vector, every N and every model below describes the SAME donor set.
  # NULL means the filter could not be applied, and that REFUSES: returning the unrestricted
  # analysis as the answer to a restricted question is a wider answer wearing the asked one's face.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")
  if(!(phenotype %in% names(m))){
    j <- match(.kgNorm(phenotype), tolower(names(m))); if(is.na(j)) return("RES-NO-PHENO")
    phenotype <- names(m)[j]
  }
  ptype <- .kgInferType(m[[phenotype]]); if(is.na(ptype)) return("RES-NO-PHENO")

  if(covariates %in% c("none","NA","")){ cov <- character(0)
  } else if(identical(covariates,"default")){ cov <- .kgMoDefaultCovs[.kgMoDefaultCovs != phenotype]
  } else { cov <- trimws(strsplit(covariates, "[,;]")[[1]]); cov <- cov[cov != phenotype & nzchar(cov)] }
  cov <- cov[cov %in% names(m)]

  if(ptype == "disc"){
    if(identical(contrast,"")) return("RES-NO-CONTRAST")
    if(identical(contrast,"anova")){ contrast_use <- "anova"; ref_use <- ref
    } else { p <- trimws(strsplit(contrast,"-",fixed=TRUE)[[1]]); if(length(p)!=2) return("RES-NO-CONTRAST"); contrast_use <- p[1]; ref_use <- p[2] }
  } else { contrast_use <- ""; ref_use <- "" }

  layer.ids <- trimws(strsplit(layers, "[,;]")[[1]]); layer.ids <- layer.ids[nzchar(layer.ids)]
  lps <- lapply(layer.ids, .kgMoLayer); names(lps) <- layer.ids
  lps <- lps[!sapply(lps, is.null)]
  if(length(lps) < 2) return("RES-NO-OMICS")

  md <- .kgCoerce(m, unique(c(phenotype, cov))); rownames(md) <- as.character(m[["record_id"]])
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite")); on.exit(dbDisconnect(con), add = TRUE)

  # per-layer: symbol -> differential effect (+ symbol -> adj-p)
  eff <- list(); padj <- list(); disp <- character(0)
  for(nm in names(lps)){
    lp <- lps[[nm]]
    de <- .kgOmicsDE(con, lp$tbl, lp$id, md, phenotype, ptype, ref_use, contrast_use)
    if(is.null(de) || nrow(de) == 0) next
    de <- de[!is.na(de$Symbol) & de$Symbol != "" & !is.na(de$Effect), ]
    de <- de[order(-abs(de$Effect)), ]; de <- de[!duplicated(de$Symbol), ]     # one row per symbol
    eff[[nm]]  <- setNames(de$Effect, de$Symbol)
    padj[[nm]] <- setNames(de$Adjusted_p_value, de$Symbol)
    disp[nm]   <- lp$disp
  }
  if(length(eff) < 2) return("RES-NO")

  nm <- names(eff)
  rows <- list()
  for(i in 1:(length(nm)-1)) for(j in (i+1):length(nm)){
    a <- nm[i]; b <- nm[j]
    shared <- intersect(names(eff[[a]]), names(eff[[b]]))
    if(length(shared) < minShared) next
    ea <- eff[[a]][shared]; eb <- eff[[b]][shared]
    cp <- .kgCorPair(ea, eb, method)
    # direction agreement among genes significant in EITHER layer
    sig <- (padj[[a]][shared] < 0.05) | (padj[[b]][shared] < 0.05); sig[is.na(sig)] <- FALSE
    conc <- if(sum(sig) > 0) mean(sign(ea[sig]) == sign(eb[sig]), na.rm = TRUE) else NA_real_
    rows[[length(rows)+1]] <- data.frame(
      LayerA = disp[a], LayerB = disp[b], N_shared_genes = cp$n,
      Effect_Pearson_r = signif(cp$pear_r,3), Effect_Spearman_rho = signif(cp$spear_rho,3),
      N_sig_either = sum(sig), Concordant_direction_pct = signif(100*conc,3),
      stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  out <- out[order(-abs(ifelse(is.na(out$Effect_Spearman_rho), out$Effect_Pearson_r, out$Effect_Spearman_rho))), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_multiomics_concordance_", safe(phenotype), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_multiomics_concordance.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}
