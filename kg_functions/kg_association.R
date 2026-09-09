# ==============================================================================
# kg_functions/kg_association.R  --  op#1 ASSOCIATION (feature <-> phenotype).
# Frames F1 (gene->phenotype), F20 (continuous), F7 (group), F30 (custom covariates),
# F26 (confounder test).
#
# Self-contained: runs the omics differential itself via .kgOmicsDE (kg_common.R) --
# limma over all features of the gene's layer, covariate-adjusted, then reports just the
# queried gene's row per layer. Nothing here calls another folder's function.
#
# DEFAULT covariates = age, sex, BMI, predistribution culture time (the primary variable is
# dropped if it is one of them).
# ⚠ PURITY REMOVED 2026-09-04, to match the precompute default, which dropped it 2026-08
# ("purity dropped 2026-08 (4 covariates for large omics)",
# local_scripts/5_precompute_associations/1_omics_meta_associations_v2_parallel.R:168).
# The on-demand default MUST name the same set as the precompute: the route decision below
# retrieves the precompute for a plain default spec and computes here otherwise, so if the
# two lists differ the same question returns a differently-adjusted answer depending on
# which path served it.
#
# Route: RETRIEVE the precompute for the plain default spec; AUTO here when the spec
# differs -- a donor SUBSET, CUSTOM covariates (F30), a confounder Z (F26), or a contrast
# the precompute lacks.
#
# SCOPE (v1): GENE entity <-> a metadata_sum_norm phenotype (continuous or a named group
# contrast). Metabolite/contaminant/flux ENTITIES are deferred (RES-NO-RESOLVE).
#
# Returns "RES-OK;<n_rows>" (kg_association.csv) | "RES-NO" | "RES-NO-RESOLVE"
#   | "RES-NO-RESOLVER" | "RES-NO-DB" | "RES-NO-META" | "RES-NO-PHENO" | "RES-NO-CONTRAST".
# ==============================================================================

.kgDefaultCovs <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

kgAssociation <- function(var_id, phenotype, var_type = "auto",
                          covariates = "default", contrast = "", ref = "",
                          subset = "all", omics_filter = "all",
                          minN = "10", mode = "tool"){
  library(RSQLite); library(DBI)
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  .kgSetPaths()
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")
  if(!file.exists(db.path)) return("RES-NO-DB")

  r <- .kgResolveFeature(res, var_id, var_type)
  if(length(r$layers) == 0 || !identical(r$type, "gene")) return("RES-NO-RESOLVE")

  # phenotype: metadata_sum_norm column + inferred type
  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  if(!(phenotype %in% names(m))){
    j <- match(.kgNorm(phenotype), tolower(names(m))); if(is.na(j)) return("RES-NO-PHENO")
    phenotype <- names(m)[j]
  }
  ptype <- .kgInferType(m[[phenotype]]); if(is.na(ptype)) return("RES-NO-PHENO")

  # covariates (default 4, drop primary; or none; or user list -- kept only if present)
  if(covariates %in% c("none", "NA", "")){
    cov <- character(0)
  } else if(identical(covariates, "default")){
    cov <- .kgDefaultCovs[.kgDefaultCovs != phenotype]
  } else {
    cov <- trimws(strsplit(covariates, "[,;]")[[1]]); cov <- cov[cov != phenotype & nzchar(cov)]
  }
  cov <- cov[cov %in% names(m)]

  # contrast: continuous -> none; group -> "A-B" (A vs B) or "anova"
  if(ptype == "disc"){
    if(identical(contrast, "")) return("RES-NO-CONTRAST")
    if(identical(contrast, "anova")){
      contrast_use <- "anova"; ref_use <- ref
    } else {
      parts <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
      if(length(parts) == 2){ contrast_use <- parts[1]; ref_use <- parts[2] } else return("RES-NO-CONTRAST")
    }
  } else { contrast_use <- ""; ref_use <- "" }

  # design frame (record_id x [phenotype + covariates]), optionally donor-subset
  md <- .kgCoerce(m, unique(c(phenotype, cov))); rownames(md) <- as.character(m[["record_id"]])
  subset_note <- "all donors"
  # the donor row filter -- ONE implementation, in kg_common.R (.kgDonorSubset). NULL means the
  # filter could not be applied, and that REFUSES rather than returning the unrestricted analysis.
  if(!identical(subset, "all") && nzchar(subset)){
    msub <- .kgDonorSubset(m, subset)
    if(is.null(msub)) return("RES-NO-SUBSET")
    keep_ids <- as.character(msub[["record_id"]])
    md <- md[rownames(md) %in% keep_ids, , drop = FALSE]
    subset_note <- subset
    if(nrow(md) < minN) return("RES-NO")
  }
  # ⚠ narrowed on the SUBSETTED design frame: a filter can leave its own column constant ("within
  # males, adjusted for sex") and the fit drops it -- see `.kgVaryingCovs` in kg_common.R. `cov` is
  # also what the `Covariates` column reports, so without this the row claims an adjustment the
  # model never made. Inert without a subset: nothing is constant across the whole cohort.
  cov <- .kgVaryingCovs(md, cov)

  allow <- if(identical(omics_filter,"all")) NULL else trimws(strsplit(omics_filter, "\\|")[[1]])
  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)

  rows <- list()
  for(la in r$layers){
    if(!is.null(allow) && !(la$display %in% allow)) next
    de <- .kgOmicsDE(con, la$table, la$read_col, md, phenotype, ptype, ref_use, contrast_use)
    if(is.null(de) || nrow(de) == 0) next
    hit <- de[as.character(de$FeatureId) == la$read_val, , drop = FALSE]
    if(nrow(hit) == 0) hit <- de[toupper(as.character(de$Symbol)) == toupper(var_id), , drop = FALSE]
    if(nrow(hit) == 0) next
    hit <- hit[1, ]
    eff <- suppressWarnings(as.numeric(hit$Effect))
    rows[[length(rows)+1]] <- data.frame(
      Feature = var_id, Symbol = as.character(hit$Symbol), Layer = la$display, N = hit$N,
      Effect = signif(eff, 3), EffectType = hit$EffectType,
      P_value = signif(as.numeric(hit$P_value), 3), Adjusted_p_value = signif(as.numeric(hit$Adjusted_p_value), 3),
      Direction = ifelse(is.na(eff), "NA", ifelse(eff > 0, "up", "down")),
      Covariates = if(length(cov)) paste(cov, collapse = ",") else "none",
      Comparison = if(ptype == "disc") contrast else "continuous",
      Subset = subset_note, stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows); out <- out[order(out$Adjusted_p_value), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_association_", safe(var_id), "_", safe(phenotype), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_association.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}
