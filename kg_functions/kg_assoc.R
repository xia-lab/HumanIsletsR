# ==============================================================================
# kg_functions/kg_assoc.R  --  kgAssoc: the unified feature <-> phenotype ASSOCIATION
# entry point (op#1). ONE interface over the full input grid; TWO internal engines chosen
# by the FEATURE axis; orthogonal modifiers (covariates / subset / layers / type) applied
# the same way in both.
#
#   features    named set ("INS" | "INS,GCG,PDX1")  -> TARGETED engine (lm)
#               "all"                                -> SCREEN engine (limma, all features)
#               ANY entity: gene | metabolite | contaminant | reaction(flux) -- one engine,
#               same model; the resolver supplies each one's table + read key.
#   phenotypes  a leaf metadata column | a PARENT concept (e.g. "insulin secretion" /
#               "gsis") that expands to its child columns | a comma list | "all"
#   layers      "all" | comma list of layer keys. LAYERS ARE ENTITY-DRIVEN -- a gene covers
#               only the gene layers, a metabolite only the metabolite layers, etc.:
#                 gene        rnaseq | protein | nanostring | pbrna_alpha | pbrna_beta
#                 metabolite  metabolite_hg | metabolite_lg | metabolite_ratio
#                 contaminant contaminant  (= TWO layers, adipose + pancreas, LOD-corrected rows,
#                             reported as "Environmental contaminants (adipose|pancreas)")
#                 flux        flux_hg | flux_lg | flux_ratio   (supported, NOT default)
#               "all" on a TARGETED query = that feature's own layers (resolver-driven);
#               "all" on a SCREEN = the GENE layers -- name a layer to screen a metabolite /
#               contaminant / flux. OUT OF SCOPE: pbrna_delta (n=3), pbrna_pp (n=7).
#   covariates  "default" (age,sex,BMI,culturetime; primary dropped) | "none" | list
#               ⚠ purity removed 2026-09-04 to match the precompute default (dropped there
#               2026-08); the two must name the same set -- see kg_association.R's header.
#   subset      "all" | "col=value" (donor subset before fitting)
#   contrast    categorical only: "A-B" (A vs B) | "anova" | "" (auto: 2-level = A-B,
#               >=3 levels = omnibus F)
#   fdr_family  "query" (BH across exactly the tests asked -> P_family; NEVER genome-wide)
#               | "none" | "genomewide" (screen only: copies limma's genome-wide adj-p)
#   method      "auto" (DEFAULT) | "lm" | "limma" | "kendall"
#               auto = KENDALL when the layer is `contaminant`, otherwise the normal split
#               (targeted -> lm, screen -> limma). This mirrors the DEPLOYED trigger, which
#               is layer-based rather than a flag: DonorRegression switches to Kendall on
#               omicsType == "proc_contaminants" (humanislets_statistics.R:1878), because
#               small n and non-normal abundances make limma's moderated-t inappropriate
#               there. Naming a method overrides the choice.
#               ⚠ Kendall's effect is TAU, so EffectType is "tau" (or "kruskalH" for a >=3
#               level omnibus, which has no per-feature effect size) -- never "log2FC".
#
# TARGETED = lm(value ~ phenotype + covariates): cont -> slope; disc 2-level/named -> mean
# difference; disc >=3 -> covariate-adjusted omnibus ANOVA F. Reads only the queried
# feature rows (one IN-query per layer). This is the fast path for one/few features.
# SCREEN  = limma (trend+robust eBayes) over ALL features, one fit per phenotype, reading
# each layer once. This is the genome-wide differential (feeds enrichment etc.).
#
# Retrieval note: kg_functions are on-demand COMPUTE. The "is this in the precompute ->
# retrieve" decision lives in the caller (it reads a different store, HI_precomputed);
# kgAssoc always COMPUTES the requested spec. all-features x all-phenotypes recompute is
# refused (RES-USE-PRECOMPUTE) unless narrowed, because that IS the precompute.
#
# Writes kg_assoc.csv (one row per feature x phenotype x layer) and kg_assoc_layers.csv (one row per layer x
# phenotype the engines reached: donors, usable donors, tested | too_few_donors). Returns "RES-OK;<n_rows>"
#   | "RES-USE-PRECOMPUTE" | "RES-NO" | "RES-NO-DB" | "RES-NO-META" | "RES-NO-PHENO"
#   | "RES-NO-FEATURE" | "RES-NO-RESOLVE" | "RES-NO-RESOLVER".
# ==============================================================================

.kgAssocDefaultCovs <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# ENTITY SCOPE: gene, metabolite, contaminant and reaction(flux) all run through the SAME
# engine -- the model is identical (lm on the feature's donor vector); only the source table
# differs, and the resolver already carries it. LAYERS ARE ENTITY-DRIVEN: a gene only ever
# covers the gene layers, a metabolite only the metabolite layers, a contaminant/reaction its
# own -- enforced because the resolver returns a feature's OWN layers, filtered to the
# supported keys below. FLUX is SUPPORTED but NOT a default layer (user, 2026-07-16): query a
# reaction token or name layers="flux_hg" to reach it; it never enters a gene/screen default.
# OUT OF SCOPE (user, 2026-07-16): pbrna_delta (n=3) and pbrna_pp (n=7) -- below minN, so
# always empty.
.kgScreenLayer <- function(name){
  switch(tolower(name),
    "rnaseq"           = c("proc_rnaseq",                  "accession"),
    "protein"          = c("proc_prot_combat",             "id"),
    "nanostring"       = c("proc_nanostring_merge",        "gene_id"),
    "pbrna_alpha"      = c("proc_pbrna_Alpha",             "symbol"),
    "pbrna_beta"       = c("proc_pbrna_Beta",              "symbol"),
    "metabolite_hg"    = c("proc_metabolite_combat_HG",    "compound"),
    "metabolite_lg"    = c("proc_metabolite_combat_LG",    "compound"),
    "metabolite_ratio" = c("proc_metabolite_combat_ratio", "compound"),
    "contaminant"      = c("proc_contaminants",            "Compound"),
    "flux_hg"          = c("proc_flux_HG",                 "rxn"),
    "flux_lg"          = c("proc_flux_LG",                 "rxn"),
    "flux_ratio"       = c("proc_flux_ratio",              "rxn"),
    NULL)
}

# SUPPORTED layer keys, grouped by entity. Note the split between SUPPORTED and DEFAULT:
# a SCREEN (features="all") defaults to the GENE layers only (what "run the differential
# across the omics" means, and what feeds enrichment); metabolite / contaminant / flux are
# SUPPORTED but never default -- reached by querying that entity's token (layers are
# entity-driven) or by naming the layer (e.g. layers="flux_hg").
.kgGeneLayerKeys        <- c("rnaseq","protein","nanostring","pbrna_alpha","pbrna_beta")
.kgMetaboliteLayerKeys  <- c("metabolite_hg","metabolite_lg","metabolite_ratio")
.kgContaminantLayerKeys <- c("contaminant")
.kgFluxLayerKeys        <- c("flux_hg","flux_lg","flux_ratio")

# ⚠ THE CONTAMINANT KEY IS TWO LAYERS -- ONE PER TISSUE (user 2026-10-01, "yes, go ahead"). MEASURED on F11 "What about
# in female donors?": 59 rows labelled `contaminant` = 57 chemicals' ADIPOSE values + the 2 measured only in pancreas,
# pancreas never tested for the rest. The key the callers send stays `contaminant`; it fans out below into both
# tissues (`.kgContaminantTissues` / `.kgTissueLayers`, kg_common.R -- shared with every other kg function).
.kgAllLayerKeys <- c(.kgGeneLayerKeys, .kgMetaboliteLayerKeys,
                     .kgContaminantLayerKeys, .kgFluxLayerKeys)

# ⚠ STAGED LAYERS -- A CONCEPT SCREEN RUNS ITS MAIN LAYERS AND OFFERS THE REST (user 2026-10-04: "in default we do only at
# most three layer Bulk RNA-seq, Protein, Metabolomics HG / LG / ratio, and tell user what is provide but not calculate
# yet ... only if both Bulk RNA-seq, Protein, has no data ... try other gene layers, if all Metabolomics layer no data, go
# Contaminants"). MEASURED on F4 "top predictors of insulin release, adjusted for islet purity" (20 measures x 9 layers):
# 1,386,571 tests, 132.8 s, of which ~108 s fitting -- pseudobulk alpha / beta alone were 792k of the tests.
# The spec rides in the existing `layers` argument (no new parameter): groups split by "|", each "<primary>/<fallback>",
# e.g. "rnaseq,protein/nanostring,pbrna_alpha,pbrna_beta|metabolite_hg,metabolite_lg,metabolite_ratio/contaminant".
# A group runs its primary layers; ONLY when no primary layer has enough donors for this call does it run the fallback
# instead. A fallback layer that is not run but has donors is OFFERED: written to the status file as `not_calculated`,
# so the caller can say it can still be calculated. A `layers` value without "/" is read exactly as before.
# Donors are counted from the layer table's donor COLUMNS (`dbListFields`, no data read) -- the same donors
# `.kgLoadLayerMatrix` keeps -- intersected with the subset, and per phenotype those with a value (= `.note_layer`).
.kgLayerDonors <- function(con, key){
  lp <- .kgScreenLayer(key); if(is.null(lp)) return(character(0))
  f <- try(DBI::dbListFields(con, lp[1]), silent = TRUE)
  if(inherits(f, "try-error")) return(character(0))
  grep("^R[0-9]+$", f, value = TRUE)
}
.kgLayerUsable <- function(con, key, keep_ids, phframe, phenos, minN){
  dc <- intersect(.kgLayerDonors(con, key), keep_ids)
  if(length(dc) < minN) return(FALSE)
  rows <- intersect(dc, rownames(phframe))
  any(vapply(phenos, function(ph){
    v <- phframe[rows, ph]; sum(!is.na(v) & !(as.character(v) %in% c("", "NA"))) >= minN }, logical(1)))
}
.kgStagedLayers <- function(con, spec, keep_ids, phframe, phenos, minN){
  run <- character(0); offered <- character(0)
  for(g in strsplit(spec, "|", fixed = TRUE)[[1]]){
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    prim <- tolower(trimws(strsplit(parts[1], "[,;]")[[1]])); prim <- prim[nzchar(prim)]
    fb <- if(length(parts) > 1) tolower(trimws(strsplit(parts[2], "[,;]")[[1]])) else character(0)
    fb <- fb[nzchar(fb)]
    ok <- vapply(prim, function(k) .kgLayerUsable(con, k, keep_ids, phframe, phenos, minN), logical(1))
    if(any(ok)) { run <- c(run, prim); offered <- c(offered, fb) } else run <- c(run, prim, fb)
  }
  run <- unique(run)
  offered <- setdiff(unique(offered), run)
  offered <- offered[vapply(offered, function(k) .kgLayerUsable(con, k, keep_ids, phframe, phenos, minN), logical(1))]
  list(run = paste(run, collapse = ","), offered = offered)
}

# short layer key from a resolver display name, so `layers` takes short codes in BOTH engines
# (e.g. "Bulk gene expression (RNA-seq)" -> "rnaseq"; "Metabolite (LG)" -> "metabolite_lg").
.kgLayerKey <- function(display){
  d <- tolower(display)
  if(grepl("methylation", d)) return(if(grepl("promoter", d)) "methylation_promoter" else "methylation")
  if(grepl("contaminant", d)) return("contaminant")
  if(grepl("metabolic flux", d))
    return(if(grepl("ratio", d)) "flux_ratio" else if(grepl("\\(lg\\)", d)) "flux_lg" else "flux_hg")
  if(grepl("metabolite", d))
    return(if(grepl("ratio", d)) "metabolite_ratio" else if(grepl("\\(lg\\)", d)) "metabolite_lg" else "metabolite_hg")
  if(grepl("pseudobulk", d)){
    if(grepl("alpha", d)) return("pbrna_alpha")
    if(grepl("beta",  d)) return("pbrna_beta")
    if(grepl("delta", d)) return("pbrna_delta")
    if(grepl("pp",    d)) return("pbrna_pp")
    return("pbrna")
  }
  if(grepl("protein", d))        return("protein")
  if(grepl("nanostring", d))     return("nanostring")
  if(grepl("rna-seq|rnaseq", d)) return("rnaseq")
  gsub("[^a-z0-9]+", "_", d)
}

# Split a feature list. Metabolite/contaminant NAMES legitimately contain commas
# (58 tokens, e.g. "1,2,3,4,6,7,8-HPCDD", "2,4,7,9-tetramethyl-5-decyne-4,7-diol"), while NO
# token contains ";" or "|". So: ";"/"|" are ALWAYS separators; "," is a separator only when
# the whole string is not itself one resolvable token. Keeps "INS,GCG,PDX1" working while
# never shredding a comma-bearing compound name.
.kgSplitFeatures <- function(s, res){
  if(grepl("[;|]", s)) return(trimws(strsplit(s, "[;|]")[[1]]))
  if(nrow(res[res$token == .kgNorm(s), , drop = FALSE]) > 0) return(s)
  trimws(strsplit(s, ",", fixed = TRUE)[[1]])
}

.kgLayerEntity <- function(key){
  if(key %in% c("metabolite_hg","metabolite_lg","metabolite_ratio")) "metabolite"
  else if(identical(key, "contaminant")) "contaminant"
  else if(key %in% c("flux_hg","flux_lg","flux_ratio")) "reaction"
  else "gene"
}

# Entity types whose donor n is too small to spend the default covariate DoF (metabolite ~37,
# contaminants ~31) -> the DEFAULT covariate set is dropped for them (an explicit user
# covariate list is still honored). Flux is NOT here: it has 371 donors, a full-size layer.
.kgSmallNTypes <- c("metabolite", "contaminant")

# lm -> limma switch for a NAMED feature set (user-set 2026-07-16). Past this many named
# features the query is a screen in all but name, so we fit the FULL layer and extract the
# named rows: moderated stats (prior from ALL features, not the named subset) instead of k
# ordinary-t lm's. NOT a statistical cliff -- nothing changes at exactly 500. Measured on
# rnaseq (1 thread): a full-layer limma is a FLAT ~3.3s at any k, while k lm's cost ~0.0036s
# each, so the two cross at k ~= 870. 500 sits below that crossover on purpose: the switch
# costs ~1.4s there but buys the better statistic, and it stays far above any real named set
# (1-50), so interactive queries keep the ~16x faster lm path.
.kgLmMaxFeatures <- 500L

# resolve the fitting method for ONE layer.
# ⚠ `auto` reproduces the DEPLOYED trigger, which is LAYER-based, not a flag: DonorRegression
# switches to Kendall on `omicsType == "proc_contaminants"`
# (humanislets_statistics.R:1878). So contaminant x phenotype gets Kendall and every other
# layer keeps the existing targeted-lm / screen-limma split -- no other frame changes.
# Scope is exactly that pairing: a collaborator decision that is owned, never derived from
# the data (deriving it as "small-n or tie-heavy" would silently capture metabolomics and
# methylation, which the decision never covered).
.kgAssocMethod <- function(method, layer_key, is_screen){
  mm <- tolower(trimws(if(is.null(method) || !nzchar(method)) "auto" else method))
  if(mm %in% c("lm", "limma", "kendall")) return(mm)
  if(identical(tolower(layer_key), "contaminant")) return("kendall")
  if(is_screen) "limma" else "lm"
}

# ---- named-set test: ONE self-contained P for the requested features as a group (set_test = "fry") -----
# Same data, donors, covariates and model as the per-feature rows (feature ~ phenotype + covariates); the
# members are then tested TOGETHER with limma::fry. Member handling is the case-study scripts' own `prep`
# (p10/p11, donor_clustering_pipeline_v2/outputs_basal_floor_case_study_20261011): a member enters when it
# is measured in >= max(9, 90%) of the donors and varies; its remaining gaps take its own mean; >= 3 members.
# A categorical phenotype is tested only as a two-level comparison (`lvl` vs `ref`, "rest" allowed).
.kgFrySet <- function(X, phv, cvd, ptype, lvl = "", ref = "", minN = 10){
  don <- intersect(colnames(X), names(phv))
  if(!is.null(cvd) && ncol(cvd)) don <- intersect(don, rownames(cvd))
  if(length(don) < minN) return(NULL)
  d <- data.frame(.ph = phv[don], row.names = don, stringsAsFactors = FALSE)
  if(!is.null(cvd) && ncol(cvd)) d <- cbind(d, cvd[don, , drop = FALSE])
  if(ptype == "cont"){
    d$.ph <- suppressWarnings(as.numeric(as.character(d$.ph)))
  } else {
    if(!nzchar(lvl) || !nzchar(ref) || identical(lvl, "anova")) return(NULL)
    v <- as.character(d$.ph); v[v %in% c("", "NA")] <- NA
    v <- if(identical(ref, "rest")) ifelse(is.na(v), NA, ifelse(v == lvl, lvl, "rest")) else ifelse(v %in% c(lvl, ref), v, NA)
    d$.ph <- factor(v, levels = c(ref, lvl))
  }
  d <- droplevels(d[stats::complete.cases(d), , drop = FALSE])     # only the levels these donors carry (full-rank design)
  if(nrow(d) < minN || (ptype != "cont" && min(table(d$.ph)) < 2)) return(NULL)
  Xs <- X[, rownames(d), drop = FALSE]; storage.mode(Xs) <- "double"
  k <- rowSums(!is.na(Xs)) >= max(9, 0.9 * ncol(Xs)) &
       apply(Xs, 1, function(v) isTRUE(stats::sd(v, na.rm = TRUE) > 0))
  Xs <- Xs[k, , drop = FALSE]; if(nrow(Xs) < 3) return(NULL)
  if(anyNA(Xs)) Xs[is.na(Xs)] <- rowMeans(Xs, na.rm = TRUE)[row(Xs)[is.na(Xs)]]
  cn <- .kgVaryingCovs(d, setdiff(names(d), ".ph"))
  D  <- stats::model.matrix(stats::as.formula(paste0("~ .ph", if(length(cn))
          paste0(" + ", paste(sprintf("`%s`", cn), collapse = " + ")) else "")), data = d)
  j  <- which(colnames(D) == if(ptype == "cont") ".ph" else paste0(".ph", lvl))
  if(length(j) != 1) return(NULL)
  fr <- try(limma::fry(Xs, index = list(set = seq_len(nrow(Xs))), design = D, contrast = j, sort = "none"),
            silent = TRUE)
  if(inherits(fr, "try-error") || !NROW(fr)) return(NULL)
  list(N = nrow(d), NGenes = nrow(Xs), Direction = as.character(fr$Direction[1]), P = as.numeric(fr$PValue[1]),
       NPred = sum(rownames(d) %in% .GlobalEnv$.kg_pred_ids))
}

.kgSetRow <- function(s, layer, ph, cvn, pt, contrast, subset_note, method = "fry")
  data.frame(Layer = layer, Phenotype = ph, N = s$N, N_members = s$NGenes, Direction = s$Direction,
             P_value = signif(s$P, 4), Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
             Comparison = if(pt == "disc") contrast else "continuous", Subset = subset_note, Method = method,
             Effect = if(is.null(s$Effect)) NA_real_ else signif(s$Effect, 4),
             N_layer = if(is.null(s$NAll)) NA_integer_ else as.integer(s$NAll),
             N_predicted = if(is.null(s$NPred)) 0L else as.integer(s$NPred),
             stringsAsFactors = FALSE)

# ---- named-set test "compared with all other features" (set_test = "relative", 2026-10-19) -----------------------------
# For a question that asks whether a set is lower / higher THAN THE OTHER genes / promoters / proteins of a layer. One score
# per donor: every feature of the layer standardised across the tested donors (z); the set's mean z minus the mean z of all
# the layer's features; scaled; then score ~ phenotype + covariates (lm), on the same donors, covariates and contrast as fry.
# The set score is the published z-score score (Lee et al. 2008, PLoS Comput Biol 4:e1000217; GSVA `zscore`) — verified
# identical to GSVA 1.52.3 on KEGG Hippo promoter methylation (r = 1.000). Subtracting the layer-wide mean is this project's
# step: it is what makes the test "than the rest", and it removes a donor-wide shift (MEASURED 2026-10-19: the donors'
# all-feature mean varies SD 0.275 z for promoter methylation, 0.059 for RNA). Supplementary Figure S5A's methylation score
# is this score without the per-feature standardisation (Hippo, 24 omics-clustered donors: -1.30 SD P 0.012 raw,
# -1.24 SD P 0.031 standardised; leave-one-gene-out all P < 0.05). Features enter as for fry (>= max(9, 90%) of the donors,
# varying; remaining gaps take the feature's own mean); >= 3 members. A CpG layer is first summarised per gene (mean of
# the gene's CpGs, S5A's step), so the set and the rest are genes: `genes` = each row's gene list, `members` = symbols.
.kgRelMethod <- "z-score set score vs all features"
.kgRelScoreSet <- function(Xall, members, phv, cvd, ptype, lvl = "", ref = "", minN = 10, genes = NULL){
  don <- intersect(colnames(Xall), names(phv))
  if(!is.null(cvd) && ncol(cvd)) don <- intersect(don, rownames(cvd))
  if(length(don) < minN) return(NULL)
  d <- data.frame(.ph = phv[don], row.names = don, stringsAsFactors = FALSE)
  if(!is.null(cvd) && ncol(cvd)) d <- cbind(d, cvd[don, , drop = FALSE])
  if(ptype == "cont"){
    d$.ph <- suppressWarnings(as.numeric(as.character(d$.ph)))
  } else {
    if(!nzchar(lvl) || !nzchar(ref) || identical(lvl, "anova")) return(NULL)
    v <- as.character(d$.ph); v[v %in% c("", "NA")] <- NA
    v <- if(identical(ref, "rest")) ifelse(is.na(v), NA, ifelse(v == lvl, lvl, "rest")) else ifelse(v %in% c(lvl, ref), v, NA)
    d$.ph <- factor(v, levels = c(ref, lvl))
  }
  d <- droplevels(d[stats::complete.cases(d), , drop = FALSE])
  if(nrow(d) < minN || (ptype != "cont" && min(table(d$.ph)) < 2)) return(NULL)
  X <- Xall[, rownames(d), drop = FALSE]; storage.mode(X) <- "double"
  if(!is.null(genes)){
    gi <- rep(seq_along(genes), lengths(genes)); gg <- toupper(trimws(unlist(genes)))
    ok <- !is.na(gg) & nzchar(gg); gi <- gi[ok]; gg <- gg[ok]
    if(!length(gg)) return(NULL)
    X <- rowsum(X[gi, , drop = FALSE], gg, na.rm = TRUE) / as.vector(rowsum(rep(1, length(gg)), gg))
    members <- toupper(trimws(as.character(members)))
  }
  k <- rowSums(!is.na(X)) >= max(9, 0.9 * ncol(X)) & apply(X, 1, function(v) isTRUE(stats::sd(v, na.rm = TRUE) > 0))
  X <- X[k, , drop = FALSE]
  s_in <- rownames(X) %in% members
  if(sum(s_in) < 3 || nrow(X) <= sum(s_in)) return(NULL)
  if(anyNA(X)) X[is.na(X)] <- rowMeans(X, na.rm = TRUE)[row(X)[is.na(X)]]
  Z <- t(scale(t(X)))
  d$.score <- as.numeric(scale(colMeans(Z[s_in, , drop = FALSE]) - colMeans(Z)))
  cn <- .kgVaryingCovs(d, setdiff(names(d), c(".ph", ".score")))
  f <- try(stats::lm(stats::as.formula(paste0(".score ~ .ph", if(length(cn))
         paste0(" + ", paste(sprintf("`%s`", cn), collapse = " + ")) else "")), data = d), silent = TRUE)
  if(inherits(f, "try-error")) return(NULL)
  cf <- summary(f)$coefficients
  j <- if(ptype == "cont") ".ph" else paste0(".ph", lvl)
  if(!(j %in% rownames(cf))) return(NULL)
  list(N = nrow(d), NGenes = sum(s_in), NAll = nrow(X), Direction = if(cf[j, 1] < 0) "Down" else "Up",
       P = as.numeric(cf[j, 4]), Effect = as.numeric(cf[j, 1]),
       NPred = sum(rownames(d) %in% .GlobalEnv$.kg_pred_ids),
       # the score the model fitted, one row per donor (2026-10-19: the IsletAgent donor plot draws them, as S5A does)
       Scores = data.frame(record_id = rownames(d), group = as.character(d$.ph), score = signif(d$.score, 5),
                           predicted = rownames(d) %in% .GlobalEnv$.kg_pred_ids, stringsAsFactors = FALSE))
}

# one row per donor of a relative set test, keyed like its set row (layer, phenotype); NULL for fry (no per-donor score)
.kgSetScores <- function(s, layer, ph){
  if(is.null(s) || is.null(s$Scores) || !nrow(s$Scores)) return(NULL)
  data.frame(Layer = layer, Phenotype = ph, s$Scores, stringsAsFactors = FALSE)
}

.kgAssocImpl <- function(features = "", phenotypes = "",
                    layers = "all", covariates = "default",
                    subset = "all", contrast = "", ref = "",
                    fdr_family = "query", minN = "10", method = "auto",
                    mode = "tool", set_test = ""){
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite"); if(!file.exists(db.path)) return("RES-NO-DB")

  # parse a group contrast "A-B" once (targeted path passes level A + ref B to .kgLmCell;
  # "" = auto, "anova" = omnibus). "A-B" means A vs B.
  ct_lvl <- contrast; ct_ref <- ref
  if(nzchar(contrast) && !identical(contrast, "anova")){
    pp <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if(length(pp) == 2){ ct_lvl <- pp[1]; ct_ref <- pp[2] }
  }

  # ⚠ "cluster_v2.<level>" = the multi-omics cluster WITH predicted donors (2026-10-19, user: "use all available donors
  # here"; the case study's donor set). The level is read off the name, the metadata is loaded at that level, and the
  # analysis runs on the plain column — rows report "cluster_v2" as before. Any other phenotype is untouched.
  cv2lv <- NULL
  if(exists(".cluster_v2_level", mode = "function")){
    .pt0 <- trimws(strsplit(phenotypes, "[,;]")[[1]])
    .lv0 <- unlist(lapply(.pt0, .cluster_v2_level))
    if(length(.lv0)){
      cv2lv <- .lv0[1]
      phenotypes <- paste(vapply(.pt0, .cluster_v2_base, character(1)), collapse = ",")
    }
  }
  m <- .kgLoadMetaFrame(level = cv2lv); if(is.null(m)) return("RES-NO-META")
  meta_names <- names(m); rec <- as.character(m[["record_id"]])

  # ---- phenotypes: expand parents -> leaves; "all" -> all catalogue leaves ----
  ptoks   <- trimws(strsplit(phenotypes, "[,;]")[[1]]); ptoks <- ptoks[nzchar(ptoks)]
  if(length(ptoks) == 0) return("RES-NO-PHENO")
  all_pheno <- (length(ptoks) == 1 && tolower(ptoks) == "all")
  phenos <- if(all_pheno) .kgAllPhenotypes(meta_names) else .kgExpandPhenotypes(ptoks, meta_names)
  phenos <- phenos[phenos %in% meta_names]; if(length(phenos) == 0) return("RES-NO-PHENO")

  # ---- feature axis -> engine (the list is split later, once the resolver is loaded) ----
  feat_raw  <- trimws(features); if(!nzchar(feat_raw)) return("RES-NO-FEATURE")
  is_screen <- identical(tolower(feat_raw), "all")

  # ---- guard (default #3): all features x all phenotypes recompute -> use precompute ----
  if(is_screen && all_pheno && identical(covariates, "default") && identical(subset, "all"))
    return("RES-USE-PRECOMPUTE")

  # ---- donor subset (applied in both engines) ----
  keep_ids <- rec; subset_note <- "all donors"
  # -- COMPARISONS, NOT JUST EQUALITY (2026-09-04), via the ONE parser in kg_common.R. This block
  #    tested grepl("=", subset) and split on "=", so a filter like "bodymassindex>30" fell straight
  #    through: `subset_note` stayed "all donors" and the analysis ran on EVERY donor while the
  #    caller believed the restriction had been applied. The thresholds come from the question
  #    ("donors over 60", "a BMI over 30"); none is invented here.
  # -- AN UNAPPLIABLE FILTER NOW REFUSES. A column absent from the frame used to skip the block and
  #    return the unrestricted analysis as the answer to a restricted question.
  if(!identical(subset, "all") && nzchar(subset)){
    msub <- .kgDonorSubset(m, subset)
    if(is.null(msub)) return("RES-NO-SUBSET")
    keep_ids <- as.character(msub[["record_id"]]); subset_note <- subset
  }
  keep_ids <- keep_ids[!is.na(keep_ids)]

  # ---- covariate universe (coerced once), and a per-phenotype selector (drop primary) ----
  covuniv <- if(covariates %in% c("none","NA","")) character(0)
             else if(identical(covariates,"default")) .kgAssocDefaultCovs
             else trimws(strsplit(covariates, "[,;]")[[1]])
  covuniv <- unique(covuniv[nzchar(covuniv) & covuniv %in% meta_names])
  covbase <- if(length(covuniv)){ cb <- .kgCoerce(m, covuniv); rownames(cb) <- rec; cb } else NULL
  covcols_for <- function(ph) setdiff(covuniv, ph)

  # ---- phenotype values (coerced once) ----
  phframe <- .kgCoerce(m, phenos); rownames(phframe) <- rec

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  rows <- list()
  set_rows <- list()
  set_scores <- list()                                              # per-donor set scores (relative test only)
  set_kind <- tolower(trimws(set_test))
  do_set <- set_kind %in% c("fry", "relative") && !is_screen       # "relative": `.kgRelScoreSet`, the set vs all features
  # ⚠ A LAYER SKIPPED FOR TOO FEW DONORS WAS SILENT (user 2026-10-01, "yes, go ahead"). Both engines `next` past a
  # layer whose donors (layer x subset) fall below `minN`, and nothing said so -- MEASURED on F4 "Which beta-cell genes
  # … Tamil Sri Lankan, in female donors only?": 6 female donors carry beta-cell data, no beta row came back, and the
  # answer called it "a MEASURED null, not an absence of data". Every layer an engine reaches is recorded here with its
  # donor count and, per phenotype, the donors that also carry a value; written to `kg_assoc_layers.csv` beside
  # `kg_assoc.csv`. Status `too_few_donors` = skipped at the `minN` gate, `tested` = went to the fit.
  lstat <- list()
  offered <- character(0)        # layers a staged spec offers but does not run (`.kgStagedLayers`)
  # ⚠ A NAMED COMPARISON COUNTS ITS TWO GROUPS (user 2026-10-07, "1b"). The targeted fit (`.kgLmCell`) keeps only the
  # donors at the two compared levels and returns nothing unless both are present and they number >= minN together
  # (one-vs-rest: >= 3 in the named group first) -- so a layer could be "tested" on 218 donors and give no row.
  # MEASURED on F2 "How does CPA1 change in type 1 diabetes, in high-purity preparations only?": RNA-seq 218 donors, 0
  # with type 1 diabetes; the answer said "not measured in any layer". With `groups = TRUE` (the targeted engine only;
  # limma has its own one-class guard) a discrete phenotype asked as a named contrast writes the donors at each level
  # (`N_level`, `N_ref`), `N_usable` = the two together, and `too_few_donors` under exactly `.kgLmCell`'s conditions.
  # Every other row is written as before, with `N_level` / `N_ref` NA.
  .note_layer <- function(disp, dc, status = NULL, groups = FALSE){
    for(ph in phenos){
      v <- phframe[intersect(dc, rownames(phframe)), ph]
      ok <- !is.na(v) & !(as.character(v) %in% c("", "NA"))
      nu <- sum(ok); nl <- NA_integer_; nr <- NA_integer_
      st <- if(!is.null(status)) status else if(length(dc) < minN) "too_few_donors" else "tested"
      if(groups && nzchar(ct_lvl) && !identical(ct_lvl, "anova") && nzchar(ct_ref) && .kgPhenoType(m, ph) == "disc"){
        vv <- as.character(v[ok])
        nl <- sum(vv == ct_lvl)
        nr <- if(identical(ct_ref, "rest")) sum(vv != ct_lvl) else sum(vv == ct_ref)
        nu <- nl + nr
        if(identical(st, "tested") && ((identical(ct_ref, "rest") && nl < 3) || nl == 0 || nr == 0 || nu < minN))
          st <- "too_few_donors"
      }
      lstat[[length(lstat)+1]] <<- data.frame(
        Layer = disp, Phenotype = ph, N_layer = length(dc),
        N_usable = nu, MinN = minN, Status = st,
        Subset = subset_note, N_level = nl, N_ref = nr, stringsAsFactors = FALSE)
    }
  }
  # ⚠ ONE CLASS IS A DATA CONDITION, NOT A CRASH. A donor subset can leave a discrete phenotype with
  # a single level ("type 1 diabetes vs pancreatitis, in lean donors only"), and limma then dies
  # inside model.matrix with "contrasts can be applied only to factors with 2 or more levels".
  # MEASURED 2026-09-22 on F10: 4 records came back RES-ERR — the guard caught the crash, which is
  # right, but a crash and "the subset has only one class" are different answers and a reader cannot
  # act on the first. Tracked here so the return can say WHICH it was.
  one.class <- FALSE
  # the fit sets this when ITS donor intersection leaves one level — see `.kgLimmaFit`. Reset per
  # call, because an Rserve session is reused across requests.
  .GlobalEnv$.kg_one_class <- FALSE

  parse_contrast <- function(){
    if(!nzchar(contrast) || identical(contrast, "anova")) return(list(c = "anova", r = ref))
    p <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if(length(p) == 2) list(c = p[1], r = p[2]) else list(c = "anova", r = ref)
  }

  # ---- resolve the named features (ANY entity) BEFORE choosing the engine ----
  feats <- list(); ftoks <- character(0)
  # DNA methylation keys this call NAMES (never a default layer; see `.kgLoadMethylation`, kg_common.R). They are fitted in
  # their own block below, after the two engines, so every other layer keeps exactly the path it had.
  # read through `.kgLayerKey`, as every display name is: the agent's trained call carries the question's own word
  # ("promoter methylation"), the fired call the key ("methylation_promoter") -- both name the same layer
  meth_keys <- if(identical(layers, "all") || grepl("/", layers, fixed = TRUE)) character(0) else {
    .lt <- trimws(strsplit(layers, "[,;]")[[1]]); .lt <- .lt[nzchar(.lt)]
    intersect(unique(vapply(.lt, .kgLayerKey, character(1), USE.NAMES = FALSE)), .kgMethLayerKeys) }
  if(!is_screen){
    res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
    ftoks <- .kgSplitFeatures(feat_raw, res); ftoks <- ftoks[nzchar(ftoks)]
    if(length(ftoks) == 0) return("RES-NO-FEATURE")
    lay_filter <- if(identical(layers, "all")) NULL else tolower(trimws(strsplit(layers, "[,;]")[[1]]))

    # resolve each feature -> its layers. ANY entity type (gene | metabolite | contaminant |
    # reaction) -- the resolver carries the source table, and the model is the same.
    # ONE vectorised pass for the whole set (a per-token scan costs ~9ms each).
    rmap <- .kgResolveFeatures(res, ftoks, "auto")
    for(tk in ftoks){
      r <- rmap[[.kgNorm(tk)]]
      if(is.null(r) || length(r$layers) == 0 || is.na(r$type)) next
      for(la in r$layers){
        lk <- .kgLayerKey(la$display)
        # entity-driven + supported-only: the resolver hands back this feature's OWN layers,
        # and anything out of scope (pbrna_delta/pbrna_pp/flux) is dropped here.
        if(!(lk %in% .kgAllLayerKeys)) next
        # the resolver's one contaminant layer -> one layer per tissue (`.kgTissueLayers`, kg_common.R)
        for(la2 in .kgTissueLayers(la)){
          if(!is.null(lay_filter) && !(lk %in% lay_filter || .kgDisplayIn(la2, lay_filter))) next
          feats[[length(feats)+1]] <- list(token = tk, etype = r$type, layer = la2)
        }
      }
    }
    # methylation features are not in the resolver (CpGs are read off the matrix's own annotation), so a call that names
    # only methylation layers goes on to the methylation block
    if(length(feats) == 0 && !length(meth_keys)) return("RES-NO-RESOLVE")
  }

  # ---- ENGINE CHOICE: named set -> lm, until it is big enough to be a screen in all but
  # name (.kgLmMaxFeatures), past which the full-layer limma fit is both better and faster.
  use_limma <- is_screen || (length(ftoks) >= .kgLmMaxFeatures)

  if(!use_limma){
    # =========================== TARGETED engine (lm) ===========================
    # group by layer table so each layer is read ONCE (IN-query over its features)
    key <- vapply(feats, function(f) paste(f$layer$display, f$layer$table, f$layer$read_col, sep = "@@"), character(1))
    for(k in unique(key)){
      idx <- which(key == k); la <- feats[[idx[1]]]$layer
      read_vals <- unique(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1)))
      FR <- .kgFeatureRows(con, la$table, la$read_col, read_vals, la$rowfilter); if(is.null(FR)) next
      dcommon <- intersect(colnames(FR$mat), keep_ids)
      .note_layer(la$display, dcommon, groups = TRUE)
      if(length(dcommon) < minN) next

      # small-n entities (metabolite/contaminant) drop the DEFAULT covariates; an explicit
      # user list is always honored. Entity is a property of the layer, so this is per-layer.
      small <- identical(covariates, "default") && (feats[[idx[1]]]$etype %in% .kgSmallNTypes)
      # The design pieces depend only on (layer, phenotype) -- NOT on the feature -- so build
      # them ONCE per layer instead of re-subsetting per feature (that churn dominated a
      # many-feature query far more than the lm's themselves).
      P <- list()
      for(ph in phenos){
        cvn <- if(small) character(0) else covcols_for(ph)
        # ⚠ A DONOR SUBSET CAN LEAVE A COVARIATE CONSTANT ON THESE DONORS ("within males, adjusted
        # for sex"), and the fit drops it -- see `.kgVaryingCovs` in kg_common.R. Narrow it HERE
        # too, because `cvn` is also what the `Covariates` column reports: otherwise the row claims
        # an adjustment the model never made.
        if(length(cvn) && !is.null(covbase)) cvn <- .kgVaryingCovs(covbase[dcommon, , drop = FALSE], cvn)
        phv <- phframe[dcommon, ph]; names(phv) <- dcommon
        P[[ph]] <- list(pt = .kgPhenoType(m, ph), cvn = cvn, phv = phv,
                        cvd = if(length(cvn) && !is.null(covbase)) covbase[dcommon, cvn, drop = FALSE] else NULL)
      }
      rpos <- match(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1)), FR$id)
      # method depends only on the LAYER, so resolve it once here. `la$display` is the display
      # name, so map it back to the short key the resolver understands.
      meth    <- .kgAssocMethod(method, .kgLayerKey(la$display), is_screen)
      cellfun <- if(identical(meth, "kendall")) .kgKendallCell else .kgLmCell

      for(j in seq_along(idx)){
        f  <- feats[[idx[j]]]; ri <- rpos[j]; if(is.na(ri)) next
        y  <- FR$mat[ri, dcommon]; names(y) <- dcommon; sym <- FR$symbol[ri]
        for(ph in phenos){
          pp <- P[[ph]]; pt <- pp$pt; cvn <- pp$cvn
          cell <- cellfun(y, pp$phv, pp$cvd, pt, ct_lvl, ct_ref, minN); if(is.null(cell)) next
          rows[[length(rows)+1]] <- data.frame(
            Feature = f$token, Symbol = sym, Entity = f$etype, Layer = la$display, Phenotype = ph,
            N = cell$N, Effect = signif(cell$Effect, 3), EffectType = cell$EffectType,
            Statistic = signif(cell$t, 3), P_value = signif(cell$P, 4),
            Direction = cell$Direction,
            Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
            Comparison = if(pt == "disc") (if(nzchar(contrast)) contrast else paste0("auto(", cell$EffectType, ")")) else "continuous",
            Subset = subset_note, Method = meth, stringsAsFactors = FALSE)
        }
      }
      if(do_set && identical(set_kind, "relative")){
        LA <- .kgLoadLayerMatrix(con, la$table, la$read_col, la$rowfilter)
        if(!is.null(LA)) for(ph in phenos){
          s <- .kgRelScoreSet(LA$mat[, intersect(colnames(LA$mat), dcommon), drop = FALSE], read_vals,
                              P[[ph]]$phv, P[[ph]]$cvd, P[[ph]]$pt, ct_lvl, ct_ref, minN)
          if(!is.null(s)) set_rows[[length(set_rows)+1]] <- .kgSetRow(s, la$display, ph, P[[ph]]$cvn, P[[ph]]$pt, contrast,
                                                                      subset_note, method = .kgRelMethod)
          set_scores[[length(set_scores)+1]] <- .kgSetScores(s, la$display, ph)
        }
      } else if(do_set){
        sx <- FR$mat[intersect(rownames(FR$mat), read_vals), dcommon, drop = FALSE]
        for(ph in phenos){
          s <- .kgFrySet(sx, P[[ph]]$phv, P[[ph]]$cvd, P[[ph]]$pt, ct_lvl, ct_ref, minN)
          if(!is.null(s)) set_rows[[length(set_rows)+1]] <- .kgSetRow(s, la$display, ph, P[[ph]]$cvn, P[[ph]]$pt, contrast, subset_note)
        }
      }
    }

  } else {
    # ===================== LIMMA engine (ALWAYS a FULL-LAYER fit) =====================
    # Two callers: (a) features="all" -> report every feature; (b) a NAMED set >= the
    # threshold -> fit the FULL layer, then EXTRACT the named rows. Never subset-then-fit:
    # eBayes estimates its variance prior + mean-variance trend from the genes in the fit, so
    # fitting only a named (usually biased) subset would make each p-value depend on which
    # other features happened to be named. Fit all -> subset the RESULTS.
    work <- list()
    if(is_screen){
      # a staged spec ("<primary>/<fallback>" groups) -> the layers to run + the layers offered (see `.kgStagedLayers`)
      if(grepl("/", layers, fixed = TRUE)){
        st <- .kgStagedLayers(con, layers, keep_ids, phframe, phenos, minN)
        layers <- st$run; offered <- st$offered
      }
      # "all" on a SCREEN = the GENE layers (entity-coherent); a metabolite/contaminant/flux
      # screen is opt-in by naming its layer (e.g. layers="metabolite_hg").
      lays <- if(identical(layers, "all")) .kgGeneLayerKeys else tolower(trimws(strsplit(layers, "[,;]")[[1]]))
      for(ln in lays){
        lp <- .kgScreenLayer(ln); if(is.null(lp)) next
        if(identical(tolower(ln), "contaminant")){
          # one work item per tissue (`.kgContaminantTissues`); `mkey` keeps the KEY for the method choice
          for(ct in .kgContaminantTissues)
            work[[length(work)+1]] <- list(disp = ct$label, tbl = lp[1], id = lp[2], keep = NULL,
                                           ent = "contaminant", mkey = "contaminant", rowfilter = ct$rowfilter)
          next
        }
        work[[length(work)+1]] <- list(disp = ln, tbl = lp[1], id = lp[2], keep = NULL,
                                       ent = .kgLayerEntity(tolower(ln)))
      }
    } else {
      key <- vapply(feats, function(f) paste(f$layer$display, f$layer$table, f$layer$read_col, sep = ""), character(1))
      for(k in unique(key)){
        idx <- which(key == k); la <- feats[[idx[1]]]$layer
        work[[length(work)+1]] <- list(
          disp = la$display, tbl = la$table, id = la$read_col,
          keep = unique(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1))),
          ent  = feats[[idx[1]]]$etype, mkey = .kgLayerKey(la$display), rowfilter = la$rowfilter)
      }
    }
    for(w in work){
      L <- .kgLoadLayerMatrix(con, w$tbl, w$id, w$rowfilter); if(is.null(L)) next
      dcommon <- intersect(colnames(L$mat), keep_ids)
      .note_layer(w$disp, dcommon)
      if(length(dcommon) < minN) next
      Lmat <- L$mat[, dcommon, drop = FALSE]
      ent   <- w$ent
      small <- identical(covariates, "default") && (ent %in% .kgSmallNTypes)
      meth <- .kgAssocMethod(method, if(!is.null(w$mkey)) w$mkey else w$disp, is_screen)
      for(ph in phenos){
        pt  <- .kgPhenoType(m, ph)
        cvn <- if(small) character(0) else covcols_for(ph)
        md  <- .kgCoerce(m, unique(c(ph, cvn))); rownames(md) <- rec
        md  <- md[intersect(rownames(md), dcommon), , drop = FALSE]
        # narrowed on the fitted rows, for the same reason as the screen path above: `cvn` is what
        # the `Covariates` column reports, and the fit has already dropped anything constant here.
        cvn <- .kgVaryingCovs(md, cvn)
        if(do_set && !is.null(w$keep)){
          s <- if(identical(set_kind, "relative"))
                 .kgRelScoreSet(Lmat, w$keep, stats::setNames(md[[ph]], rownames(md)),
                                if(length(cvn)) md[, cvn, drop = FALSE] else NULL, pt, ct_lvl, ct_ref, minN)
               else
                 .kgFrySet(Lmat[intersect(rownames(Lmat), w$keep), , drop = FALSE],
                           stats::setNames(md[[ph]], rownames(md)), if(length(cvn)) md[, cvn, drop = FALSE] else NULL,
                           pt, ct_lvl, ct_ref, minN)
          if(!is.null(s)) set_rows[[length(set_rows)+1]] <- .kgSetRow(s, w$disp, ph, cvn, pt, contrast, subset_note,
                                                                      method = if(identical(set_kind, "relative")) .kgRelMethod else "fry")
          set_scores[[length(set_scores)+1]] <- .kgSetScores(s, w$disp, ph)
        }
        # ONE choice of fit; the two functions return the same columns (kg_common.R).
        fitfun <- if(identical(meth, "kendall")) .kgKendallFit else .kgLimmaFit
        if(pt == "disc"){
          pc <- parse_contrast()
          # omnibus needs a real reference level; default to the first level present
          rr <- pc$r
          # ⚠ CHECKED FOR EVERY DISCRETE FIT, NOT ONLY THE OMNIBUS. This guard used to sit inside
          # the `anova` branch, so a NAMED contrast ("Type1-None") walked straight into limma with a
          # one-level factor and crashed. The level count is a property of the DATA after subsetting
          # — it has nothing to do with which contrast was asked for.
          lv <- levels(factor(as.character(md[[ph]])))
          if(length(lv) < 2){ one.class <- TRUE; next }
          if(identical(pc$c, "anova") && !nzchar(rr)) rr <- lv[1]
          de <- if(identical(meth, "kendall"))
                  .kgKendallFit(Lmat, L$info, md, ph, pt, ref = rr, contrast = pc$c, minN = minN)
                else fitfun(Lmat, L$info, md, ph, pt, ref = rr, contrast = pc$c)
        } else {
          de <- if(identical(meth, "kendall"))
                  .kgKendallFit(Lmat, L$info, md, ph, pt, "", "", minN = minN)
                else fitfun(Lmat, L$info, md, ph, pt, "", "")
        }
        if(is.null(de) || nrow(de) == 0) next
        # EXTRACT the named rows AFTER the full-layer fit (prior/trend stay correct).
        if(!is.null(w$keep)){
          de <- de[as.character(de$FeatureId) %in% w$keep, , drop = FALSE]
          if(nrow(de) == 0) next
        }
        de <- data.frame(
          Feature = de$FeatureId, Symbol = de$Symbol, Entity = ent, Layer = w$disp, Phenotype = ph,
          N = de$N, Effect = signif(de$Effect, 3), EffectType = de$EffectType,
          Statistic = signif(de$T_statistic, 3), P_value = signif(de$P_value, 4),
          Direction = ifelse(is.na(de$Effect), "NA", ifelse(de$Effect > 0, "up", "down")),
          Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
          Comparison = if(pt == "disc") (if(nzchar(contrast)) contrast else "anova") else "continuous",
          Subset = subset_note, Method = meth,
          Adjusted_p_value = signif(de$Adjusted_p_value, 4), stringsAsFactors = FALSE)
        rows[[length(rows)+1]] <- de
      }
    }
  }

  # ===================== DNA METHYLATION (per CpG; ALWAYS the full-layer limma fit) =====================
  # The same fit the limma engine runs (`.kgLimmaFit`: trend + robust eBayes, the question's covariates, constant ones
  # dropped), on EVERY CpG of the key, then the named CpGs extracted -- never subset-then-fit, so the variance prior comes
  # from the whole layer. A named set's CpGs also get the group test (fry), exactly as on the other layers.
  for(mk in meth_keys){
    L <- .kgLoadMethylation(mk); if(is.null(L)) next
    keep <- if(is_screen) NULL else .kgMethRows(L, ftoks)
    dcommon <- intersect(colnames(L$mat), keep_ids)
    if(!is_screen && !length(keep)){ .note_layer(L$disp, dcommon, status = "no_feature"); next }
    .note_layer(L$disp, dcommon)
    if(length(dcommon) < minN) next
    Lmat <- L$mat[, dcommon, drop = FALSE]
    for(ph in phenos){
      pt  <- .kgPhenoType(m, ph)
      cvn <- covcols_for(ph)
      md  <- .kgCoerce(m, unique(c(ph, cvn))); rownames(md) <- rec
      md  <- md[intersect(rownames(md), dcommon), , drop = FALSE]
      cvn <- .kgVaryingCovs(md, cvn)
      if(do_set && identical(set_kind, "relative") && length(keep)){
        s <- .kgRelScoreSet(Lmat, ftoks, stats::setNames(md[[ph]], rownames(md)),
                            if(length(cvn)) md[, cvn, drop = FALSE] else NULL, pt, ct_lvl, ct_ref, minN, genes = L$genes)
        if(!is.null(s)) set_rows[[length(set_rows)+1]] <- .kgSetRow(s, L$disp, ph, cvn, pt, contrast, subset_note,
                                                                    method = .kgRelMethod)
        set_scores[[length(set_scores)+1]] <- .kgSetScores(s, L$disp, ph)
      } else if(do_set && length(keep)){
        s <- .kgFrySet(Lmat[keep, , drop = FALSE], stats::setNames(md[[ph]], rownames(md)),
                       if(length(cvn)) md[, cvn, drop = FALSE] else NULL, pt, ct_lvl, ct_ref, minN)
        if(!is.null(s)) set_rows[[length(set_rows)+1]] <- .kgSetRow(s, L$disp, ph, cvn, pt, contrast, subset_note)
      }
      if(pt == "disc"){
        pc <- parse_contrast(); rr <- pc$r
        lv <- levels(factor(as.character(md[[ph]])))
        if(length(lv) < 2){ one.class <- TRUE; next }
        if(identical(pc$c, "anova") && !nzchar(rr)) rr <- lv[1]
        de <- .kgLimmaFit(Lmat, L$info, md, ph, pt, ref = rr, contrast = pc$c)
      } else {
        de <- .kgLimmaFit(Lmat, L$info, md, ph, pt, "", "")
      }
      if(is.null(de) || nrow(de) == 0) next
      if(!is.null(keep)){
        de <- de[as.character(de$FeatureId) %in% keep, , drop = FALSE]
        if(nrow(de) == 0) next
      }
      # a SCREEN names its layer by KEY, as the gene-layer screens do (`disp = ln`): kgEnrichment matches the user folder's
      # differential on that key (`.kgEnrichFromUserFolder`); a targeted run keeps the display name
      rows[[length(rows)+1]] <- data.frame(
        Feature = de$FeatureId, Symbol = de$Symbol, Entity = "cpg", Layer = if(is_screen) mk else L$disp, Phenotype = ph,
        N = de$N, Effect = signif(de$Effect, 3), EffectType = de$EffectType,
        Statistic = signif(de$T_statistic, 3), P_value = signif(de$P_value, 4),
        Direction = ifelse(is.na(de$Effect), "NA", ifelse(de$Effect > 0, "up", "down")),
        Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
        Comparison = if(pt == "disc") (if(nzchar(contrast)) contrast else "anova") else "continuous",
        Subset = subset_note, Method = "limma",
        Adjusted_p_value = signif(de$Adjusted_p_value, 4), stringsAsFactors = FALSE)
    }
  }

  # an OFFERED layer (staged spec, not run, has donors) is written with its donor counts and Status `not_calculated`, so
  # the caller can say it can still be calculated; contaminants are named per tissue, as the run would name them
  for(k in offered){
    dc <- intersect(.kgLayerDonors(con, k), keep_ids)
    disps <- if(identical(k, "contaminant")) vapply(.kgContaminantTissues, function(ct) ct$label, character(1)) else k
    for(dsp in disps) .note_layer(dsp, dc, status = "not_calculated")
  }
  # the per-layer status (see `.note_layer`), written on EVERY return from here on -- RES-OK and the empty ones alike
  utils::write.csv(if(length(lstat)) do.call(rbind, lstat) else
                     data.frame(Layer = character(0), Phenotype = character(0), N_layer = integer(0),
                                N_usable = integer(0), MinN = integer(0), Status = character(0),
                                Subset = character(0), N_level = integer(0), N_ref = integer(0)),
                   "kg_assoc_layers.csv", row.names = FALSE)
  # the named-set test, when asked for: its OWN file (empty when no set row could be computed), so
  # kg_assoc.csv and the RES-OK count stay exactly as before. BH across this call's set rows (fdrFamily "query").
  if(do_set){
    so <- if(length(set_rows)) do.call(rbind, set_rows) else
            data.frame(Layer = character(0), Phenotype = character(0), N = integer(0), N_members = integer(0),
                       Direction = character(0), P_value = numeric(0), Covariates = character(0),
                       Comparison = character(0), Subset = character(0), Method = character(0),
                       Effect = numeric(0), N_layer = integer(0), N_predicted = integer(0))
    so$P_family <- if(nrow(so)) signif(stats::p.adjust(so$P_value, method = "BH"), 4) else numeric(0)
    utils::write.csv(so, "kg_assoc_set.csv", row.names = FALSE)
    # the relative test's per-donor scores (2026-10-19), its OWN file so nothing above changes; empty when none was computed
    if(identical(set_kind, "relative")){
      set_scores <- Filter(Negate(is.null), set_scores)
      sc <- if(length(set_scores)) do.call(rbind, set_scores) else
              data.frame(Layer = character(0), Phenotype = character(0), record_id = character(0), group = character(0),
                         score = numeric(0), predicted = logical(0))
      utils::write.csv(sc, "kg_assoc_setscore.csv", row.names = FALSE)
    }
  }

  # ⚠ "NOTHING TO COMPARE" IS NOT "COMPARED AND FOUND NOTHING". A bare RES-NO means the test RAN and
  # came back empty, which IS an answer; RES-ONE-CLASS means it could not run at all. Merging them
  # would teach a reader that a subset with one diagnosis level is a measured null.
  if(length(rows) == 0)
    return(if(one.class || isTRUE(.GlobalEnv$.kg_one_class)) "RES-ONE-CLASS" else "RES-NO")
  # one call can mix targeted rows (no Adjusted_p_value) with methylation rows (limma, with it): give every row block the
  # same columns, NA where a block has none. Blocks that already match are unchanged.
  .allc <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(d){ for(cn in setdiff(.allc, names(d))) d[[cn]] <- NA; d[, .allc, drop = FALSE] })
  out <- do.call(rbind, rows)

  # ---- family FDR (default #1: BH across exactly the tests asked; never genome-wide) ----
  if(identical(fdr_family, "query")){
    out$P_family <- signif(stats::p.adjust(out$P_value, method = "BH"), 4)
  } else if(identical(fdr_family, "genomewide") && "Adjusted_p_value" %in% names(out)){
    out$P_family <- out$Adjusted_p_value
  }

  out <- out[order(out$P_value), ]
  ts  <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 30))
  utils::write.csv(out, paste0("kg_assoc_", safe(features), "_", safe(phenotypes), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_assoc.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgAssoc <- function(...) .kgGuard("kgAssoc", .kgAssocImpl, ...)

# ---- SEVERAL NAMED SETS IN ONE CALL (2026-10-22, user: "go ahead with running several sets in one R call") -------------
# A question naming several sets ("Are Hippo and Wnt pathway promoters less methylated than other gene promoters in Cluster 2
# compared with Cluster 1 …?") was answered by one kgAssoc call per set, and each call paid the same fixed cost again: the
# scripts sourced (~5 s), the layer read (methylation .qs, 2.7 s) and the full-layer fit (limma, 1.6 s per phenotype) --
# MEASURED 2026-10-19: 14.4 s per call at load 5.6, 71-76 s each under load. Here every set runs the SAME guarded call
# (`kgAssoc`) with its own members, in ITS OWN user folder (a sibling of this one, made by /login), so each set's files,
# rows, set test, adjustments and RES code are exactly those of its own call. Only what does not depend on the members is
# shared, through a cache that lives for this call only: the layer reads and the full-layer fits.
# `sets` = "<user folder>=<members>||<user folder>=<members>…" (members as `features` takes them; the folders must exist).
# Returns "RES-SETS;<set 1's code>||<set 2's code>…", one code per set, in the order given.
.kgMemoKey <- function(a){
  part <- function(x){
    if(is.null(x)) return("NULL")
    if(is.matrix(x) || is.data.frame(x)){
      rn <- rownames(x); cn <- colnames(x)
      return(paste0("M", nrow(x), "x", ncol(x), "|", paste(utils::head(rn, 3), collapse = ","), "|",
                    paste(utils::tail(rn, 3), collapse = ","), "|", paste(cn, collapse = ",")))
    }
    if(is.list(x)) return(paste0("L(", paste(names(x), vapply(x, part, character(1)), sep = "=", collapse = ";"), ")"))
    if(is.atomic(x)) return(if(length(x) <= 64) paste0("A", paste(as.character(x), collapse = ","))
                            else paste0("A", length(x), "|", paste(utils::head(x, 3), collapse = ","), "|",
                                        paste(utils::tail(x, 3), collapse = ",")))
    NA_character_                                   # a DB connection: not part of what is computed
  }
  k <- vapply(a, part, character(1)); k <- k[!is.na(k)]
  paste(names(a), k, sep = "=", collapse = "&")
}
.kgAssocSetsImpl <- function(features = "", phenotypes = "", layers = "all", covariates = "default", subset = "all",
                             contrast = "", ref = "", fdr_family = "query", minN = "10", method = "auto", mode = "tool",
                             set_test = "", sets = ""){
  parts <- strsplit(as.character(sets), "||", fixed = TRUE)[[1]]
  parts <- trimws(parts[nzchar(trimws(parts))])
  if(!length(parts) || !all(grepl("=", parts, fixed = TRUE))) return("RES-NO-FEATURE")
  dirs <- trimws(sub("=.*$", "", parts)); mems <- trimws(sub("^[^=]*=", "", parts))
  wd0 <- getwd(); root <- dirname(wd0)
  if(any(!grepl("^[A-Za-z0-9_]+$", dirs)) || any(!dir.exists(file.path(root, dirs))) || any(!nzchar(mems)))
    return("RES-NO-FEATURE")
  memo <- c(".kgLoadLayerMatrix", ".kgLoadMethylation", ".kgLimmaFit", ".kgKendallFit")
  memo <- memo[vapply(memo, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)]
  orig <- mget(memo, envir = .GlobalEnv)
  on.exit({ setwd(wd0); for(n in names(orig)) assign(n, orig[[n]], envir = .GlobalEnv) }, add = TRUE)
  cache <- new.env(parent = emptyenv())
  for(n in memo) local({
    f0 <- orig[[n]]; nm <- n
    assign(nm, function(...){
      a <- list(...); k <- paste0(nm, "#", .kgMemoKey(a))
      if(exists(k, envir = cache, inherits = FALSE)) return(get(k, envir = cache, inherits = FALSE))
      v <- f0(...); assign(k, v, envir = cache); v
    }, envir = .GlobalEnv)
  })
  codes <- character(length(parts))
  for(i in seq_along(parts)){
    setwd(file.path(root, dirs[i]))
    codes[i] <- as.character(kgAssoc(features = mems[i], phenotypes = phenotypes, layers = layers, covariates = covariates,
                                     subset = subset, contrast = contrast, ref = ref, fdr_family = fdr_family, minN = minN,
                                     method = method, mode = mode, set_test = set_test))[1]
  }
  paste0("RES-SETS;", paste(codes, collapse = "||"))
}
kgAssocSets <- function(...) .kgGuard("kgAssocSets", .kgAssocSetsImpl, ...)
