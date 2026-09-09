# ==============================================================================
# kg_functions/kg_enrichment.R  --  pathway ENRICHMENT (op#5).
# Frames: F12 (disease pathways) - F27 (mechanism) - F4/F7/F14 (enrichment after a set).
#
# ⚠ THIS FUNCTION NEVER COMPUTES A DIFFERENTIAL -- IT CONSUMES ONE THAT ALREADY EXISTS.
# (user ruling 2026-07-20). It is the same pattern the deployed performGSEA follows: that
# function reads dea_results.csv -- the file DonorRegression already wrote -- and refuses
# with "run the regression first" when it is not there. Enrichment is downstream by nature;
# re-fitting a differential inside it duplicates work and can silently disagree with the
# numbers the user was just shown.
#
# TWO SOURCES, resolved by `input`:
#   1. USER FOLDER   kg_assoc.csv     -- a previous kgAssoc SCREEN run (custom cov / subset)
#   2. PRECOMPUTE    omics_outcomes   -- HI_precomputed_v2.sqlite, the default spec
# Neither holds it -> RES-NO-DIFFERENTIAL, and the CALLER runs kgAssoc first. The frames
# already emit A(1) then A(5) as two separate calls, so the chaining belongs to the caller.
#
# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------
# phenotype     A metadata_sum_norm column (diagnosis, hba1c, donorage, ...).
#               MATCHING KEY ONLY -- nothing is fitted here, so this selects WHICH stored
#               differential to read, it is never a model term.
#
# contrast      ""      -> continuous phenotype
#               "A-B"   -> that group contrast (A vs B)
#               "anova" -> omnibus
#               MATCHING KEY. ⚠ An omnibus carries NO signed effect, so analysisType is
#               forced to "ora" -- there is nothing to rank. (Same rule as before; it is
#               §K's method-by-INPUT, not a preference.)
#
# ref           Reference level. The "A-B" contrast string already carries it, so `ref` is
#               only consulted when `contrast` is not in A-B form. Kept for call
#               compatibility.
#
# omics_layer   rnaseq | nanostring | protein | pbrna_alpha | pbrna_beta
#               | metabolite_hg | metabolite_lg | metabolite_ratio
#               MATCHING KEY, and it also selects the annotation table that supplies the
#               enrichment identifier (see "THE ENRICHMENT KEY" below).
#               ⚠ protein and metabolite resolve to the _combat tables -- per the user
#               ruling (2026-07-20) those are the only correct ones for ALL downstream
#               analysis, not merely for pathways.
#
# library       gene layers:       kegg | reactome | go_bp | go_mf
#               metabolite layers: hsa_kegg  -- THE ONLY ONE (user ruling 2026-07-20; §P says
#                                  the same: metabolite sets are KEGG-only). It loads
#                                  kegg_hsa_met.qs; hsaGem_metabolite.* is NOT reachable.
#               ⚠ ONE LIBRARY PER CALL, TESTED ALONE. §P: every library has its OWN
#               background and its OWN multiple-testing correction and is NEVER pooled --
#               the caller loops over libraries, and the result NAMES which one it is
#               (the Library column below).
#               ⚠ go_cc is deliberately NOT offered: CellularComponent describes where a
#               protein sits, not what it does, and §P excludes it from enrichment.
#
# analysisType  "gsea" -> ranked test; needs a signed effect
#               "ora"  -> over-representation of the significant set against the universe
#               Forced to "ora" for an omnibus contrast.
#
# covariates    MATCHING KEY ONLY -- nothing is adjusted here; the differential was already
#               fitted with whatever covariates it used. "default" means the standard 5
#               minus the tested variable.
#               Used to accept or reject a USER-FOLDER result.
#               ⚠ The PRECOMPUTE cannot be checked this way -- omics_outcomes carries no
#               covariates column -- so on that path the documented default set is trusted,
#               and the output records Source="precompute" so the answer can say so.
#
# fdr           Significance threshold. Used for the reported count and, for ORA, to pick
#               the significant set. It never filters the GSEA input (§K: GSEA needs no
#               significance cutoff, and the whole ranked universe is what it consumes).
#
# input         auto (DEFAULT) | userfolder | precompute
#
#               auto        USER FOLDER FIRST, then the precompute. Rationale: a user who
#                           has just run a custom limma is asking about THEIR result; only
#                           when no matching prior result exists do we fall back to the
#                           default-spec precompute.
#               userfolder  Force the prior result. Typed refusal if it is absent, does not
#                           match the question, or is not a valid universe.
#               precompute  Force the default-spec result.
#
#               ⚠ A USER-FOLDER RESULT MUST BE A FULL-LAYER SCREEN (Method == "limma").
#               A targeted kgAssoc run writes only the named features' rows. Enriching a
#               handful of named genes against a whole library is a meaningless test that
#               still returns confident-looking numbers, so it is refused with
#               RES-NO-UNIVERSE rather than run. (kg_assoc.R applies the same principle from
#               the other side: it always fits the FULL layer and then extracts the named
#               rows, because eBayes borrows its variance prior across features.)
#
# ------------------------------------------------------------------------------
# THE ENRICHMENT KEY -- read from the omics table, not from the differential
# ------------------------------------------------------------------------------
# Neither source carries the identifier a library is keyed on: kg_assoc.csv holds
# Feature/Symbol, and omics_outcomes holds Feature/ID. So the key is looked up from the
# layer's own annotation columns (a light 2-3 column read, no donor matrix):
#     gene layers        -> gene_id   (entrez)     for kegg / reactome / go_bp / go_mf
#     metabolite layers  -> kegg_id                for hsa_kegg
# Verified present on proc_metabolite_combat_* (compound, hmdb_id, kegg_id, gem_id, ...),
# so metabolite enrichment needs NO external id mapping.
#
# RANKING  GSEA ranks by Effect (the coefficient / log2FC), matching the deployed
#          performGSEA default rank.stat = "coef" (user ruling 2026-07-20). Both stores
#          also carry the t-statistic, so the alternative remains available without a refit.
#
# Result:  kg_enrichment.csv + a timestamped copy. Returns "RES-OK;<n_sig>".
# Refusals: RES-NO-OMICS | RES-NO-LIBRARY | RES-NO-DIFFERENTIAL | RES-NO-UNIVERSE
#         | RES-NO-KEY | RES-NO-ENRICH
# ==============================================================================

.kgEnrichDefaultCovs <- c("donorage", "donorsex", "bodymassindex",
                          "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# layer key -> precompute Omics_ID + annotation table + id column + entity class.
# ⚠ protein and metabolite are the _combat tables (user ruling 2026-07-20).
#
# ⚠ `omics_id` AND `tbl` ARE TWO DIFFERENT KEYS AND USED TO BE ONE FIELD. `omics_id` is what
# HI_precomputed_v2.omics_outcomes stores in `Omics_ID`; `tbl` is the annotation table to read in
# HI_omics_v2.sqlite. For rnaseq and pbrna they happen to be the same string, which is why the
# single field survived — but for protein the id is `proc_prot_v2` while the table is
# `proc_prot_combat` (the id retrieves the combat table), and every metabolite layer is stored as
# `proc_metabolite_*` while its annotation lives in `proc_metabolite_combat_*`. Using one field for
# both made the precompute lookup miss and the layer answered RES-NO-DIFFERENTIAL for every
# phenotype. MEASURED 2026-08-08 on hba1c: Omics_ID proc_prot_v2 = 7,822 rows, proc_prot_combat = 0.
#
# ⚠ `id` MUST NAME THE COLUMN THAT HOLDS THE PRECOMPUTE'S OWN `ID` VALUES — it is the merge key, and
# a wrong column empties the join and answers RES-NO-KEY. MEASURED on hba1c:
#   protein     ID = "A1CF"                        -> `symbol`   (NOT `id`, which is UniProt Q9NQ94)
#   nanostring  ID = "ABCC8"                       -> `symbol`   (NOT `gene_id`, which is entrez 6833)
#   metabolite  ID = "AGPKZVBTJJNPAG-WHFBIAKZSA-N" -> `inchikey` (NOT `compound`, "Isoleucine")
#   rnaseq      ID = "ENSG00000000003"             -> `accession` (unchanged, already correct)
#   pbrna_*     ID = "A1BG"                        -> `symbol`    (unchanged, already correct)
.kgEnrichLayer <- function(omics_layer){
  switch(tolower(omics_layer),
    "rnaseq"           = list(omics_id="proc_rnaseq",           tbl="proc_rnaseq",                  id="accession", entity="gene"),
    "nanostring"       = list(omics_id="proc_nanostring_merge", tbl="proc_nanostring_merge",        id="symbol",    entity="gene"),
    "protein"          = list(omics_id="proc_prot_v2",          tbl="proc_prot_combat",             id="symbol",    entity="gene"),
    "pbrna_alpha"      = list(omics_id="proc_pbrna_Alpha",      tbl="proc_pbrna_Alpha",             id="symbol",    entity="gene"),
    "pbrna_beta"       = list(omics_id="proc_pbrna_Beta",       tbl="proc_pbrna_Beta",              id="symbol",    entity="gene"),
    "metabolite_hg"    = list(omics_id="proc_metabolite_HG",    tbl="proc_metabolite_combat_HG",    id="inchikey",  entity="metabolite"),
    "metabolite_lg"    = list(omics_id="proc_metabolite_LG",    tbl="proc_metabolite_combat_LG",    id="inchikey",  entity="metabolite"),
    "metabolite_ratio" = list(omics_id="proc_metabolite_ratio", tbl="proc_metabolite_combat_ratio", id="inchikey",  entity="metabolite"),
    NULL)
}

# library -> the annotation column it is keyed on, and the entity class it belongs to.
# ⚠ hsa_kegg (file kegg_hsa_met.qs) is the ONLY metabolite library -- user ruling 2026-07-20,
# and §P says the same: metabolite sets are KEGG-only.
.kgEnrichKey <- function(library){
  switch(tolower(library),
    "hsa_kegg" = list(col="kegg_id", entity="metabolite"),
    list(col="gene_id", entity="gene"))     # kegg / reactome / go_bp / go_mf
}

# FeatureId -> (Symbol, Key) from the layer's annotation columns. Reads 2-3 columns of the
# proc_* table -- never the donor matrix.
.kgEnrichAnnot <- function(con, tbl, id_col, key_col){
  have <- try(DBI::dbListFields(con, tbl), silent = TRUE)
  if(inherits(have, "try-error")) return(NULL)
  if(!(id_col %in% have) || !(key_col %in% have)) return(NULL)
  cols <- unique(c(id_col, key_col, if("symbol" %in% have) "symbol"))
  q <- sprintf('SELECT %s FROM %s', paste(sprintf('"%s"', cols), collapse = ", "), tbl)
  a <- try(DBI::dbGetQuery(con, q), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) == 0) return(NULL)
  data.frame(FeatureId = as.character(a[[id_col]]),
             Key       = as.character(a[[key_col]]),
             Symbol    = if("symbol" %in% names(a)) as.character(a$symbol) else as.character(a[[id_col]]),
             stringsAsFactors = FALSE)
}

# are two covariate specifications the same SET? (order and spacing are irrelevant)
.kgCovSetEq <- function(written, expected){
  w <- trimws(strsplit(ifelse(is.na(written) | written == "", "none", written), "[,;]")[[1]])
  w <- sort(unique(tolower(w[nzchar(w)])))
  e <- sort(unique(tolower(expected[nzchar(expected)])))
  if(length(e) == 0) e <- "none"
  if(length(w) == 0) w <- "none"
  identical(w, e)
}

# ---- source 1: the USER FOLDER (a previous kgAssoc SCREEN run) ----------------
# Returns data.frame(FeatureId, Effect, P_value, Adjusted_p_value) or a refusal string.
.kgEnrichFromUserFolder <- function(phenotype, layer_key, contrast, cov_expected){
  f <- "kg_assoc.csv"
  if(!file.exists(f)) return(NULL)
  d <- try(utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if(inherits(d, "try-error") || nrow(d) == 0) return(NULL)
  need <- c("Feature", "Layer", "Phenotype", "Effect", "P_value", "Method", "Comparison")
  if(!all(need %in% names(d))) return(NULL)

  # match the QUESTION first -- this is what stops a stale result from an earlier,
  # different question in the same folder being enriched by accident.
  d <- d[.kgNorm(d$Phenotype) == .kgNorm(phenotype), , drop = FALSE]
  if(nrow(d) == 0) return(NULL)
  d <- d[.kgNorm(d$Layer) == .kgNorm(layer_key), , drop = FALSE]
  if(nrow(d) == 0) return(NULL)
  want.cmp <- if(!nzchar(contrast)) "continuous" else contrast
  d <- d[.kgNorm(d$Comparison) == .kgNorm(want.cmp), , drop = FALSE]
  if(nrow(d) == 0) return(NULL)
  if("Covariates" %in% names(d)){
    d <- d[vapply(d$Covariates, .kgCovSetEq, logical(1), expected = cov_expected), , drop = FALSE]
    if(nrow(d) == 0) return(NULL)
  }

  # ⚠ it must be a FULL-LAYER SCREEN. A targeted lm result is not a universe.
  if(!any(.kgNorm(d$Method) == "limma")) return("RES-NO-UNIVERSE")
  d <- d[.kgNorm(d$Method) == "limma", , drop = FALSE]

  adj <- if("Adjusted_p_value" %in% names(d)) d$Adjusted_p_value
         else if("P_family" %in% names(d))    d$P_family
         else rep(NA_real_, nrow(d))
  data.frame(FeatureId = as.character(d$Feature),
             Effect    = suppressWarnings(as.numeric(d$Effect)),
             P_value   = suppressWarnings(as.numeric(d$P_value)),
             Adjusted_p_value = suppressWarnings(as.numeric(adj)),
             stringsAsFactors = FALSE)
}

# ---- source 2: the PRECOMPUTE (omics_outcomes) --------------------------------
# Anchored on Metadata_ID_col (index ix_meta) -- an indexed range scan, never a table scan.
.kgEnrichFromPrecompute <- function(phenotype, omics_id, contrast){
  .kgSetPaths()
  db <- paste0(sqlite.path, "HI_precomputed_v2.sqlite")
  if(!file.exists(db)) return(NULL)
  con <- try(DBI::dbConnect(RSQLite::SQLite(), db), silent = TRUE)
  if(inherits(con, "try-error")) return(NULL)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
  # sqlite stores '' for a continuous contrast (neo4j uses 'NA' -- the two stores encode
  # this differently, so normalise here rather than at the caller).
  ct <- if(identical(contrast, "anova")) "anova" else contrast
  q  <- paste('SELECT ID AS FeatureId, Coefficient AS Effect,',
              'P__value AS P_value, Adjusted_p__value AS Adjusted_p_value',
              'FROM omics_outcomes',
              'WHERE Metadata_ID_col = ? AND Omics_ID = ? AND COALESCE(Contrast, "") = ?')
  d <- try(DBI::dbGetQuery(con, q, params = list(phenotype, omics_id, ct)), silent = TRUE)
  if(inherits(d, "try-error") || nrow(d) == 0) return(NULL)
  d$FeatureId <- as.character(d$FeatureId)
  d
}

.kgEnrichmentImpl <- function(phenotype, contrast = "", ref = "",
                         omics_layer = "rnaseq", library = "kegg",
                         analysisType = "gsea", covariates = "default",
                         fdr = "0.05", input = "auto", mode = "tool"){
  library(RSQLite); library(DBI)
  .kgSetPaths()
  fdr.n <- suppressWarnings(as.numeric(fdr)); if(is.na(fdr.n)) fdr.n <- 0.05
  input <- tolower(trimws(ifelse(is.null(input) || !nzchar(input), "auto", input)))

  lp <- .kgEnrichLayer(omics_layer); if(is.null(lp)) return("RES-NO-OMICS")
  kp <- .kgEnrichKey(library)

  # a gene library on a metabolite layer (or the reverse) can only produce nonsense
  if(!identical(kp$entity, lp$entity)) return("RES-NO-LIBRARY")

  # normalise the contrast the same way the stored results do
  if(nzchar(contrast) && !identical(contrast, "anova")){
    pp <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if(length(pp) != 2 && nzchar(ref)) contrast <- paste0(contrast, "-", ref)
  }
  # an omnibus has no signed effect -- nothing to rank, so ORA is the only valid test
  if(identical(contrast, "anova")) analysisType <- "ora"

  # the covariate SET this question implies (matching key only -- nothing is fitted here)
  cov.expected <- if(covariates %in% c("none", "NA", "")) character(0)
                  else if(identical(covariates, "default")) .kgEnrichDefaultCovs[.kgEnrichDefaultCovs != phenotype]
                  else { cv <- trimws(strsplit(covariates, "[,;]")[[1]]); cv[nzchar(cv) & cv != phenotype] }

  # ---- get the differential; NEVER compute one ----
  de <- NULL; src <- NA_character_
  if(input %in% c("auto", "userfolder")){
    uf <- .kgEnrichFromUserFolder(phenotype, omics_layer, contrast, cov.expected)
    if(is.character(uf)) return(uf)                 # RES-NO-UNIVERSE
    if(!is.null(uf)){ de <- uf; src <- "userfolder" }
  }
  if(is.null(de) && input %in% c("auto", "precompute")){
    # ⚠ THE PRECOMPUTE IS ONLY VALID FOR THE DEFAULT COVARIATE SET.
    # It is fitted with the standard covariates (minus the tested variable). Serving it for a
    # question that asked for a CUSTOM set -- or for none -- would answer with numbers adjusted
    # differently from what was requested, and nothing downstream could tell: the enrichment
    # would look complete and be wrong. So a custom set falls through to RES-NO-DIFFERENTIAL and
    # the CALLER runs kgAssoc with those covariates first. This holds even when `precompute` was
    # forced, because a covariate mismatch is a correctness problem, not a preference.
    cov.is.default <- .kgCovSetEq(paste(cov.expected, collapse = ","),
                                  .kgEnrichDefaultCovs[.kgEnrichDefaultCovs != phenotype])
    if(cov.is.default){
      pc <- .kgEnrichFromPrecompute(phenotype, lp$omics_id, contrast)
      if(!is.null(pc)){ de <- pc; src <- "precompute" }
    }
  }
  if(is.null(de) || nrow(de) == 0) return("RES-NO-DIFFERENTIAL")

  # ---- attach the library's identifier from the layer's annotation columns ----
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  ann <- .kgEnrichAnnot(con, lp$tbl, lp$id, kp$col)
  dbDisconnect(con)
  if(is.null(ann)) return("RES-NO-KEY")

  de <- merge(de, ann, by = "FeatureId", all.x = TRUE)
  de <- de[!is.na(de$Key) & de$Key != "" & de$Key != "NA", , drop = FALSE]
  if(nrow(de) == 0) return("RES-NO-KEY")

  # ---- enrichment ----
  if(tolower(analysisType) == "ora"){
    sig <- de$Key[!is.na(de$Adjusted_p_value) & de$Adjusted_p_value < fdr.n]
    er  <- .kgOra(sig, de$Key, library, fdr.n)
  } else {
    d2 <- de[!is.na(de$Effect), , drop = FALSE]
    # one row per identifier, keeping the largest |effect|
    d2 <- d2[order(-abs(d2$Effect)), , drop = FALSE]
    d2 <- d2[!duplicated(d2$Key), , drop = FALSE]
    ranks <- setNames(as.numeric(d2$Effect), d2$Key)
    # the Key -> Symbol map the annotation join already produced, so the leading edge comes back
    # as gene symbols rather than the library's internal ids
    .sym <- setNames(as.character(d2$Symbol), as.character(d2$Key))
    er <- .kgFgsea(ranks, library, fdr.n, sym = .sym)
  }
  if(is.null(er) || nrow(er) == 0) return("RES-NO-ENRICH")

  # §P: the result must NAME its library, and it carries where the differential came from
  er$Library <- library
  er$Source  <- src
  er <- er[order(er$Adjusted_p_value), , drop = FALSE]
  n.sig <- sum(er$Adjusted_p_value < fdr.n, na.rm = TRUE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(er, paste0("kg_enrichment_", safe(phenotype), "_", safe(library), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(er, "kg_enrichment.csv", row.names = FALSE)
  paste0("RES-OK;", n.sig)
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgEnrichment <- function(...) .kgGuard("kgEnrichment", .kgEnrichmentImpl, ...)
