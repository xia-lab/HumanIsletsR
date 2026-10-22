# ==============================================================================
# build_kg_resolver.R  --  build the Knowledge Search FEATURE resolver (qs)
# ------------------------------------------------------------------------------
# Replaces the stale, 9.2 MB knowledge_search_index_v2.json (pre proteomics-v2)
# with one small, fast qs built DIRECTLY from the current proc_* omics tables, so
# the resolver is always fresh and loads in ~0.2s instead of ~3.9s (fromJSON).
#
# Produces (in other.tables.path/display_kg/):
#   kg_feature_resolver.qs  - one row per (token, table). Columns:
#       token     normalized alias a user query may carry (symbol / accession(ENSG) /
#                 gene_id(entrez) / ensembl / UniProt id / Protein_Group / compound /
#                 kegg / hmdb / inchikey / ik_conn / gem / rxn)
#       etype     gene | metabolite | contaminant | reaction  -- disambiguates
#                 cross-type token collisions (e.g. gene INS vs metabolite Inosine)
#       layer     omics display name
#       tbl       proc_* table
#       read_col / read_val   the UNIQUE key to SELECT the feature row back from tbl
#
# PHENOTYPES are NOT built here: everything kgPhenoCorr needs (the phenotype list,
# the per-donor values, and cont/categorical inferable from the data) already lives
# in display_data/metadata_sum_norm.csv, so it is resolved live -- no qs, and no
# dependency on the stale proc_variable_summary.
#
# Donor namespace: every analysis proc_* table's joinable donors are R### (bulk
# tables are all-R###; pbrna is the R### subset of its mixed columns), so the
# feature reader just selects ^R[0-9]+$ columns -- no per-table donor flag.
#
# Identifier semantics VERIFIED against current data + processing scripts
# (2_proc_proteomics_v2.R, 2_process_pseudobulk.R):
#   - proteomics (proc_prot_combat) is SYMBOL/UniProt keyed; gene_id(entrez) has 76
#     NA rows -> match symbol + id + Protein_Group + entrez(if present), read by the
#     unique `id`; composite "GATD3A;GATD3B" / "A0A0B4J2D5;P0DPI2" split into tokens.
#   - pseudobulk (proc_pbrna_*) is SYMBOL keyed (rownames=input_symbol); gene_id is
#     entrez, ensembl is ENSG, both annotations that can be NA -> read by symbol.
#   - RNA-seq stores ENSG in `accession` (gene_id there is entrez).
#
# RUN (from anywhere):
#   Rscript /Users/lzy/NetBeansProjects/restxialab/src/main/webapp/resources/rscripts/HumanIsletsR/kg_functions/build_kg_resolver.R
# ==============================================================================

suppressMessages({ library(RSQLite); library(DBI); library(qs) })

# ---- locate + source filepath_utils.R relative to THIS script (deploy-safe) --
.args <- commandArgs(trailingOnly = FALSE)
.sp   <- sub("^--file=", "", .args[grep("^--file=", .args)])
.dir  <- if (length(.sp)) dirname(normalizePath(.sp)) else getwd()
source(file.path(.dir, "..", "web_tool_functions", "filepath_utils.R"))
setPaths()

norm    <- function(x) tolower(trimws(as.character(x)))
out.dir <- paste0(other.tables.path, "display_kg/")
db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")

# read_col = the UNIQUE column used to SELECT the feature row.
# alias    = every column whose value a user query may carry (all point to read_val).
cfg <- list(
  list(tbl="proc_rnaseq",                 layer="Bulk gene expression (RNA-seq)",    etype="gene",
       read_col="accession", alias=c("symbol","accession","gene_id")),
  list(tbl="proc_nanostring_merge",       layer="Bulk gene expression (Nanostring)", etype="gene",
       read_col="gene_id",   alias=c("symbol","gene_id")),
  list(tbl="proc_prot_combat",            layer="Bulk protein expression",           etype="gene",
       read_col="id",        alias=c("symbol","id","Protein_Group","gene_id")),
  list(tbl="proc_pbrna_Alpha",  layer="Pseudobulk gene expression (alpha cells)", etype="gene",
       read_col="symbol",    alias=c("symbol","gene_id","ensembl")),
  list(tbl="proc_pbrna_Beta",   layer="Pseudobulk gene expression (beta cells)",  etype="gene",
       read_col="symbol",    alias=c("symbol","gene_id","ensembl")),
  list(tbl="proc_pbrna_Delta",  layer="Pseudobulk gene expression (delta cells)", etype="gene",
       read_col="symbol",    alias=c("symbol","gene_id","ensembl")),
  list(tbl="proc_pbrna_PP",     layer="Pseudobulk gene expression (PP cells)",    etype="gene",
       read_col="symbol",    alias=c("symbol","gene_id","ensembl")),
  list(tbl="proc_metabolite_combat_HG",    layer="Metabolite (HG)",    etype="metabolite",
       read_col="compound", alias=c("compound","kegg_id","hmdb_id","inchikey","ik_conn","gem_id")),
  list(tbl="proc_metabolite_combat_LG",    layer="Metabolite (LG)",    etype="metabolite",
       read_col="compound", alias=c("compound","kegg_id","hmdb_id","inchikey","ik_conn","gem_id")),
  list(tbl="proc_metabolite_combat_ratio", layer="Metabolite (ratio)", etype="metabolite",
       read_col="compound", alias=c("compound","kegg_id","hmdb_id","inchikey","ik_conn","gem_id")),
  list(tbl="proc_contaminants", layer="Environmental contaminants", etype="contaminant",
       read_col="Compound", alias=c("Compound")),
  list(tbl="proc_flux_HG",    layer="Metabolic flux (HG)",    etype="reaction", read_col="rxn", alias=c("rxn")),
  list(tbl="proc_flux_LG",    layer="Metabolic flux (LG)",    etype="reaction", read_col="rxn", alias=c("rxn")),
  list(tbl="proc_flux_ratio", layer="Metabolic flux (ratio)", etype="reaction", read_col="rxn", alias=c("rxn"))
)

# alias columns that can carry composite ";"-joined ids -> split into separate tokens
split.cols <- c("symbol", "Protein_Group")

con   <- dbConnect(SQLite(), db.path)
parts <- list()
for (c in cfg) {
  cols <- unique(c(c$read_col, c$alias))
  df   <- dbGetQuery(con, sprintf("SELECT %s FROM %s",
                                   paste(sprintf('"%s"', cols), collapse=","), c$tbl))
  rv <- as.character(df[[c$read_col]])
  for (a in c$alias) {                          # vectorized per alias column
    av   <- as.character(df[[a]])
    keep <- !is.na(av) & nzchar(trimws(av)) & !is.na(rv)
    if (!any(keep)) next
    avk <- av[keep]; rvk <- rv[keep]
    if (a %in% split.cols && any(grepl(";", avk, fixed=TRUE))) {
      sp  <- strsplit(avk, ";", fixed=TRUE)
      tok <- norm(unlist(sp)); rvv <- rep(rvk, lengths(sp))
    } else {
      tok <- norm(avk); rvv <- rvk
    }
    ok <- nzchar(tok)
    parts[[length(parts)+1]] <- data.frame(token=tok[ok], etype=c$etype, layer=c$layer,
                                           tbl=c$tbl, read_col=c$read_col, read_val=rvv[ok],
                                           stringsAsFactors=FALSE)
  }
}
dbDisconnect(con)

feat.df <- do.call(rbind, parts)
feat.df <- feat.df[!duplicated(feat.df[, c("token","tbl")]), ]

qs::qsave(feat.df, paste0(out.dir, "kg_feature_resolver.qs"), preset="high")
cat(sprintf("kg_feature_resolver.qs : %d rows, %d tokens, %.2f MB  (%s)\n",
            nrow(feat.df), length(unique(feat.df$token)),
            file.info(paste0(out.dir,"kg_feature_resolver.qs"))$size/1e6,
            format(Sys.time())))
