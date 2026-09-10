# ==============================================================================
# HumanIslets web-tool data processing   
#   * De-duplicates the 3x seahorse transpose block and the 3x perifusion blocks.
#   * Seahorse comes from REDCap. Alice (ocr_calculations) is used ONLY for the
#     3 columns REDCap lacks but the UI shows: calc_atp_resp_gluc, calc_bhi,
#     calc_bhi_gluc. Everything else = REDCap.
#   * metadata_sum_norm / _raw are built IN MEMORY and each written ONCE
#     (no write-then-read-back-then-rewrite).
#   * Every output file is written exactly once.
#
# NOT YET RUN against live REDCap. Diff the outputs against the current
# display_data/*.csv before deploying.
#
# Paths are explicit arguments (were global before) for reproducibility.
# ==============================================================================
other.tables.path <- "/Users/lzy/humanislet/humanislets/"
sqlite.path ="/Users/lzy/humanislet/sqlite/"
alice.path   = "/Users/lzy/humanislet/ocr_calculations_2026-04-02.csv"
states.path  = "/Users/lzy/humanislet/humanislets/donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/donor_states.csv"
demo.path    = "/Users/lzy/humanislet/humanislets/updated_data/ADIIsletCoreHumanIsl-AgeSexBMIHbA1CDiabet_DATA_2026-08-01_2009.csv"
## ---- helpers ----------------------------------------------------------------

# Wide seahorse csv (variables in rows, "R###_<rep>" columns) -> tidy per-sample.
read_seahorse_wide <- function(path) {
  x <- read.csv(path, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(x)[1] <- "variable"
  x <- x[x$variable != "" & !is.na(x$variable), , drop = FALSE]
  rownames(x) <- x$variable
  x <- as.data.frame(t(x[, -1, drop = FALSE]), stringsAsFactors = FALSE)
  x$record_id <- substr(rownames(x), 1, 4)
  x$replicate <- substr(rownames(x), 6, 6)
  rownames(x) <- NULL
  meas <- setdiff(colnames(x), c("record_id", "replicate"))
  x[meas] <- lapply(x[meas], function(v) {
    v[v %in% c("No DNA", "", " ")] <- NA
    suppressWarnings(as.numeric(v))
  })
  x[, c("record_id", "replicate", meas)]
}

# One perifusion pzfx pair (non-diabetic + T2D tables) -> raw wide traces.
process_peri_trace <- function(nd_tbl, t2d_tbl, drop_na_time = FALSE) {
  if (drop_na_time) {
    nd_tbl  <- nd_tbl[!is.na(nd_tbl$`Time (min)`), ]
    t2d_tbl <- t2d_tbl[!is.na(t2d_tbl$`Time (min)`), ]
  }
  w <- cbind(nd_tbl, t2d_tbl[, -1])
  rownames(w) <- paste0("time_", w$`Time (min)`)
  w <- as.data.frame(t(w[, -1]))
  meta <- data.frame(record_id = trimws(gsub("_.*", "", rownames(w))),
                     replicate = gsub("R...._", "", gsub("R..._", "", rownames(w))))
  w <- cbind(meta, w)
  rownames(w) <- NULL
  w[order(w$record_id), ]
}

# Raw wide trace -> long, replicate-averaged, baseline-normalised (baseline = time<25).
normalize_peri_trace <- function(wide_df) {
  long <- reshape2::melt(wide_df, id.vars = c("record_id", "replicate"))
  long <- na.omit(long)
  long <- as.data.table(long)
  long <- long[, .(insulin = mean(value)), by = .(record_id, variable)]
  long[, time := as.numeric(gsub("time_", "", variable))]
  long <- long[, .(record_id, insulin, time)]
  bl <- long[time < 25, .(baseline = mean(insulin)), by = .(record_id)]
  long <- merge(long, bl, by = "record_id")
  long[, insulin_norm := insulin / baseline]
  long
}

# Trapezoidal AUC + peak per donor per perifusion phase.
calc_peri_auc <- function(peri_dt) {
  peri_dt <- as.data.table(peri_dt)[order(record_id, time)]
  peri_dt[, phase := fcase(
    time <= 15,                  "baseline",
    time >= 20 & time <= 60,     "stim1",
    time >= 65 & time <= 90,     "wash1",
    time >= 92.5 & time <= 130,  "stim2",
    time >= 135 & time <= 160,   "wash2",
    time >= 162.5 & time <= 180, "stim3",
    time >= 185 & time <= 200,   "wash3",
    default = NA_character_
  )]
  peri_dt[!is.na(phase), .(
    auc  = if (.N >= 2) trapz(time, insulin) else insulin * 5,
    peak = max(insulin, na.rm = TRUE)
  ), by = .(record_id, phase)]
}

# Spread ephys_donor (cell_type x glucose) into wide "_<celltype>_<gluc>" columns.
join_ephys <- function(df, ephys) {
  for (ct in c("Alpha", "Beta", "Delta", "PP")) {
    for (gluc in c(1, 5, 10)) {
      df1 <- ephys[ephys$cell_type == ct & ephys$glucose_mM == gluc, -c(2:3)]
      names(df1)[-1] <- paste0(names(df1)[-1], "_", tolower(ct), "_", gluc)
      df <- dplyr::left_join(df, df1, by = "record_id")
    }
  }
  df
}

# Mitochondrial calc_ outcomes, donor-averaged, from a per-replicate df.
SEAHORSE_REDCAP_COLS <- c(
  "calc_nonmito_oc", "calc_basal_resp", "calc_atp_resp", "calc_proton_leak",
  "calc_max_gluc_resp", "calc_stim_gluc_resp", "calc_max_resp", "calc_spare_cap",
  "calc_gluc_stim_oci", "calc_gluc_stim_sparecap"
)               # the 10 REDCap provides
ALICE_ONLY_COLS <- c("calc_atp_resp_gluc", "calc_bhi", "calc_bhi_gluc")  # UI needs these

seahorse_outcomes <- function(dna_df, cols = SEAHORSE_REDCAP_COLS) {
  d <- dna_df[, c("record_id", cols)]
  aggregate(d[cols], by = list(record_id = d$record_id), FUN = mean, na.rm = TRUE)
}

# Long OCR time course, replicate-averaged (record_id, insulin, time).
seahorse_timecourse <- function(df) {
  tc <- grep("^time", colnames(df), value = TRUE)
  long <- reshape2::melt(df[, c("record_id", tc)], id.vars = "record_id")
  long <- as.data.table(long)
  long[, time := as.numeric(gsub("time_", "", variable))]
  long <- long[, .(insulin = mean(as.numeric(value), na.rm = TRUE)), by = .(record_id, time)]
  as.data.frame(long[, .(record_id, insulin, time)])
}

# Negative/floored-skew log transform (hand-tuned to avoid over-correcting positive,
# high-floor, long-tailed columns like the perifusion KCl AUCs): log10(x + |min| + shift).
# Right-skew columns instead get plain log10(x). Both lists are curated (below), NOT
# auto-detected — auto over-corrects some columns into negative skew.
log_negskew <- function(x, shift) log10(x + abs(min(x, na.rm = TRUE)) + shift)

# Curated skew lists (from the original pipeline). right.skew -> log10(x);
# neg.skew -> log10(x + |min| + 1000). Display CSVs back-transform the right.skew ones (10^x).
RIGHT_SKEW <- c("hba1c","pdinsulinperieq","pdinsulindnaratio","insulindnaratio","auc_baseline_3mmgluc",
                "peak_gluc_15mmgluc","auc_gluc_15mmgluc","ieqperpancreasweight","insulincontent",
                "auc_gluc_6mmgluc","peak_leu_5mmleu","auc_leu_5mmleu","auc_leu_5mmleu_6mmgluc","peak_olp_1p5mmolp",
                "gsis_index_1_10","gsis_index_1_16p7","gsis_index_2p8_16p7","total_insulin_content",
                "insulin_secretion_1","insulin_secretion_10","insulin_secretion_16p7","insulin_secretion_2p8",
                "dnacontent","insulinperieq","pdinsulincontent","avglgcp","avglgpi","avghgcp","avghgpi",
                "avglysatecp","avglysatepi","tg_weight_recovery","calc_max_resp","calc_max_gluc_resp")
NEG_SKEW   <- c("auc_gluc_30mmkcl","auc_leu_30mmkcl","auc_olp_1p5mmolp","auc_olp_1p5mmolp_6mmgluc",
                "auc_olp_30mmkcl","calc_basal_resp","calc_atp_resp","calc_proton_leak",
                "calc_stim_gluc_resp","calc_spare_cap")


## ---- main pipeline ----------------------------------------------------------

redcap_export_fun <- function(api_token,
                              other.tables.path,
                              sqlite.path ="/Users/lzy/humanislet/sqlite/",
                              alice.path   = "/Users/lzy/humanislet/ocr_calculations_2026-04-02.csv",
                              states.path  = "/Users/lzy/humanislet/humanislets/donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/donor_states.csv",
                              demo.path    = "/Users/lzy/humanislet/humanislets/updated_data/ADIIsletCoreHumanIsl-AgeSexBMIHbA1CDiabet_DATA_2026-08-01_2009.csv") {

  library(RSQLite); library(dplyr); library(data.table)
  library(pzfx); library(REDCapR); library(readxl); library(caTools)

  api_url <- "https://redcap.ualberta.ca/api/"
  cur.dir <- getwd()

  demographics <- read.csv(demo.path)   # ADIIsletCoreHumanIsl age/sex/BMI/HbA1c/diabetes

  ## ---- 1. REDCap downloads --------------------------------------------------
  dl <- function(field, file) redcap_file_download_oneshot(
    directory = cur.dir, redcap_uri = api_url, token = api_token,
    record = "R000", field = field, overwrite = TRUE, file_name = file)
  dl("seahorse_data_normalized_dna",      "seahorse.csv")
  dl("perifusion_data",                   "perifusion.pzfx")
  dl("seahorse_data_normalized_baseline", "seahorse_norm.csv")

  report.ids   <- c(32427, 32428, 32429, 32431, 32438, 37603, 40523)
  report.names <- c("donor", "distribution", "exocytosis", "isolation", "gsis", "ephys", "inventory")
  tables <- list()
  for (i in seq_along(report.ids)) {
    dat <- redcap_report(redcap_uri = api_url, token = api_token,
                         report_id = report.ids[i], guess_type = FALSE)
    while (sink.number() > 0) sink()
    dat <- as.data.frame(dat$data[-1, ])
    dat[] <- lapply(dat, function(x) { x[x %in% c("", "no data", "INF", "Inf", "no image available.jpg")] <- NA; x })
    tables[[report.names[i]]] <- dat
  }

  ## ---- 2a. isolation --------------------------------------------------------
  iso <- tables$isolation
  char.cols <- c("record_id", "isolationdate", "collagenase", "collagenasetype",
                 "dithizonestaining", "dithizonestaining2", "dithizonestaining3",
                 "dithizonestaining4", "fatinfiltration", "pancreasconsistency")
  num.cols <- setdiff(colnames(iso), char.cols)
  iso[, num.cols] <- apply(iso[, num.cols], 2, as.numeric)
  # NOTE: index with as.character() — some of these columns were coerced to numeric
  # above, and numeric indexing is positional (v[0] drops rows). as.character() forces
  # name-matching so lengths are preserved and unknown/NA values map to NA.
  iso$collagenase     <- unname(c("1"="Roche","2"="VitaCyte","3"="Serva","4"="Sigma")[as.character(iso$collagenase)])
  iso$collagenasetype <- unname(c("1"="Liberase","2"="Clzyme","3"="Gold","4"="NB1","5"="Sigma","6"="Other","7"="Recombinant")[as.character(iso$collagenasetype)])
  iso$purifiedtissue      <- unname(c("0"="No","1"="Yes")[as.character(iso$purifiedtissue)])
  iso$fatinfiltration     <- unname(c("0"="No","1"="Yes")[as.character(iso$fatinfiltration)])
  iso$pancreasconsistency <- unname(c("1"="Soft","2"="Normal","3"="Fibrotic","4"="Inconsistent")[as.character(iso$pancreasconsistency)])
  inventory <- tables$inventory
  inventory$embeddedbiopsy <- unname(c("0"="No","1"="Yes")[as.character(inventory$embeddedbiopsy)])
  iso <- merge(iso, inventory, by = "record_id", all = TRUE)
  iso <- iso[grepl("^R", iso$record_id), ]
  tables$inventory <- NULL
  tables$isolation <- iso

  ## ---- 2b. donor ------------------------------------------------------------
  donor <- tables$donor
  colnames(donor)[colnames(donor) == "medicalconditions___1"] <- "T1_diabetes"
  colnames(donor)[colnames(donor) == "medicalconditions___2"] <- "T2_diabetes"
  donor <- donor[, -grep("medicalconditions", colnames(donor))]
  donor <- donor[grepl("^R", donor$record_id), ]
  donor$hla_a2 <- ifelse(is.na(donor$hlaa), NA, 0)
  hlaa <- strsplit(donor$hlaa, ", "); names(hlaa) <- donor$record_id
  hlaa <- hlaa[lengths(hlaa) > 1] %>% as.data.frame() %>% t() %>% as.data.frame()
  hlaa.donors <- rownames(hlaa)[hlaa$V1 == "2" | hlaa$V2 == "2"]
  donor$hla_a2[donor$record_id %in% hlaa.donors] <- 1
  donor$diagnosis <- NA
  donor$diagnosis[donor$T1_diabetes == "1"] <- "Type1"
  donor$diagnosis[donor$T2_diabetes == "1"] <- "Type2"
  donor$diagnosis[donor$T1_diabetes == "0" & donor$T2_diabetes == "0"] <- "None"
  donor$diagnosis_computed <- donor$diagnosis
  h <- as.numeric(donor$hba1c)
  donor$diagnosis_computed[donor$diagnosis == "None" & h > 6.5] <- "Type2"
  donor$diagnosis_computed[donor$diagnosis == "None" & h < 6.5 & h > 5.7] <- "Pre.T2D"
  donor$yearsdiabetic <- demographics$yearsdiabetic[match(donor$record_id, demographics$record_id)]
  donor$yearsdiabetic[donor$yearsdiabetic == "unknown"] <- NA
  donor <- merge(donor, tables$distribution, by = "record_id")
  donor$donorsex     <- unname(c("1"="Male","2"="Female","3"="Unknown")[as.character(donor$donorsex)])
  donor$donationtype <- unname(c("1"="NDD","2"="DCD","3"="MAID")[as.character(donor$donationtype)])
  donor$hla_a2       <- unname(c("0"="Negative","1"="Positive")[as.character(donor$hla_a2)])
  tables$distribution <- NULL
  char.cols <- c("record_id","rrid","donorsex","donationtype","hlaa","hlab","hlabw","hlac",
                 "hlacw","hladrb1","hladr","hladqb1","hladqa1","hladpb1","dpa1","hlaother",
                 "hla_a2","T1_diabetes","T2_diabetes","diagnosis","diagnosis_computed","yearsdiabetic")
  num.cols <- setdiff(colnames(donor), char.cols)
  donor[, char.cols] <- lapply(donor[, char.cols], as.character)
  donor[, num.cols]  <- lapply(donor[, num.cols], as.numeric)
  donor[donor == ""] <- NA
  tables$donor <- donor

  ## ---- 2c. exocytosis -------------------------------------------------------
  exo <- tables$exocytosis
  exo <- exo[, grep("exocytosis", colnames(exo), invert = TRUE)]
  exo <- reshape2::melt(exo, id.vars = "record_id", measure.vars = colnames(exo)[-1])
  exo <- na.omit(exo); colnames(exo) <- c("record_id", "exposure", "insulin_exo")
  exo <- exo[order(exo$record_id), ]
  exo$cell_id <- ave(exo$record_id, exo$record_id, FUN = function(x) paste0(x, "_", seq_along(x)))
  exo <- exo[, c("record_id", "cell_id", "exposure", "insulin_exo")]
  exo$exposure <- c(lgexo = "gluc_1", mgexo = "gluc_5", hgexo = "gluc_10")[
    sub(".*(lgexo|mgexo|hgexo).*", "\\1", exo$exposure)]
  exo$insulin_exo <- pmax(as.numeric(exo$insulin_exo), 0)
  tables$exocytosis <- exo

  ## ---- 2d. gsis -------------------------------------------------------------
  gsis <- tables$gsis
  gsis.ct <- gsis[, c("record_id", "culturetime2")]
  gsis <- data.table::melt(
    setDT(gsis), id.vars = "record_id",
    measure.vars = list(
      c("lowmedcontent12","lowmedcontent22","lowmedcontent32"),
      c("lowglucose12_60a3b7","lowglucose22","lowglucose32"),
      c("medglucose12","medglucose22","medglucose32"),
      c("highcontent12","highcontent22","highcontent32"),
      c("lowglucose72","lowglucose82","lowglucose92"),
      c("highglucose12","highglucose22","highglucose32"),
      c("content28167_1","content28167_2","content28167_3"),
      c("value28_1","value28_2","value28_3"),
      c("value167_1","value167_2","value167_3")),
    value.name = c("insulin.content.1","first.1","second.1","insulin.content.2","first.2",
                   "second.2","insulin.content.3","first.3","second.3"),
    variable.name = "replicate")
  mk <- function(cols, g1, g2) {
    d <- gsis[, c(1, 2, cols), with = FALSE]
    d <- d[apply(d, 1, function(x) sum(is.na(x)) != 3), ]
    setnames(d, 3:5, c("total_insulin_content", "first_insulin_secretion", "second_insulin_secretion"))
    d$first_gluc_conc <- g1; d$second_gluc_conc <- g2; d
  }
  gsis <- rbindlist(list(mk(3:5, 1, 10), mk(6:8, 1, 16.7), mk(9:11, 2.8, 16.7)))
  gsis <- gsis[order(record_id)]
  gsis[, replicate := as.character(rowid(record_id))]
  gsis <- as.data.frame(gsis)
  num.cols <- setdiff(colnames(gsis), c("record_id", "replicate"))
  gsis[, num.cols] <- lapply(gsis[, num.cols], as.numeric)
  gsis$first_insulin_percent  <- gsis$first_insulin_secretion  / gsis$total_insulin_content * 100
  gsis$second_insulin_percent <- gsis$second_insulin_secretion / gsis$total_insulin_content * 100
  gsis$gluc_conc_group <- paste0(gsis$first_gluc_conc, "_", gsis$second_gluc_conc)
  gsis$stim_index <- gsis$second_insulin_secretion / gsis$first_insulin_secretion
  gsis <- merge(gsis, gsis.ct, by = "record_id")
  tables$gsis <- gsis

  ## ---- 2e. ephys ------------------------------------------------------------
  ephys <- tables$ephys
  ephys <- ephys[ephys$epnumber > 0, ]
  ephys <- reshape2::melt(ephys, id.vars = "record_id", measure.vars = colnames(ephys)[-1])
  ephys <- na.omit(ephys); ephys$variable <- as.character(ephys$variable)
  ephys <- ephys[!ephys$variable %in% c("epnumber", "electrophysiology_complete"), ]
  ephys$cell_id  <- gsub("[^0-9.-]", "", ephys$variable)
  ephys$variable <- gsub("[0-9.+]", "", ephys$variable)
  ephys <- reshape2::dcast(ephys, record_id + cell_id ~ variable, value.var = "value")
  ephys <- ephys[ephys$eptype %in% c("1", "2"), ]
  ephys$eptype <- c("1"="Beta","2"="Alpha")[ephys$eptype]
  ephys <- ephys[, c("record_id","cell_id","eptype","epsize","epntc","epnfdc","epnldc",
                     "nci","epnpsca","ephisc","epnepcca","epnlcca")]
  ephys[, 4:12] <- lapply(ephys[, 4:12], as.numeric)
  colnames(ephys) <- c("record_id","cell_id","eptype","cell_size_pF","total_exocytosis_fF_pF",
                       "early_exocytosis_fF_pF","late_exocytosis_fF_pF","calcium_entry_pC_pF",
                       "na_current_amp_pA_pF","na_half_inactivation_mV","early_ca_current_pA_pF",
                       "late_ca_current_pA_pF")
  for (c in c("total_exocytosis_fF_pF","early_exocytosis_fF_pF","late_exocytosis_fF_pF"))
    ephys[[c]] <- pmax(ephys[[c]], 0)
  for (c in c("calcium_entry_pC_pF","na_current_amp_pA_pF","early_ca_current_pA_pF","late_ca_current_pA_pF"))
    ephys[[c]] <- -ephys[[c]]
  tables$ephys <- ephys

  ## ---- 2f. seahorse (REDCap) ------------------------------------------------
  dna      <- read_seahorse_wide("seahorse.csv")        # DNA-normalized
  baseline <- read_seahorse_wide("seahorse_norm.csv")   # baseline-OC-normalized
  hybrid <- baseline %>%                                 # baseline + DNA calc_basal_resp
    select(-any_of("calc_basal_resp")) %>%
    left_join(dna %>% select(record_id, replicate, calc_basal_resp),
              by = c("record_id", "replicate")) %>%
    arrange(record_id)
  tables$seahorse_norm_dna            <- dna
  tables$seahorse_norm_dna_baselineoc <- baseline
  tables$seahorse                     <- hybrid

  ## ---- 2g. perifusion traces (raw wide, for sqlite) -------------------------
  pf <- function(t) read_pzfx("perifusion.pzfx", table = t)
  tables$peri_gluc <- process_peri_trace(pf(1), pf(4))
  tables$peri_leu  <- process_peri_trace(pf(2), pf(5))
  tables$peri_olp  <- process_peri_trace(pf(3), pf(6), drop_na_time = TRUE)

  ## ---- 2h. externally-provided tables ---------------------------------------
  read_num <- function(f) {
    d <- read.csv(paste0(other.tables.path, f))
    d[, -1] <- apply(d[, -1], 2, as.numeric)   # keep ALL donor rows, incl. all-NA
    d
  }
  prohormone <- read.csv(paste0(other.tables.path, "updated_data/ADIIsletCoreHumanIsl-ProhormoneExpression_DATA_2026-08-01_2224.csv"))
  prohormone$redcap_data_access_group <- NULL
  prohormone$prohormone_expression_complete <- NULL
  prohormone[] <- lapply(prohormone, function(x) { x[x %in% c("", "N/a")] <- NA; x })
  prohormone[, -1] <- apply(prohormone[, -1], 2, as.numeric)   # keep ALL donor rows, incl. all-NA
  prohormone <- prohormone[grepl("^R", prohormone$record_id), , drop = FALSE]  # cohort donors only: drop M### mouse rows (varValues ^R filter doesn't reach this table)
  tables$prohormone  <- prohormone
  tables$grs         <- read_num("display_data/numerical_grs.csv")
  tables$ancestry    <- read_num("display_data/numerical_ancestry.csv")
  lip_extract <- read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))
  lip_extract[, -1] <- apply(lip_extract[, -1], 2, as.numeric)
  tables$lip_extract <- lip_extract

  raw_variable_summary  <- read.csv(paste0(other.tables.path, "outcomes_processing_input/raw_variable_summary.csv"))
  proc_variable_summary <- read.csv(paste0(other.tables.path, "display_interface/proc_variable_summary_v2.csv"))
  patchseq              <- readRDS(paste0(other.tables.path, "outcomes_processing_input/patchseq_metadata_v2.rds"))
  computed              <- read.csv(paste0(other.tables.path, "outcomes_processing_input/composition_proteomics_v2.csv"))

  ## ---- 2i. combine ephys (redcap + patchseq + exocytosis) -------------------
  ephys.rc <- tables$ephys
  ephys.rc$cell_id <- paste0(ephys.rc$record_id, "_", ephys.rc$cell_id, "_rc")
  ephys.rc$glucose_mM <- "5"; colnames(ephys.rc)[3] <- "cell_type"
  ephys.exo <- tables$exocytosis
  ephys.exo$cell_id <- paste0(ephys.exo$cell_id, "_exo")
  ephys.exo$exposure <- gsub("gluc_", "", ephys.exo$exposure)
  ephys.exo$cell_type <- "Beta"
  colnames(ephys.exo)[4] <- "total_exocytosis_fF_pF"; colnames(ephys.exo)[3] <- "glucose_mM"
  ephys.all <- merge(ephys.rc, patchseq, by = colnames(ephys.rc), all = TRUE)
  ephys.all <- merge(ephys.all, ephys.exo, by = colnames(ephys.exo), all = TRUE)
  ephys.all <- ephys.all[grep("R", ephys.all$record_id), ]
  ephys.all <- ephys.all[, c(1, 2, 5, 3, 6, 4, 7:13)]
  ephys.all[, 5:13] <- apply(ephys.all[, 5:13], 2, as.numeric)
  ephys.dt <- as.data.table(ephys.all[, -2])[, lapply(.SD, mean), by = .(record_id, cell_type, glucose_mM)]
  ephys.dt <- as.data.frame(ephys.dt)
  colnames(ephys.dt)[4:12] <- paste0(colnames(ephys.dt)[4:12], "_donor")
  colnames(patchseq)[5:13] <- paste0(colnames(patchseq)[5:13], "_cell")
  tables$ephys_donor <- ephys.dt
  tables$ephys_cell  <- patchseq
  tables$computed    <- computed
  tables$proc_variable_summary <- proc_variable_summary
  tables$exocytosis <- NULL; tables$ephys <- NULL

  ## ---- 3. mitochondrial outcomes: REDCap (10) + Alice (3) -------------------
  alice <- read.csv(alice.path)
  alice_wide <- reshape2::dcast(
    alice[alice$normalisation_type == "Norm DNA", c("donor_id", "measurement", "value")],
    donor_id ~ measurement, value.var = "value")
  colnames(alice_wide)[1] <- "record_id"
  # Mito outcomes: REDCap DNA-normalized is PRIMARY; Alice "Norm DNA" only FILLS GAPS.
  # As of 2026-08, 22 donors (R525,R526,R536,R541,R542,R543,R546,R547,R551,R552,R554,
  # R555,R556,R583,R585,R587,R591,R599,R608,R612,R614,R623) have NO DNA-normalized
  # seahorse in the REDCap export (seahorse.csv is empty for them; only baseline exists),
  # but Alice's ocr_calculations "Norm DNA" has them, on the SAME scale as REDCap
  # (verified: Alice/REDCap ratio = 1.000 across all shared donors). We therefore fill
  # the 10 REDCap mito columns from Alice ONLY where REDCap is NA; a real REDCap value is
  # never overridden. (Alice also supplies the 3 UI-only cols REDCap doesn't compute.)
  #
  # TODO / self-deactivating: this gap-fill touches NAs ONLY. Once these donors'
  # DNA-normalized seahorse is available in REDCap, their Alice values stop being used
  # automatically — no code change needed. i.e. when the data can be pulled from REDCap
  # we no longer add ANY Alice mito data; this block can then be deleted for clarity.
  redcap   <- seahorse_outcomes(dna)                              # 10 REDCap, donor-mean
  all_ids  <- union(redcap$record_id, alice_wide$record_id)       # keep Alice-only donors (e.g. R599)
  ri       <- match(all_ids, redcap$record_id)
  ali      <- match(all_ids, alice_wide$record_id)
  mito_out <- data.frame(record_id = all_ids, stringsAsFactors = FALSE)
  for (c in SEAHORSE_REDCAP_COLS) {                               # REDCap primary, Alice fallback
    rv <- if (c %in% colnames(redcap))     redcap[[c]][ri]                                  else NA
    av <- if (c %in% colnames(alice_wide)) suppressWarnings(as.numeric(alice_wide[[c]][ali])) else NA
    mito_out[[c]] <- ifelse(is.na(rv), av, rv)
  }
  for (c in ALICE_ONLY_COLS)                                      # 3 Alice-only (not in REDCap)
    mito_out[[c]] <- if (c %in% colnames(alice_wide)) suppressWarnings(as.numeric(alice_wide[[c]][ali])) else NA
  MITO_COLS <- c(SEAHORSE_REDCAP_COLS, ALICE_ONLY_COLS)

  ## ---- 3b. perifusion AUC summary (function_summary) — needed by proc_metadata ----
  gluc.norm <- normalize_peri_trace(tables$peri_gluc)
  leu.norm  <- normalize_peri_trace(tables$peri_leu)
  olp.norm  <- normalize_peri_trace(tables$peri_olp)
  auc_cols <- function(nrm) {
    a <- calc_peri_auc(nrm[, .(record_id, insulin, time)])
    list(wide = dcast(a, record_id ~ phase, value.var = "auc"),
         peak = dcast(a, record_id ~ phase, value.var = "peak"))
  }
  g <- auc_cols(gluc.norm); l <- auc_cols(leu.norm); o <- auc_cols(olp.norm)
  function_summary <- Reduce(function(a, b) merge(a, b, by = "record_id", all = TRUE), list(
    data.frame(record_id = g$wide$record_id, auc_baseline_3mmgluc = g$wide$baseline,
               peak_gluc_15mmgluc = g$peak$stim1, auc_gluc_15mmgluc = g$wide$stim1,
               auc_gluc_6mmgluc = g$wide$stim2, auc_gluc_30mmkcl = g$wide$stim3),
    data.frame(record_id = l$wide$record_id, peak_leu_5mmleu = l$peak$stim1,
               auc_leu_5mmleu = l$wide$stim1, auc_leu_5mmleu_6mmgluc = l$wide$stim2,
               auc_leu_30mmkcl = l$wide$stim3),
    data.frame(record_id = o$wide$record_id, peak_olp_1p5mmolp = o$peak$stim1,
               auc_olp_1p5mmolp = o$wide$stim1, auc_olp_1p5mmolp_6mmgluc = o$wide$stim2,
               auc_olp_30mmkcl = o$wide$stim3)))
  tables$function_summary <- function_summary

  ## ---- 4. build proc_metadata (varValues) -----------------------------------
  varValues <- NULL
  for (tb in unique(raw_variable_summary$table)) {
    # ephys is single-cell / multi-outcome-per-donor — it's added later via join_ephys,
    # NOT here. Skip any table not present in `tables` (e.g. tables$ephys is NULL).
    if (is.null(tables[[tb]])) next
    cols <- raw_variable_summary$column[raw_variable_summary$table == tb]
    # only pull columns that actually exist (e.g. numerical_grs.csv dropped qu22/
    # sharp21/udler18, but raw_variable_summary may still list them) — avoids a crash.
    have <- intersect(cols, colnames(tables[[tb]]))
    if (length(have) < length(cols))
      message(sprintf("proc_metadata: %s missing %d/%d cols from raw_variable_summary: %s",
                      tb, length(cols) - length(have), length(cols),
                      paste(setdiff(cols, have), collapse = ", ")))
    dat  <- tables[[tb]][, c("record_id", have)]
    if (tb == "seahorse") {
      dat <- aggregate(dat[setdiff(colnames(dat), "record_id")],
                       by = dat["record_id"], FUN = mean)
    } else if (tb == "gsis") {
      d <- tables$gsis
      secr <- rbind(
        setNames(d[, c("record_id","first_gluc_conc","first_insulin_secretion")],  c("record_id","gluc_conc","insulin_secretion")),
        setNames(d[, c("record_id","second_gluc_conc","second_insulin_secretion")], c("record_id","gluc_conc","insulin_secretion")))
      secr$gluc_conc <- paste0("insulin_secretion_", gsub("\\.", "p", secr$gluc_conc))
      med.secr <- reshape2::dcast(aggregate(secr["insulin_secretion"], secr[c("record_id","gluc_conc")], median, na.rm = TRUE),
                                  record_id ~ gluc_conc, value.var = "insulin_secretion")
      st <- aggregate(d["stim_index"], d[c("record_id","gluc_conc_group")], median, na.rm = TRUE)
      st$gluc_conc_group <- paste0("gsis_index_", gsub("\\.", "p", st$gluc_conc_group))
      med.stim <- reshape2::dcast(st, record_id ~ gluc_conc_group, value.var = "stim_index")
      med.ins  <- aggregate(d["total_insulin_content"], d["record_id"], median, na.rm = TRUE)
      ct       <- distinct(d[, c("record_id", "culturetime2")])
      dat <- Reduce(function(a, b) merge(a, b, by = "record_id", all = TRUE),
                    list(med.secr, med.stim, med.ins, ct))
      dat <- dat %>% mutate_all(~ ifelse(is.nan(.) | is.infinite(.), NA, .))
      dat$culturetime2 <- as.numeric(dat$culturetime2)
      # drop donors > 4 SD from mean (on log scale), per secretion/index column
      gout <- dat[, setdiff(colnames(dat), c("record_id", "culturetime2"))]
      for (cc in colnames(gout)) {
        lx <- log10(gout[[cc]]); m <- mean(lx, na.rm = TRUE); s <- sd(lx, na.rm = TRUE)
        dat[which(lx > m + 4*s | lx < m - 4*s), cc] <- NA
      }
    }
    # normalize the merge key: stray whitespace (e.g. "R307 ") would otherwise
    # create a phantom, unmerged duplicate row for that donor.
    dat$record_id <- trimws(as.character(dat$record_id))
    varValues <- if (is.null(varValues)) dat else merge(varValues, dat, by = "record_id", all = TRUE)
  }
  # Cohort donors only (R000-R650). Non-R ids (e.g. M### rows) leak in from source
  # tables through the all=TRUE merges above; drop them here so proc_metadata and every
  # varValues-derived output (numerical_donor_info, metadata norm/raw, etc.) stay R-only.
  varValues <- varValues[grepl("^R", varValues$record_id), , drop = FALSE]
  varValues_raw <- varValues   # pre-log copy

  # Curated skew normalization: right.skew -> log10(x); neg.skew -> log10(x + |min| + 1000).
  for (c in intersect(RIGHT_SKEW, colnames(varValues))) varValues[[c]] <- log10(varValues[[c]])
  for (c in intersect(NEG_SKEW,   colnames(varValues))) varValues[[c]] <- log_negskew(varValues[[c]], 1000)
  varValues <- varValues %>% mutate_all(~ ifelse(is.nan(.) | is.infinite(.), NA, .))
  tables$proc_metadata <- varValues

  ## ---- 5. numerical_* display CSVs ------------------------------------------
  wcsv <- function(x, f) write.csv(x, paste0(other.tables.path, "display_data/", f), row.names = FALSE)
  # back-transform the right.skew (plain log10) columns for raw-unit display CSVs.
  # neg.skew columns are never in the display var lists, so they need no reversal.
  unskew <- function(df, cols = names(df)) {
    for (c in intersect(cols, RIGHT_SKEW)) df[[c]] <- 10^df[[c]]
    df
  }

  donor.info <- unskew(varValues[, c("record_id","donorage","bodymassindex","hba1c","yearsdiabetic")])
  wcsv(donor.info, "numerical_donor_info.csv")

  iso.vars <- c("record_id","coldischemiatime","puritypercentage","trappedpercentage","pancreasweight",
                "digesttime","totalieq","isletparticleindex","ieqperpancreasweight","insulincontent",
                "insulinperieq","predistributionculturetime","percentieqrecoverypostculture","pdisletparticleindex",
                "pdinsulinperieq","cryotubesremaining","sftubesremaining")
  iso.info <- unskew(varValues[, iso.vars], iso.vars)
  iso.info$insulincontent <- iso.info$insulincontent / 1000
  wcsv(iso.info, "numerical_isolation_info.csv")

  gsis.vars <- c("record_id","culturetime2","total_insulin_content","insulin_secretion_1","insulin_secretion_2p8",
                 "insulin_secretion_10","insulin_secretion_16p7","gsis_index_1_10","gsis_index_1_16p7","gsis_index_2p8_16p7")
  gsis.info <- varValues[, gsis.vars]
  gsis.info$culturetime2 <- as.numeric(gsis.info$culturetime2)
  gsis.info <- unskew(gsis.info, gsis.vars)
  for (c in c("total_insulin_content","insulin_secretion_1","insulin_secretion_10","insulin_secretion_16p7","insulin_secretion_2p8"))
    gsis.info[[c]] <- gsis.info[[c]] / 1e6      # pg/mL -> ug/mL
  wcsv(gsis.info, "numerical_gsis.csv")

  wcsv(seahorse_timecourse(dna),    "mito_func.csv")       # DNA trace
  wcsv(seahorse_timecourse(hybrid), "mito_func_norm.csv")  # baseline/hybrid trace
  wcsv(mito_out[, c("record_id", MITO_COLS)], "numerical_mito.csv")  # REDCap 10 + Alice 3

  # prohormone: raw avg/ratio columns — the donor-view prohormone histogram loads this
  proh_cols <- grep("avg|ratio", colnames(prohormone), value = TRUE)
  wcsv(prohormone[, c("record_id", proh_cols)], "numerical_prohormone.csv")

  # lipid extraction: raw values — numerical_ version for the donor-view lipid histogram
  wcsv(lip_extract, "numerical_lip_extract.csv")

  ## ---- 6. write perifusion traces + function_summary (computed in 3b) --------
  wcsv(gluc.norm, "peri_gluc.csv"); wcsv(leu.norm, "peri_leu.csv"); wcsv(olp.norm, "peri_olp.csv")
  write.csv(function_summary, paste0(other.tables.path, "outcomes_processing_input/function_summary.csv"), row.names = FALSE)

  ## ---- 7. ephys numerical -----------------------------------------------------
  wcsv(ephys.dt, "numerical_ephys.csv")
  ephys.all.out <- ephys.all
  colnames(ephys.all.out)[5:13] <- paste0(colnames(ephys.all.out)[5:13], "_donor")
  wcsv(ephys.all.out, "numerical_ephys_cell.csv")

  ## ---- 8. metadata_sum (built in memory, written once) ----------------------
  # NORM mito calc_: curated right.skew -> log10(x); neg.skew -> log10(x + |min| + 1).
  mito_norm <- mito_out
  for (c in intersect(c("calc_max_resp","calc_max_gluc_resp","calc_nonmito_oc"), colnames(mito_norm)))
    mito_norm[[c]] <- log10(mito_norm[[c]])
  for (c in intersect(c("calc_basal_resp","calc_atp_resp","calc_proton_leak","calc_stim_gluc_resp",
                        "calc_spare_cap","calc_atp_resp_gluc","calc_bhi","calc_bhi_gluc"), colnames(mito_norm)))
    mito_norm[[c]] <- log_negskew(mito_norm[[c]], 1)
  mito_norm <- mito_norm %>% mutate_all(~ ifelse(is.nan(.) | is.infinite(.), NA, .))
  # NORM function_summary: log10 all AUC/peak columns (they are all right-skewed).
  fs_norm  <- function_summary
  fs_cols  <- setdiff(colnames(function_summary), "record_id")
  fs_norm[fs_cols] <- lapply(fs_norm[fs_cols], log10)
  fs_norm  <- fs_norm %>% mutate_all(~ ifelse(is.nan(.) | is.infinite(.), NA, .))

  set_cols <- function(df, src, cols) {
    for (c in cols) df[[c]] <- src[[c]][match(df$record_id, src$record_id)]
    df
  }
  metadata_norm <- join_ephys(varValues,     ephys.dt)
  metadata_norm <- set_cols(metadata_norm, mito_norm, intersect(MITO_COLS, colnames(metadata_norm)))
  metadata_norm <- set_cols(metadata_norm, fs_norm,   intersect(fs_cols,   colnames(metadata_norm)))
  # ensure the 3 Alice-only calc_ columns are present even if not in raw_variable_summary
  for (c in setdiff(MITO_COLS, colnames(metadata_norm)))
    metadata_norm[[c]] <- mito_norm[[c]][match(metadata_norm$record_id, mito_norm$record_id)]
  wcsv(metadata_norm, "metadata_sum_norm.csv")

  metadata_raw <- join_ephys(varValues_raw, ephys.dt)
  metadata_raw <- set_cols(metadata_raw, mito_out,        intersect(MITO_COLS, colnames(metadata_raw)))
  metadata_raw <- set_cols(metadata_raw, function_summary, intersect(fs_cols,  colnames(metadata_raw)))
  for (c in setdiff(MITO_COLS, colnames(metadata_raw)))
    metadata_raw[[c]] <- mito_out[[c]][match(metadata_raw$record_id, mito_out$record_id)]
  wcsv(metadata_raw, "metadata_sum_raw.csv")

  ## ---- 9. Spearman correlation ----------------------------------------------
  library(Hmisc)
  df_clean <- metadata_norm[, -1]
  df_clean[] <- lapply(df_clean, function(x) if (is.character(x) || is.factor(x)) as.numeric(factor(x)) else x)
  df_clean <- df_clean[, colSums(!is.na(df_clean)) >= 6]
  rc <- rcorr(as.matrix(df_clean), type = "spearman")
  cor_df <- reshape2::melt(rc$r)
  cor_df$pval   <- reshape2::melt(rc$P)$value
  cor_df$npairs <- reshape2::melt(rc$n)$value
  cor_df$Var1 <- as.character(cor_df$Var1); cor_df$Var2 <- as.character(cor_df$Var2)
  metadata_corr <- cor_df %>% filter(Var1 != Var2) %>% rowwise() %>%
    mutate(pair_id = paste(sort(c(Var1, Var2)), collapse = "_")) %>% ungroup() %>%
    distinct(pair_id, .keep_all = TRUE) %>% select(Var1, Var2, value, pval, npairs)
  metadata_corr <- metadata_corr[!is.na(metadata_corr$pval) & metadata_corr$npairs > 5, ]
  wcsv(metadata_corr, "metadata_corr.csv")

  ## ---- 10. write sqlite -----------------------------------------------------
  # seahorse_norm_dna: numeric-coerce; seahorse_outcome = mito_out (REDCap 10 + Alice 3)
  tables$seahorse_norm_dna[, -c(1:2)] <- lapply(tables$seahorse_norm_dna[, -c(1:2)], as.numeric)
  tables$seahorse_outcome <- mito_out
  tables$raw_variable_summary <- raw_variable_summary
  # proc_metadata = normalized metadata + donor-state label (final_cluster) from the
  # endotype deliverable's `cluster` column, matched by donor_id. This is the ONLY
  # place final_cluster is attached (the metadata_sum_*.csv files never carry it).
  donor_states <- read.csv(states.path)
  proc_metadata <- metadata_norm
  # the interface expects final_cluster as "C0".."C4" (atlas/kg-search filter on that);
  # donor_states$cluster is 0..4, so prefix "C". NA stays NA (donors with no cluster).
  cl <- donor_states$cluster[match(proc_metadata$record_id, donor_states$donor_id)]
  proc_metadata$final_cluster <- ifelse(is.na(cl), NA, paste0("C", cl))
  tables$proc_metadata <- proc_metadata

  # also carry final_cluster on the donor table (the old pipeline had it there, and the
  # interface may read the label from `donor`). Same 0..4 -> "C0".."C4" mapping.
  cld <- donor_states$cluster[match(tables$donor$record_id, donor_states$donor_id)]
  tables$donor$final_cluster <- ifelse(is.na(cld), NA, paste0("C", cld))

  sqlite.tables <- c("computed","donor","ephys_cell","ephys_donor","gsis","isolation",
                     "proc_metadata","proc_variable_summary","raw_variable_summary",
                     "seahorse","seahorse_norm_dna","seahorse_norm_dna_baselineoc","seahorse_outcome",
                     "function_summary","peri_gluc","peri_leu","peri_olp",
                     "grs","ancestry","prohormone","lip_extract")
  con <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  for (nm in sqlite.tables) dbWriteTable(con, nm, tables[[nm]], overwrite = TRUE)
  dbDisconnect(con)

  invisible(tables)
}
