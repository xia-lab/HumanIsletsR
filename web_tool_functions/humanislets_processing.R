# R Functions for HumanIslets web tool
# Author: Jessica Ewald

################################################################################

redcap_export_fun <- function(api_token){
  
  library(RSQLite)
  library(dplyr)
  library(data.table)
  library(pzfx)
  library(REDCapR)
  library(readxl)
  
  api_url <- "https://redcap.ualberta.ca/api/";
  cur.dir <- getwd()
  
  report.ids <- c(32427, 32428, 32429, 32431, 32438, 37603, 40523)
  report.names <- c("donor", "distribution", "exocytosis", "isolation", "gsis", "ephys", "inventory")
  
  # download whole files
  redcap_download_file_oneshot(directory = cur.dir,
                               redcap_uri = api_url,
                               token = api_token,
                               record = "R000",
                               field = "seahorse_data_normalized_dna",
                               overwrite = TRUE,
                               file_name = "seahorse.csv")
  
  redcap_download_file_oneshot(directory = cur.dir,
                               redcap_uri = api_url,
                               token = api_token,
                               record = "R000",
                               field = "perifusion_data",
                               overwrite = TRUE,
                               file_name = "perifusion.pzfx")
  
  redcap_download_file_oneshot(directory = cur.dir,
                               redcap_uri = api_url,
                               token = api_token,
                               record = "R000",
                               field = "seahorse_data_normalized_baseline",
                               overwrite = TRUE,
                               file_name = "seahorse_norm.csv")
  
  # get tables from report
  tables <- list()
  for(i in c(1:length(report.ids))){
    
    dat <- redcap_report(redcap_uri = api_url,
                         token = api_token, 
                         report_id = report.ids[i],
                         guess_type = FALSE)
    dat <- dat$data
    dat <- dat[-1, ] # remove donor R000
    
    # handle missing values
    dat[is.na(dat)] <- NA
    dat[dat == ""] <- NA
    dat[dat == "no data"] <- NA
    dat[dat == "INF"] <- NA
    dat[dat == "Inf"] <- NA
    dat[dat == "no image available.jpg"] <- NA
    tables[[i]] <- dat
  }
  names(tables) <- report.names
  
  
  #### TABLE-SPECIFIC FORMATTING ####
  
  ### isolation
  iso <- tables$isolation
  char.cols <- c("record_id", "isolationdate", "collagenase", "collagenasetype",
                 "dithizonestaining", "dithizonestaining2", "dithizonestaining3",
                 "dithizonestaining4", "fatinfiltration", "pancreasconsistency")
  num.cols <- colnames(iso)[!(colnames(iso) %in% char.cols)]
  
  iso[,num.cols] <- apply(iso[,num.cols], 2, as.numeric)
  
  # give factors meaningful names
  iso$collagenase[iso$collagenase == "1"] <- "Roche"
  iso$collagenase[iso$collagenase == "2"] <- "VitaCyte"
  iso$collagenase[iso$collagenase == "3"] <- "Serva"
  iso$collagenase[iso$collagenase == "4"] <- "Sigma"
  
  iso$collagenasetype[iso$collagenasetype == "1"] <- "Liberase"
  iso$collagenasetype[iso$collagenasetype == "2"] <- "Clzyme"
  iso$collagenasetype[iso$collagenasetype == "3"] <- "Gold"
  iso$collagenasetype[iso$collagenasetype == "4"] <- "NB1"
  iso$collagenasetype[iso$collagenasetype == "5"] <- "Sigma"
  iso$collagenasetype[iso$collagenasetype == "6"] <- "Other"
  iso$collagenasetype[iso$collagenasetype == "7"] <- "Recombinant"
  
  iso$purifiedtissue[iso$purifiedtissue == "0"] <- "No"
  iso$purifiedtissue[iso$purifiedtissue == "1"] <- "Yes"
  
  iso$fatinfiltration[iso$fatinfiltration == "0"] <- "No"
  iso$fatinfiltration[iso$fatinfiltration == "1"] <- "Yes"
  
  iso$pancreasconsistency[iso$pancreasconsistency == "1"] <- "Soft"
  iso$pancreasconsistency[iso$pancreasconsistency == "2"] <- "Normal"
  iso$pancreasconsistency[iso$pancreasconsistency == "3"] <- "Fibrotic"
  iso$pancreasconsistency[iso$pancreasconsistency == "4"] <- "Inconsistent"
  
  # merge with inventory
  inventory <- tables$inventory
  
  inventory$embeddedbiopsy[inventory$embeddedbiopsy == "0"] <- "No"
  inventory$embeddedbiopsy[inventory$embeddedbiopsy == "1"] <- "Yes"
  
  iso <- merge(iso, inventory, by = "record_id", all = TRUE)
  # remove unneeded distribution table
  tables <- tables[names(tables) != "inventory"]
  
  tables[[which(names(tables) == "isolation")]] <- iso
  
  ### donor
  
  # have to get donor info from Pat/Joss/Aliya now
  donor <- tables$donor #read.csv(paste0(other.tables.path, "outcomes_processing_input/donor_info.csv"))

  # rework medical diagnosis
  colnames(donor)[colnames(donor) == "medicalconditions___1"] <- "T1_diabetes"
  colnames(donor)[colnames(donor) == "medicalconditions___2"] <- "T2_diabetes"

  inds.cond <- grep("medicalconditions", colnames(donor))
  donor <- donor[,-inds.cond]

  # column for HLA typing
  donor$hla_a2 <- rep(0, dim(donor)[1])
  donor$hla_a2[is.na(donor$hlaa)] <- NA

  # parse hlaa typing info
  hlaa <- donor$hlaa
  hlaa <- strsplit(hlaa, ", ")
  names(hlaa) <- donor$record_id
  hlaa <- hlaa[lapply(hlaa, length) > 1] %>% as.data.frame() %>% t() %>% as.data.frame()
  hlaa$A2 <- hlaa$V1 == "2" | hlaa$V2 == "2"
  hlaa.donors <- rownames(hlaa)[hlaa$A2]

  donor$hla_a2[donor$record_id %in% hlaa.donors] <- 1
  donor$hla_a2 <- as.character(donor$hla_a2)

  # add diagnosis summary column
  donor$diagnosis <- rep(NA, dim(donor)[1])
  donor$diagnosis[donor$T1_diabetes == "1"] <- "Type1"
  donor$diagnosis[donor$T2_diabetes == "1"] <- "Type2"
  donor$diagnosis[donor$T1_diabetes == "0" & donor$T2_diabetes == "0"] <- "None"

  # add hba1c corrected column
  donor$diagnosis_computed <- donor$diagnosis
  donor$diagnosis_computed[donor$diagnosis == "None" & as.numeric(donor$hba1c) > 6.5] <- "Type2"
  donor$diagnosis_computed[donor$diagnosis == "None" & as.numeric(donor$hba1c) < 6.5 & as.numeric(donor$hba1c) > 5.7] <- "Pre.T2D"

  # merge with distribution
  donor <- merge(donor, tables[["distribution"]], by = "record_id")

  # give factors meaningful names
  donor$donorsex[donor$donorsex == "1"] <- "Male"
  donor$donorsex[donor$donorsex == "2"] <- "Female"
  donor$donorsex[donor$donorsex == "3"] <- "Unknown"

  donor$donationtype[donor$donationtype == "1"] <- "NDD"
  donor$donationtype[donor$donationtype == "2"] <- "DCD"
  donor$donationtype[donor$donationtype == "3"] <- "MAID"

  donor$hla_a2[donor$hla_a2 == "0"] <- "Negative"
  donor$hla_a2[donor$hla_a2 == "1"] <- "Positive"

  # remove unneeded distribution table
  tables <- tables[names(tables) != "distribution"]

  # check column types
  char.cols <- c("record_id", "rrid", "donorsex", "donationtype", "hlaa", "hlab",
                 "hlabw", "hlac", "hlacw", "hladrb1", "hladr", "hladqb1", "hladqa1",
                 "hladpb1", "dpa1", "hlaother", "hla_a2", "T1_diabetes", "T2_diabetes",
                 "diagnosis", "diagnosis_computed")
  num.cols <- colnames(donor)[!(colnames(donor) %in% char.cols)]

  donor[,char.cols] <- apply(donor[,char.cols], 2, as.character)
  donor[,num.cols] <- apply(donor[,num.cols], 2, as.numeric)
  donor[donor == ""] <- NA
  
  tables[[which(names(tables) == "donor")]] <- donor
  
  ### exocytosis
  exo <- tables$exocytosis
  exo <- exo[,grep("exocytosis", colnames(exo), invert = TRUE)]
  exo <- data.table::melt(exo, id.vars = "record_id", measure.vars = colnames(exo)[-1])
  exo <- na.omit(exo)
  colnames(exo) <- c("record_id", "exposure", "insulin_exo")
  exo$cell_id <- NA
  
  don.uniq <- unique(exo$record_id)
  for(i in c(1:length(don.uniq))){
    exo[exo$record_id == don.uniq[i], ]$cell_id <- paste0(don.uniq[i], "_", 1:sum(exo$record_id == don.uniq[i]))
  }
  
  exo <- exo[order(exo$record_id), ]
  rownames(exo) <- NULL
  exo <- exo[,c(1,4,2,3)]
  exo$exposure <- as.character(exo$exposure)
  exo$exposure[grep("lgexo", exo$exposure)] <- "gluc_1"
  exo$exposure[grep("mgexo", exo$exposure)] <- "gluc_5"
  exo$exposure[grep("hgexo", exo$exposure)] <- "gluc_10"
  exo$insulin_exo <- as.numeric(exo$insulin_exo)
  
  # set exocytosis values below zero to zero (match patchseq processing)
  exo$insulin_exo[exo$insulin_exo < 0] <- 0
  
  tables[[which(names(tables) == "exocytosis")]] <- exo
  
  ### reformat gsis
  gsis <- tables[["gsis"]]
  gsis <- as.data.table(gsis)
  gsis.ct <- gsis[,c("record_id", "culturetime2")]
  gsis <- data.table::melt(gsis, id.vars = 'record_id', 
                           measure.vars = list(
                             c("lowmedcontent12", "lowmedcontent22", "lowmedcontent32"), 
                             c("lowglucose12_60a3b7", "lowglucose22", "lowglucose32"),
                             c("medglucose12", "medglucose22", "medglucose32"),
                             c("highcontent12", "highcontent22", "highcontent32"),
                             c("lowglucose72", "lowglucose82", "lowglucose92"),
                             c("highglucose12", "highglucose22", "highglucose32"),
                             c("content28167_1", "content28167_2", "content28167_3"),
                             c("value28_1", "value28_2", "value28_3"),
                             c("value167_1", "value167_2", "value167_3")),
                           value.name = c('insulin.content.1', 'first.1', 'second.1',
                                          'insulin.content.2', 'first.2', 'second.2',
                                          'insulin.content.3', 'first.3', 'second.3'),
                           variable.name = "replicate")
  
  # split into three
  gsis1 <- gsis[,c(1:5)]
  gsis2 <- gsis[,c(1,2,6:8)]
  gsis3 <- gsis[,c(1,2,9:11)]
  
  gsis1 <- gsis1[apply(gsis1, 1, function(x){sum(is.na(x)) != 3}),]
  gsis2 <- gsis2[apply(gsis2, 1, function(x){sum(is.na(x)) != 3}),]
  gsis3 <- gsis3[apply(gsis3, 1, function(x){sum(is.na(x)) != 3}),]
  
  colnames(gsis1)[3:5] <- c("total_insulin_content", "first_insulin_secretion", "second_insulin_secretion")
  colnames(gsis2)[3:5] <- c("total_insulin_content", "first_insulin_secretion", "second_insulin_secretion")
  colnames(gsis3)[3:5] <- c("total_insulin_content", "first_insulin_secretion", "second_insulin_secretion")
  
  gsis1$first_gluc_conc <- 1
  gsis1$second_gluc_conc <- 10
  gsis2$first_gluc_conc <- 1
  gsis2$second_gluc_conc <- 16.7
  gsis3$first_gluc_conc <- 2.8
  gsis3$second_gluc_conc <- 16.7
  
  gsis <- rbind(gsis1, gsis2)
  gsis <- rbind(gsis, gsis3)
  gsis <- gsis[,c(1:2,6:7,3:5)]
  gsis <- gsis[order(gsis$record_id),]
  
  gsis$rows <- rownames(gsis)
  gsis$replicate <- gsis[ ,rank(rows), record_id]$V1
  gsis <- as.data.frame(gsis)
  gsis <- gsis[,colnames(gsis) != "rows"]
  
  # convert columns to correct type
  char.cols <- c("record_id", "replicate")
  num.cols <- colnames(gsis)[!(colnames(gsis) %in% char.cols)]
  gsis[,char.cols] <- apply(gsis[,char.cols], 2, as.character)
  gsis[,num.cols] <- apply(gsis[,num.cols], 2, as.numeric)
  
  # Compute composite variables
  gsis$first_insulin_percent <- (gsis$first_insulin_secretion/gsis$total_insulin_content)*100
  gsis$second_insulin_percent <- (gsis$second_insulin_secretion/gsis$total_insulin_content)*100
  
  gsis$gluc_conc_group <- paste0(gsis$first_gluc_conc, "_", gsis$second_gluc_conc)
  gsis$stim_index <- gsis$second_insulin_secretion/gsis$first_insulin_secretion
  
  gsis <- merge(gsis, gsis.ct, by = "record_id")
  
  tables[[which(names(tables) == "gsis")]] <- gsis
  
  ### electrophysiology
  ephys <- tables$ephys
  
  # get into better format
  ephys <- ephys[ephys$epnumber > 0, ]
  ephys <- data.table::melt(ephys, id.vars = "record_id", measure.vars = colnames(ephys)[-1])
  ephys <- na.omit(ephys)
  ephys$variable <- as.character(ephys$variable)
  ephys <- ephys[ephys$variable != "epnumber", ]
  ephys <- ephys[ephys$variable != "electrophysiology_complete", ]
  ephys$cell_id <- gsub("[^0-9.-]", "", ephys$variable)
  ephys$variable <- gsub("[0-9.+]", "", ephys$variable)
  ephys <- dcast(ephys, record_id+cell_id ~ variable, value.var = "value")
  
  # filter out unneeded cell types
  ephys <- ephys[ephys$eptype %in% c("1", "2"), ]
  ephys$eptype[ephys$eptype == "1"] <- "Beta"
  ephys$eptype[ephys$eptype == "2"] <- "Alpha" 
  
  # remove unneeded columns
  ephys <- ephys[,!c(colnames(ephys) %in% c("epnscc", "eprpbr", "epvspc", "epname"))]
  
  # re-order columns
  ephys <- ephys[,c("record_id", "cell_id", "eptype", "epdays", "epsize", "epntc", "epnfdc", "epnldc", "nci", "epnpsca", "ephisc", "epnepcca", "epnlcca")]
  ephys[,4:13] <- apply(ephys[,4:13], 2, as.numeric)
  colnames(ephys) <- c("record_id", "cell_id", "eptype", "epdays", "cell_size_pF", "total_exocytosis_fF_pF", "early_exocytosis_fF_pF", "late_exocytosis_fF_pF", "calcium_entry_pC_pF", "na_current_amp_pA_pF", "na_half_inactivation_mV", "early_ca_current_pA_pF", "late_ca_current_pA_pF")
  
  # set exocytosis values below zero to zero
  ephys$total_exocytosis_fF_pF[ephys$total_exocytosis_fF_pF < 0] <- 0
  ephys$early_exocytosis_fF_pF[ephys$early_exocytosis_fF_pF < 0] <- 0
  ephys$late_exocytosis_fF_pF[ephys$late_exocytosis_fF_pF < 0] <- 0
  
  # flip sign of the calcium and sodium variables
  ephys$calcium_entry_pC_pF <- -1*ephys$calcium_entry_pC_pF
  ephys$na_current_amp_pA_pF <- -1*ephys$na_current_amp_pA_pF
  ephys$early_ca_current_pA_pF <- -1*ephys$early_ca_current_pA_pF
  ephys$late_ca_current_pA_pF <- -1*ephys$late_ca_current_pA_pF
  
  tables[[which(names(tables) == "ephys")]] <- ephys
  
  ### inventory
  inventory <- tables$inventory
  
  ### seahorse
  mito <- read.csv("seahorse.csv", skip=1)
  mito.norm <- read.csv("seahorse_norm.csv", skip = 1)
  
  # process first df
  mito <- as.data.frame(mito)
  colnames(mito)[1] <- "X"
  blank.inds <- which(mito$X == "")
  if(length(blank.inds) > 0){mito <- mito[-blank.inds,]}
  rownames(mito) <- mito$X
  mito <- mito[,-1]
  
  mito <- t(mito) %>% as.data.frame()
  mito$record_id <- substr(rownames(mito), 1, 4)
  mito$replicate <- substr(rownames(mito), 6, 6)
  rownames(mito) <- NULL
  mito <- mito[,c(45:46,1:44)]
  
  # process second df
  mito.norm <- as.data.frame(mito.norm)
  colnames(mito.norm)[1] <- "X"
  blank.inds <- which(mito.norm$X == "")
  if(length(blank.inds) > 0){mito.norm <- mito.norm[-blank.inds,]}
  rownames(mito.norm) <- mito.norm$X
  mito.norm <- mito.norm[,-1]
  
  mito.norm <- t(mito.norm) %>% as.data.frame()
  mito.norm$record_id <- substr(rownames(mito.norm), 1, 4)
  mito.norm$replicate <- substr(rownames(mito.norm), 6, 6)
  rownames(mito.norm) <- NULL
  mito.norm <- mito.norm[,c(45:46,1:44)]
  
  proc.mito <- mito.norm[,colnames(mito.norm) != "calc_basal_resp"]
  proc.mito <- merge(proc.mito, mito[,c("record_id", "replicate", "calc_basal_resp")], by = c("record_id", "replicate"))
  proc.mito <- proc.mito[order(proc.mito$record_id), ]
  
  tables[["seahorse_norm_dna"]] <- mito
  tables[["seahorse_norm_dna_baselineoc"]] <- mito.norm
  tables[["seahorse"]] <- proc.mito
  
  ### do the same for perifusion once it's set up
  gluc.nd <- read_pzfx("perifusion.pzfx", table = 1)
  gluc.t2d <- read_pzfx("perifusion.pzfx", table = 4)
  
  leu.nd <- read_pzfx("perifusion.pzfx", table = 2)
  leu.t2d <- read_pzfx("perifusion.pzfx", table = 5)
  
  olp.nd <- read_pzfx("perifusion.pzfx", table = 3)
  olp.t2d <- read_pzfx("perifusion.pzfx", table = 6)
  
  # process tables
  gluc <- cbind(gluc.nd, gluc.t2d[,-1])
  rownames(gluc) <- paste0("time_", gluc$`Time (min)`)
  gluc <- gluc[,-1] %>% t() %>% as.data.frame()
  donor <- data.frame(record_id = gsub("_.*", "", rownames(gluc)), replicate = gsub("R..._", "", rownames(gluc)))
  donor$replicate <- gsub("R...._", "", donor$replicate)
  gluc <- cbind(donor, gluc)
  rownames(gluc) <- NULL
  gluc <- gluc[order(gluc$record_id), ]
  tables[["peri_gluc"]] <- gluc
  
  leu <- cbind(leu.nd, leu.t2d[,-1])
  rownames(leu) <- paste0("time_", leu$`Time (min)`)
  leu <- leu[,-1] %>% t() %>% as.data.frame()
  donor <- data.frame(record_id = gsub("_.*", "", rownames(leu)), replicate = gsub("R..._", "", rownames(leu)))
  donor$replicate <- gsub("R...._", "", donor$replicate)
  leu <- cbind(donor, leu)
  rownames(leu) <- NULL
  leu <- leu[order(leu$record_id), ]
  tables[["peri_leu"]] <- leu
  
  olp.nd <- olp.nd[!is.na(olp.nd$`Time (min)`),]
  olp.t2d <- olp.t2d[!is.na(olp.t2d$`Time (min)`),]
  olp <- cbind(olp.nd, olp.t2d[,-1])
  rownames(olp) <- paste0("time_", olp$`Time (min)`)
  olp <- olp[,-1] %>% t() %>% as.data.frame()
  donor <- data.frame(record_id = gsub("_.*", "", rownames(olp)), replicate = gsub("R..._", "", rownames(olp)))
  donor$replicate <- gsub("R...._", "", donor$replicate)
  olp <- cbind(donor, olp)
  rownames(olp) <- NULL
  olp <- olp[order(olp$record_id), ]
  tables[["peri_olp"]] <- olp
  
  ### Add additional summary tables
  raw_variable_summary <- read.csv(paste0(other.tables.path, "outcomes_processing_input/raw_variable_summary.csv"));
  function_summary <- read.csv(paste0(other.tables.path, "outcomes_processing_input/function_summary.csv"));
  disc_groups <- read.csv(paste0(other.tables.path, "display_interface/disc_groups.csv"));
  proc_variable_summary <- read.csv(paste0(other.tables.path, "display_interface/proc_variable_summary.csv"));
  
  patchseq <- readRDS(paste0(other.tables.path, "outcomes_processing_input/patchseq_metadata.rds"));
  
  ### combine ephys data
  ephys.rc <- tables[["ephys"]]
  ephys.ps <- patchseq
  ephys.exo <- tables[["exocytosis"]]
  
  # make columns uniform
  ephys.rc <- ephys.rc[,-c(4)] # remove cell ID and days in culture
  ephys.rc$cell_id <- paste0(ephys.rc$record_id, "_", ephys.rc$cell_id, "_rc")
  ephys.rc$glucose_mM <- "5"
  colnames(ephys.rc)[3] <- "cell_type"
  
  ephys.exo$cell_id <- paste0(ephys.exo$cell_id, "_exo")
  ephys.exo$exposure <- gsub("gluc_", "", ephys.exo$exposure)
  ephys.exo$cell_type <- "Beta"
  colnames(ephys.exo)[4] <- "total_exocytosis_fF_pF"
  colnames(ephys.exo)[3] <- "glucose_mM"
  
  # merge by columns
  ephys.all <- merge(ephys.rc, ephys.ps, by = colnames(ephys.rc), all = TRUE)
  ephys.all <- merge(ephys.all, ephys.exo, by = colnames(ephys.exo), all = TRUE)
  
  # remove non-UofA donors
  ephys.all <- ephys.all[grep("R", ephys.all$record_id), ]
  ephys.all <- ephys.all[,c(1,2,5,3,6,4,7:13)]
  ephys.all[,c(5:13)] <- apply(ephys.all[,5:13], 2, as.numeric)
  ephys.dt <- as.data.table(ephys.all[,-c(2)])
  ephys.dt <- ephys.dt[,lapply(.SD, mean), by = .(record_id, cell_type, glucose_mM)]
  ephys.dt <- as.data.frame(ephys.dt)
  colnames(ephys.dt)[4:12] <- paste0(colnames(ephys.dt)[4:12], "_donor")
  colnames(patchseq)[5:13] <- paste0(colnames(patchseq)[5:13], "_cell")
  
  tables[["ephys_donor"]] <- ephys.dt
  tables[["exocytosis"]] <- NULL
  tables[["ephys"]] <- NULL
  
  #### Add computed values
  computed <- read.csv(paste0(other.tables.path, "outcomes_processing_input/composition.csv"));
  tables[["computed"]] <- computed;
  
  tables[["function_summary"]] <- function_summary;
  tables[["proc_variable_summary"]] <- proc_variable_summary;
  tables[["ephys_cell"]] <- patchseq
  
  ### Process metadata variables
  # we don't add ephys here because it is single-cell (patchseq) or there are multiple outcomes per donor (gluc, cell type)
  # patchseq also does not need transformation because of spearman correlation
  metaVar <- raw_variable_summary
  uniq.tables <- unique(metaVar$table)
  varValues <- data.frame()
  for(i in c(1:length(uniq.tables))){
    dat <- tables[[uniq.tables[i]]]
    dat <- dat[,c("record_id", metaVar$column[metaVar$table == uniq.tables[i]])]
    
    if(uniq.tables[i] == "donor"){
      varValues <- dat
    } else if(uniq.tables[i] == "isolation"){
      dat$cryotubesremaining <- as.numeric(dat$cryotubesremaining)
      dat$sftubesremaining <- as.numeric(dat$sftubesremaining)
      varValues <- merge(varValues, dat, by = "record_id", all = TRUE)
    } else if(uniq.tables[i] == "seahorse"){
      dat <- aggregate(dat[which(colnames(dat) != "record_id")], dat[which(colnames(dat) == "record_id")], mean)
      varValues <- merge(varValues, dat, by = "record_id", all = TRUE)
    } else if(uniq.tables[i] == "gsis"){
      first.secr <- dat[,c("record_id", "first_gluc_conc", "first_insulin_secretion")]
      second.secr <- dat[,c("record_id", "second_gluc_conc", "second_insulin_secretion")]
      ct <- distinct(dat[,c("record_id", "culturetime2")])
      colnames(first.secr) <- c("record_id", "gluc_conc", "insulin_secretion")
      colnames(second.secr) <- c("record_id", "gluc_conc", "insulin_secretion")
      ins.secr <- rbind(first.secr, second.secr)
      ins.secr$gluc_conc <- paste0("insulin_secretion_", gsub("\\.", "p", ins.secr$gluc_conc))
      
      med.secr <- aggregate(ins.secr[c("insulin_secretion")], ins.secr[c("record_id", "gluc_conc")], median, na.rm = TRUE)
      med.secr <- dcast(med.secr, record_id ~ gluc_conc, value.var = "insulin_secretion")
      
      med.stim <- aggregate(dat[c("stim_index")], dat[c("record_id", "gluc_conc_group")], median, na.rm = TRUE)
      med.stim$gluc_conc_group <- paste0("gsis_index_", gsub("\\.", "p", med.stim$gluc_conc_group))
      med.stim <- dcast(med.stim, record_id ~ gluc_conc_group, value.var = "stim_index")
      
      med.ins <- aggregate(dat[c("total_insulin_content")], dat[c("record_id")], median, na.rm = TRUE)
      
      dat <- merge(med.secr, med.stim, by = "record_id", all = TRUE)
      dat <- merge(dat, med.ins, by = "record_id", all = TRUE)
      dat <- merge(dat, ct, by = "record_id", all = TRUE)
      dat <- dat %>% mutate_all(~ifelse(is.nan(.), NA, .))
      dat <- dat %>% mutate_all(~ifelse(is.infinite(.), NA, .))
      dat$culturetime2 <- as.numeric(dat$culturetime2)
      
      # detect and remove donors that are four standard deviations beyond the mean (after log transformation)
      gsis.out <- dat[, c(colnames(dat) != "culturetime2")]
      rownames(gsis.out) <- dat$record_id
      gsis.out <- gsis.out[,-1]
      gsis.outlier <- apply(gsis.out, 2, function(x){
        log.x <- log10(x)
        col.mean <- mean(log.x, na.rm = T)
        col.sd <- sd(log.x, na.rm = T)
        inds <- which(log.x > col.mean + (4*col.sd) | log.x < col.mean - (4*col.sd))
      })
      for(i in c(1:length(gsis.outlier))){
        dat[gsis.outlier[[i]], names(gsis.outlier)[i]] <- NA
      }
      
      varValues <- merge(varValues, dat, by = "record_id", all = TRUE)
    } else if(uniq.tables[i] == "function_summary"){
      varValues <- merge(varValues, dat, by = "record_id", all = TRUE)
    } else if (uniq.tables[i] == "computed"){
      varValues <- merge(varValues, dat, by = "record_id", all = TRUE)
    }
  }
  
  # transform metadata
  right.skew <- c("hba1c", "pdinsulinperieq", "pdinsulindnaratio", "insulindnaratio", 
                  "auc_baseline_3mmgluc", "peak_gluc_15mmgluc", "auc_gluc_15mmgluc", "ieqperpancreasweight", "insulincontent",
                  "auc_gluc_6mmgluc", "peak_leu_5mmleu", "auc_leu_5mmleu", "auc_leu_5mmleu_6mmgluc", "peak_olp_1p5mmolp",
                  "gsis_index_1_10", "gsis_index_1_16p7", "gsis_index_2p8_16p7", "total_insulin_content",
                  "insulin_secretion_1", "insulin_secretion_10", "insulin_secretion_16p7", "insulin_secretion_2p8",
                  "dnacontent", "insulinperieq", "pdinsulincontent")
  
  for(i in c(1:length(right.skew))){
    varValues[,right.skew[i]] <- log10(varValues[,right.skew[i]])
  }
  
  neg.skew.funct <- c("auc_gluc_30mmkcl", "auc_leu_30mmkcl", "auc_olp_1p5mmolp", "auc_olp_1p5mmolp_6mmgluc", "auc_olp_30mmkcl")
  for(i in c(1:length(neg.skew.funct))){
    varValues[,neg.skew.funct[i]] <- log10(varValues[,neg.skew.funct[i]] + abs(min(varValues[,neg.skew.funct[i]], na.rm=TRUE)) + 1000)
  }
  
  varValues[varValues == Inf] <- NA
  varValues <- varValues %>% mutate_all(~ifelse(is.nan(.), NA, .))
  varValues[varValues == -Inf] <- NA
  
  tables[["proc_metadata"]] <- varValues;
  
  ### Write out data we need for donor view page
  donor.info <- varValues[,c("record_id", "donorage", "bodymassindex", "hba1c")]
  donor.info$hba1c <- 10^donor.info$hba1c
  write.csv(donor.info, paste0(other.tables.path, "display_data/numerical_donor_info.csv"), row.names = FALSE)
  
  ### Write out isolation data for donor view page
  iso.vars <- c("record_id", "coldischemiatime", "puritypercentage", "trappedpercentage", "pancreasweight",
                "digesttime", "totalieq", "isletparticleindex", "ieqperpancreasweight", "insulincontent", 
                "insulinperieq", "predistributionculturetime", "percentieqrecoverypostculture", "pdisletparticleindex",
                "pdinsulinperieq", "cryotubesremaining", "sftubesremaining")
  
  iso.info <- varValues[, iso.vars]
  iso.info[, iso.vars[iso.vars %in% right.skew]] <- apply(iso.info[, iso.vars[iso.vars %in% right.skew]], 2, function(x){10^x})
  iso.info$insulincontent <- iso.info$insulincontent/1000
  write.csv(iso.info, paste0(other.tables.path, "display_data/numerical_isolation_info.csv"), row.names = FALSE)
  
  ### Write out GSIS data for donor view page
  gsis.vars <- c("record_id", "culturetime2", "total_insulin_content", "insulin_secretion_1", "insulin_secretion_2p8", 
                 "insulin_secretion_10", "insulin_secretion_16p7", "gsis_index_1_10", "gsis_index_1_16p7", "gsis_index_2p8_16p7")
  gsis.info <- varValues[, gsis.vars]
  gsis.info$culturetime2 <- as.numeric(gsis.info$culturetime2)
  gsis.info[, gsis.vars[gsis.vars %in% right.skew]] <- apply(gsis.info[, gsis.vars[gsis.vars %in% right.skew]], 2, function(x){10^x})
  # RedCap says values are pg/mL, we want ug/mL, so divide by 1,000,000
  gsis.info$total_insulin_content <- gsis.info$total_insulin_content/1000000
  gsis.info$insulin_secretion_1 <- gsis.info$insulin_secretion_1/1000000
  gsis.info$insulin_secretion_10 <- gsis.info$insulin_secretion_10/1000000
  gsis.info$insulin_secretion_16p7 <- gsis.info$insulin_secretion_16p7/1000000
  gsis.info$insulin_secretion_2p8 <- gsis.info$insulin_secretion_2p8/1000000
  write.csv(gsis.info, paste0(other.tables.path, "display_data/numerical_gsis.csv"), row.names = FALSE)
  
  ### write out mito func time series data
  seahorse <- read.csv("seahorse.csv", skip = 1)
  seahorse <- as.data.frame(seahorse)
  colnames(seahorse)[1] <- "X"
  blank.inds <- which(seahorse$X == "")
  if(length(blank.inds) > 0){seahorse <- seahorse[-blank.inds,]}
  rownames(seahorse) <- seahorse$X
  seahorse <- seahorse[,-1]
  
  seahorse <- t(seahorse) %>% as.data.frame()
  seahorse$record_id <- substr(rownames(seahorse), 1, 4)
  seahorse$replicate <- substr(rownames(seahorse), 6, 6)
  rownames(seahorse) <- NULL
  seahorse <- seahorse[,c(45:46,1:44)]
  
  mito.cols <- grep("time", colnames(seahorse))
  mito.time <- seahorse[,c(1,mito.cols)]
  mito.time <- reshape2::melt(mito.time)
  colnames(mito.time) <- c("record_id", "time", "insulin")
  mito.time$time <- gsub("time_", "", mito.time$time)
  mito.time$time <- as.numeric(mito.time$time)
  mito.time <- as.data.table(mito.time)
  mito.time <- mito.time[,.(insulin = mean(insulin)), by = .(record_id, time)]
  mito.time <- mito.time[,c("record_id", "insulin", "time")]
  write.csv(mito.time, paste0(other.tables.path, "display_data/mito_func.csv"), row.names = FALSE)
  
  ### write out normalized mito func time series data
  mito.norm <- tables[["seahorse"]] 
  mito.cols <- grep("time", colnames(mito.norm))
  mito.time <- mito.norm[,c(1,mito.cols)]
  mito.time <- reshape2::melt(mito.time)
  colnames(mito.time) <- c("record_id", "time", "insulin")
  mito.time$time <- gsub("time_", "", mito.time$time)
  mito.time$time <- as.numeric(mito.time$time)
  mito.time <- as.data.table(mito.time)
  mito.time <- mito.time[,.(insulin = mean(insulin)), by = .(record_id, time)]
  mito.time <- mito.time[,c("record_id", "insulin", "time")]
  write.csv(mito.time, paste0(other.tables.path, "display_data/mito_func_norm.csv"), row.names = FALSE)
  
  ### write out mito numerical values based on normalized matrix
  mito.outcomes <- proc.mito[,c("record_id", "calc_nonmito_oc", "calc_basal_resp", "calc_atp_resp", "calc_proton_leak", "calc_max_gluc_resp", 
                                "calc_stim_gluc_resp", "calc_max_resp", "calc_spare_cap", "calc_gluc_stim_oci", "calc_gluc_stim_sparecap")]
  mito.outcomes <- aggregate(mito.outcomes[,2:11], by = list(mito.outcomes$record_id), FUN = mean)
  colnames(mito.outcomes)[1] <- "record_id"
  write.csv(mito.outcomes, paste0(other.tables.path, "display_data/numerical_mito.csv"), row.names = FALSE)
  
  ### write out time series perifusion data
  gluc2 <- reshape2::melt(gluc)
  gluc2 <- na.omit(gluc2)
  gluc2 <- as.data.table(gluc2)
  gluc2 <- gluc2[,.(insulin = mean(value)), by = .(record_id, variable)]
  gluc2$variable <- as.numeric(gsub("time_", "", gluc2$variable))
  colnames(gluc2) <- c("record_id", "time", "insulin")
  gluc2 <- gluc2[,c("record_id", "insulin", "time")]
  gluc.bl <- gluc2[time < 25, .(baseline = mean(insulin)), by = .(record_id)]
  gluc.norm <- merge(gluc2, gluc.bl, by = "record_id")
  gluc.norm$insulin_norm <- gluc.norm$insulin/gluc.norm$baseline
  write.csv(gluc.norm, paste0(other.tables.path, "display_data/peri_gluc.csv"), row.names = FALSE)
  
  leu2 <- reshape2::melt(leu)
  leu2 <- na.omit(leu2)
  leu2 <- as.data.table(leu2)
  leu2 <- leu2[,.(insulin = mean(value)), by = .(record_id, variable)]
  leu2$variable <- as.numeric(gsub("time_", "", leu2$variable))
  colnames(leu2) <- c("record_id", "time", "insulin")
  leu2 <- leu2[,c("record_id", "insulin", "time")]
  leu.bl <- leu2[time < 25, .(baseline = mean(insulin)), by = .(record_id)]
  leu.norm <- merge(leu2, leu.bl, by = "record_id")
  leu.norm$insulin_norm <- leu.norm$insulin/leu.norm$baseline
  write.csv(leu.norm, paste0(other.tables.path, "display_data/peri_leu.csv"), row.names = FALSE)
  
  olp2 <- reshape2::melt(olp)
  olp2 <- na.omit(olp2)
  olp2 <- as.data.table(olp2)
  olp2 <- olp2[,.(insulin = mean(value)), by = .(record_id, variable)]
  olp2$variable <- as.numeric(gsub("time_", "", olp2$variable))
  colnames(olp2) <- c("record_id", "time", "insulin")
  olp2 <- olp2[,c("record_id", "insulin", "time")]
  olp.bl <- olp2[time < 25, .(baseline = mean(insulin)), by = .(record_id)]
  olp.norm <- merge(olp2, olp.bl, by = "record_id")
  olp.norm$insulin_norm <- olp.norm$insulin/olp.norm$baseline
  write.csv(olp.norm, paste0(other.tables.path, "display_data/peri_olp.csv"), row.names = FALSE)
  
  ### write out ephys numerical values
  write.csv(ephys.dt, paste0(other.tables.path, "display_data/numerical_ephys.csv"), row.names = FALSE)
  
  colnames(ephys.all)[5:13] <- paste0(colnames(ephys.all)[5:13], "_donor")
  write.csv(ephys.all, paste0(other.tables.path, "display_data/numerical_ephys_cell.csv"), row.names = FALSE)
  
  ### download images
  dat <- redcap_report(redcap_uri = api_url,
                       token = api_token, 
                       report_id = 40576,
                       guess_type = FALSE)
  
  dat <- dat$data
  
  # there are two types: post-isolation and pre-distribution, with low-res and high-res
  # we want low res of each type
  dat <- dat[,c("record_id", "dithizonestaining", "dithizonestaining3")]
  
  dat$dithizonestaining[grep("no image available", dat$dithizonestaining)] <- NA
  dat$dithizonestaining3[grep("no image available", dat$dithizonestaining3)] <- NA
  dat$dithizonestaining[grep("tif", dat$dithizonestaining)] <- NA
  dat$dithizonestaining3[grep("tif", dat$dithizonestaining3)] <- NA
  
  postiso_path = paste0(other.tables.path, "images/islets/postiso")
  predist_path = paste0(other.tables.path, "images/islets/predist")
  
  records <- dat$record_id
  for(i in c(1:length(records))){
    
    if(!is.na(dat[i,"dithizonestaining"])){ # check if image in RedCap
      imgNm <- paste0(records[i], "_postiso_32.png")
      
      if(!file.exists(paste0(postiso_path, "/", imgNm))){ # check if image previously downloaded
        redcap_download_file_oneshot(
          file_name = imgNm,
          directory = postiso_path,
          overwrite = TRUE,
          redcap_uri = api_url,
          token = api_token,
          record = records[i],
          field = "dithizonestaining"
        )
      }
    }
    
    if(!is.na(dat[i,"dithizonestaining3"])){
      imgNm <- paste0(records[i], "_predist_32.png")
      
      if(!file.exists(paste0(predist_path, "/", imgNm))){
        redcap_download_file_oneshot(
          file_name = paste0(records[i], "_predist_32.png"),
          directory = predist_path,
          overwrite = TRUE,
          redcap_uri = api_url,
          token = api_token,
          record = records[i],
          field = "dithizonestaining3"
        )
      }
      
    }
  }
  
  #### Write out to sqlite ####
  sqlite.tables <- c("computed", "donor", "ephys_cell", "ephys_donor", "gsis", 
                     "isolation", "proc_metadata", "proc_variable_summary", "seahorse",
                     "seahorse_norm_dna", "seahorse_norm_dna_baselineoc",
                     "function_summary",'peri_gluc', 'peri_leu', 'peri_olp')
  con <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  for(i in c(1:length(sqlite.tables))){
    dat <- tables[[sqlite.tables[i]]]
    table.name <- sqlite.tables[i]
    
    dbWriteTable(con, table.name, dat, overwrite = TRUE)
  }
  dbDisconnect(con)
  
}
