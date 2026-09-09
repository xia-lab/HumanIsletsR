# R Functions for HumanIslets web tool
# Author: Jessica Ewald, Yao Lu

################################################################################

getDonors_fun <- function(version, donorListKey = ""){
  require(RSQLite)
  require(dplyr)
  require(data.table)
  require(stringr)
  require(jsonlite)

  # set sqlite file path
if(version=="v2"){
 omics.path <- paste0(sqlite.path, "HI_omics_v2.sqlite");
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite");
}else{
 omics.path <- paste0(sqlite.path, "HI_omics.sqlite");
  outcomes.path <- paste0(sqlite.path, "HI_tables_v1.sqlite");

}
 
  # parse filters
  filter <- read_json("filter.json");
  
  if(filter$omics == "bulk"){
    # get donor information table
    con <- dbConnect(RSQLite::SQLite(), outcomes.path)
    donor.info <- dbReadTable(con, "donor")
    iso.info <- dbReadTable(con, "isolation")
    prohormone.info <- dbReadTable(con, "prohormone")
    dbDisconnect(con)
    donor.table <- as.data.table(donor.info)
    iso.table <- as.data.table(iso.info)
    prohormone.table <-  as.data.table(prohormone.info)
    lipid.table <- NULL
  } else if(filter$omics == "singlecell"){
    if(version=="v2"){
    donor.info <- readRDS(paste0(other.tables.path, "analysis_input/patchseq_donors_v2.rds"));
    }else{
    donor.info <- readRDS(paste0(other.tables.path, "analysis_input/patchseq_donors.rds"));
   }

    donor.table <- as.data.table(donor.info);
  }

  filter <- filter[!(names(filter) %in% c("usrDir", "omics"))]
  filter <- lapply(filter, function(x){str_split(x, "\\|")[[1]]}) %>% as.data.frame() %>% t() %>% as.data.frame();
  filter <- filter[filter$V1 == "true", ]

  filter.vars <- rownames(filter)
  filter <- as.list(filter$V2)
  names(filter) <- filter.vars
 print(filter)
  donors <- donor.info$record_id
  
  # for loop with switch statement inside to progressively filter donor list
  if (length(filter.vars) > 0) {
    for (i in c(1:length(filter.vars))) {
      if (length(donors) > 0) {
        
        filter.var <- filter.vars[i];
        
        switch (filter.var,
                
               avail = {
                 filt <- str_split(filter$avail, ",")[[1]];
                 if (filter$avail != "") {

                   # HI_tables tables
                   filt.outcomes <- filt[filt %in% c('ephys_donor', 'gsis', 'genes', 'seahorse', 'peri_gluc', 'peri_leu', 'peri_olp', 'ephys_cell')]
                   con <- dbConnect(RSQLite::SQLite(), outcomes.path)
                   filt.list <- c()
                   if(length(filt.outcomes) > 0){
                   for (j in c(1:length(filt.outcomes))) {
                     if (length(donors) > 1){
                       table <- dbReadTable(con, name = filt.outcomes[j])
   
                       switch (filt.outcomes[j],
                               ephys_donor = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               ephys_cell = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               gsis = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               seahorse = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               peri_gluc = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               peri_leu = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               peri_olp = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               prohormone = {
                                 filt.list <- unique(table$record_id);
                                 donors <- intersect(donors, filt.list);
                               },
                               { print("table name not found"); });
                     } else {
                       break;
                     }
                   }
                     }
                   dbDisconnect(con);

                   # Omics tables
                   filt.omics <- filt[filt %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta','proc_contaminants')]
                   if(length(filt.omics) > 0){
                     con <- dbConnect(RSQLite::SQLite(), omics.path)
                     filt.list <- c()
                     for (j in c(1:length(filt.omics))) {
                       if (length(donors) > 1){
                         # v2 proteomics spans two acquisition batches; the donor
                         # cohort is the merged proc_prot_combine matrix (234
                         # donors). v1 proc_prot is batch 1 only (134 donors) and
                         # is the sole proteomics table in HI_omics.sqlite, hence
                         # the version gate -- without it "has proteomics data"
                         # would silently cap every selection at the old 134.
                         omics.nm <- filt.omics[j]
                         if(omics.nm == "proc_prot" && version == "v2" && dbExistsTable(con, "proc_prot_combine")){
                           omics.nm <- "proc_prot_combine"
                         }
                         table <- dbReadTable(con, name = omics.nm)
                         switch (filt.omics[j],
                                 proc_nanostring_merge = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 proc_prot = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 proc_rnaseq = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 proc_pbrna_Alpha = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 proc_pbrna_Beta = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 proc_contaminants = {
                                   donors <- intersect(donors, colnames(table));
                                 },
                                 { print("table name not found"); });
                       } else {
                         break;
                       }
                     }
                     dbDisconnect(con);
                   };

                   # CSV-based tables (lipid extraction)
                   if('lipid' %in% filt){
                     if(is.null(lipid.table)){
                       lipid.table <- as.data.table(read.csv(paste0(other.tables.path, "display_data/lip_extract.csv")))
                     }
                     lipid.donors <- unique(lipid.table$record_id);
                     donors <- intersect(donors, lipid.donors);
                   }

                 }
               },
               
               sex = {
                 filt <- str_split(filter$sex, ",")[[1]];
                 filt.list <- donor.table[(donorsex %in% filt), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               diagnosis = {
                 filt <- str_split(filter$diagnosis, ",")[[1]];
                 if (sum(filt %in% c("Type1", "Type2", "None")) != 3) {
                   
                   filt.list <- c();
                   
                   if ("Type1" %in% filt) {
                     temp <- donor.table[diagnosis == "Type1", record_id];
                     filt.list <- unique(c(filt.list, temp));
                   }
                   
                   if ("Type2" %in% filt) {
                     temp <- donor.table[diagnosis == "Type2", record_id];
                     filt.list <- unique(c(filt.list, temp));
                   }
                   
                   if ("None" %in% filt) {
                     temp <- donor.table[diagnosis == "None", record_id];
                     filt.list <- unique(c(filt.list, temp));
                   }
                 }
                 donors <- intersect(donors, filt.list);
               },
               
               age = {
                 filt <- str_split(filter$age, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(donorage > filt[1] & donorage < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               hba1c = {
                 filt <- str_split(filter$hba1c, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(hba1c > filt[1] & hba1c < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               bmi = {
                 filt <- str_split(filter$bmi, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(bodymassindex > filt[1] & bodymassindex < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               digTime = {
                 filt <- str_split(filter$digTime, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(digesttime > filt[1] & digesttime < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               hla2Type = {
                 filt <- str_split(filter$hla2Type, ",")[[1]];
                 filt.list <- donor.table[(hla_a2 %in% filt), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               distInsPerIEQ = {
                 filt <- str_split(filter$distInsPerIEQ, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(pdinsulinperieq > filt[1] & pdinsulinperieq < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               cit = {
                 filt <- str_split(filter$cit, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(coldischemiatime > filt[1] & coldischemiatime < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               isoInsDNA = {
                 filt <- str_split(filter$isoInsDNA, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(insulindnaratio > filt[1] & insulindnaratio < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               corDiag = {
                filt <- str_split(filter$corDiag, ",")[[1]];
                filt.list <- donor.table[(diagnosis_computed %in% filt), record_id];
                donors <- intersect(donors, filt.list);
               },

               cluster = {
                filt <- str_split(filter$cluster, ",")[[1]];
                filt.list <- donor.table[(final_cluster %in% filt), record_id];
                donors <- intersect(donors, filt.list);
               },

               donType = {
                 filt <- str_split(filter$donType, ",")[[1]];
                 filt.list <- donor.table[(donationtype %in% filt), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               isletParticles = {
                 filt <- str_split(filter$isletParticles, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(isletparticles > filt[1] & isletparticles < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               distInsDNA = {
                 filt <- str_split(filter$distInsDNA, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(pdinsulindnaratio > filt[1] & pdinsulindnaratio < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               perIEQRec = {
                 filt <- str_split(filter$perIEQRec, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(percentieqrecoverypostculture > filt[1] & percentieqrecoverypostculture < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               isoIPI = {
                 filt <- str_split(filter$isoIPI, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(isletparticleindex > filt[1] & isletparticleindex < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               isoInsPerIEQ = {
                 filt <- str_split(filter$isoInsPerIEQ, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(insulinperieq > filt[1] & insulinperieq < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               pancWeight = {
                 filt <- str_split(filter$pancWeight, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(pancreasweight > filt[1] & pancreasweight < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               purityPer = {
                 filt <- str_split(filter$purityPer, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(puritypercentage > filt[1] & puritypercentage < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               IEQPerPanc = {
                 filt <- str_split(filter$IEQPerPanc, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(ieqperpancreasweight > filt[1] & ieqperpancreasweight < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               insContent = {
                 filt <- str_split(filter$insContent, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(insulincontent > filt[1] & insulincontent < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               totIEQ = {
                 filt <- str_split(filter$totIEQ, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(totalieq > filt[1] & totalieq < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               distIPI = {
                 filt <- str_split(filter$distIPI, ",")[[1]] %>% as.numeric();
                 filt.list <- donor.table[(pdisletparticleindex > filt[1] & pdisletparticleindex < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               DNAContent = {
                 filt <- str_split(filter$DNAContent, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(dnacontent > filt[1] & dnacontent < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
                 cryotubesremaining = {
                 filt <- str_split(filter$cryotubesremaining, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(as.numeric(cryotubesremaining) > filt[1] & as.numeric(cryotubesremaining) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               sftubesremaining = {
                 filt <- str_split(filter$sftubesremaining, ",")[[1]] %>% as.numeric();
                 filt.list <- iso.table[(as.numeric(sftubesremaining) > filt[1] & as.numeric(sftubesremaining) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               lgcp = {
                 filt <- str_split(filter$lgcp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avglgcp) > filt[1] & as.numeric(avglgcp) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
                
                 lgpi = {
                 filt <- str_split(filter$lgpi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avglgpi) > filt[1] & as.numeric(avglgpi) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
 
            lgcppi = {
                 filt <- str_split(filter$lgcppi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(lgcppiratio) > filt[1] & as.numeric(lgcppiratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               lgpicp = {
                 filt <- str_split(filter$lgpicp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(lgpicpratio) > filt[1] & as.numeric(lgpicpratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               
               hgcp = {
                 filt <- str_split(filter$hgcp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avghgcp) > filt[1] & as.numeric(avghgcp) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
                
                 hgpi = {
                 filt <- str_split(filter$hgpi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avghgpi) > filt[1] & as.numeric(avghgpi) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
 
            hgcppi = {
                 filt <- str_split(filter$hgcppi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(hgcppiratio) > filt[1] & as.numeric(hgcppiratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               hgpicp = {
                 filt <- str_split(filter$hgpicp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(hgpicpratio) > filt[1] & as.numeric(hgpicpratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

                 lysatecp = {
                 filt <- str_split(filter$lysatecp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avglysatecp) > filt[1] & as.numeric(avglysatecp) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
                
                 lysatepi = {
                 filt <- str_split(filter$lysatepi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(avglysatepi) > filt[1] & as.numeric(avglysatepi) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
 
            lysatecppi = {
                 filt <- str_split(filter$lysatecppi, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(lysatecppiratio) > filt[1] & as.numeric(lysatecppiratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },
               
               lysatepicp = {
                 filt <- str_split(filter$lysatepicp, ",")[[1]] %>% as.numeric();
                 filt.list <- prohormone.table[(as.numeric(lysatepicpratio) > filt[1] & as.numeric(lysatepicpratio) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               tcWeightRecovery = {
                 if(is.null(lipid.table)){ lipid.table <- as.data.table(read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))) }
                 filt <- str_split(filter$tcWeightRecovery, ",")[[1]] %>% as.numeric();
                 filt.list <- lipid.table[(as.numeric(tc_weight_recovery) > filt[1] & as.numeric(tc_weight_recovery) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               tgWeightRecovery = {
                 if(is.null(lipid.table)){ lipid.table <- as.data.table(read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))) }
                 filt <- str_split(filter$tgWeightRecovery, ",")[[1]] %>% as.numeric();
                 filt.list <- lipid.table[(as.numeric(tg_weight_recovery) > filt[1] & as.numeric(tg_weight_recovery) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               fcWeightRecovery = {
                 if(is.null(lipid.table)){ lipid.table <- as.data.table(read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))) }
                 filt <- str_split(filter$fcWeightRecovery, ",")[[1]] %>% as.numeric();
                 filt.list <- lipid.table[(as.numeric(fc_weight_recovery) > filt[1] & as.numeric(fc_weight_recovery) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               ceValue = {
                 if(is.null(lipid.table)){ lipid.table <- as.data.table(read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))) }
                 filt <- str_split(filter$ceValue, ",")[[1]] %>% as.numeric();
                 filt.list <- lipid.table[(as.numeric(ce) > filt[1] & as.numeric(ce) < filt[2]), record_id];
                 donors <- intersect(donors, filt.list);
               },

               { print('filter name not found'); })
      } else {
        break
      }
    }
  }
  
  donors <- unique(donors);
  if (length(donors) > 0) {
    donor_filename <- if (nzchar(donorListKey)) paste0("donors_", donorListKey, ".rds") else "donors.rds"
    saveRDS(donors, donor_filename);
  }

  num.donors <- as.character(length(donors));
  res <- paste0(num.donors, "RES-OK");

  return(res)
}

################################################################################

downloadDonors_fun <- function (){

    donors <- readRDS("donors.rds");
    write.table(donors, "donors.txt", quote = FALSE, row.names = FALSE, col.names = FALSE);
    return("donors.txtRES-OK")
}

################################################################################

createHITables_fun <- function(tables, filetype = 'csv', version = 'v2'){

  require(RSQLite)
  require(stringr)
  # xlsx imports rJava and STARTS A JVM -- measured at 3.58s, paid on every
  # single request because RCenter.getRConnection opens a fresh Rserve
  # connection and re-source()s this file each time. Its only use anywhere in
  # the R backend is the write.xlsx2 in the final else branch of the write
  # chain below, so load it only when that branch can actually be reached.
  #
  # The guard is deliberately the EXACT COMPLEMENT of the csv/txt/rds chain
  # rather than `filetype == "xlsx"`: HITables.java:88 lets the literal string
  # "null" through validation and RCenter.java:104 interpolates it straight
  # into the call, so any unrecognised filetype falls into the else branch and
  # still needs xlsx. isTRUE() keeps the guard total for NA / zero-length.
  # Kept here (not down in the write loop) so JVM init still happens before the
  # setwd() into download_<timestamp>, exactly as it does today.
  if(!isTRUE(filetype %in% c("csv", "txt", "rds"))) require(xlsx)
  require(data.table)   # fwrite: ~30x faster than write.csv on the omics matrices

  # set sqlite file path
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite");
  omics.db <- ifelse(version=="v2", "HI_omics_v2.sqlite", "HI_omics.sqlite")
  omics.path <- paste0(sqlite.path, omics.db);
  # RNA-seq acquisition batch per donor (record_id, batch; 371 rows + header).
  # Converted once from the source batch_rnaseq.xlsx so the download ships a plain CSV
  # like batch_info_prot.csv -- reading .xlsx here would boot rJava (~3.6s) per request.
  # NOTE: this used to point at omics_processing_input/raw/rnaseq_batch_info.txt,
  # which does NOT exist on the server -- file.copy() only warns on a missing
  # source, so the RNA-seq download silently shipped with no batch file at all.
  batch.path <- paste0(other.tables.path, "display_interface/batch_rnaseq.csv");
  # Proteomics v2 is acquired in two batches; this file maps every donor to its
  # acquisition batch (record_id, sample, batch = B1/B2). Shipped with the
  # proteomics download the same way rnaseq_batch_info.txt is (see below).
  prot.batch.path <- paste0(other.tables.path, "display_interface/batch_info_prot.csv");

  donors <- readRDS("donors.rds")
  tables <- str_split(tables, ",")[[1]]

  # split into different table types
  outcome.tables <- tables[tables %in% c('donor', 'isolation', 'gsis', 'peri_gluc', 'peri_leu', 'peri_olp', 'ephys_donor', 'seahorse','grs','cell_pro','prohormone', 'lip_extract', 'ancestry')]
  # The download dropdown for Metabolomics sends `proc_metabo` (legacy
  # singular `proc_metabolite` also accepted). v2 stores the metabolite
  # matrix split across 6 per-glucose-condition tables, all of which are
  # bundled by the proc_metabolite branch below.
  omics.tables <- tables[tables %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta','proc_contaminants', 'proc_metabolite', 'proc_metabo')]
  # Methylation lives in its own database (HI_methylation.sqlite), handled separately below.
  methylation.tables <- tables[tables %in% c('proc_methylation')]

  sqlite.tables <- list()
  
  # get outcome tables
  if(length(outcome.tables) > 0){
    con <- dbConnect(RSQLite::SQLite(), outcomes.path)
    for(i in c(1:length(outcome.tables))){
      
      if(outcome.tables[i] == "seahorse"){
        temp1 <- dbReadTable(con, "seahorse_norm_dna")
        temp1 <- temp1[temp1$record_id %in% donors, ]
        
        temp2 <- dbReadTable(con, "seahorse_norm_dna_baselineoc")
        temp2 <- temp2[temp2$record_id %in% donors, ]
        
        if(dim(temp1)[1] > 0) {
          sqlite.tables[["seahorse_norm_dna"]] <- temp1
          sqlite.tables[["seahorse_norm_dna_baselineoc"]] <- temp2
        } else {
          print("No table entries for these donors!")
          outcome.tables <- outcome.tables[-i]
        }
        
      } else if(outcome.tables[i] == "cell_pro"){
 
          temp <- dbReadTable(con, "computed")
          temp <- temp[temp$record_id %in% donors, ]
          if(dim(temp)[1] > 0) {
            sqlite.tables[["cell_pro"]] <- temp
          }
    
       } else if(outcome.tables[i] == "lip_extract"){
          temp <- read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))
          temp <- temp[temp$record_id %in% donors, ]
          if(dim(temp)[1] > 0) {
            sqlite.tables[["lip_extract"]] <- temp
          }

       } else if(outcome.tables[i] == "ancestry"){
          # Ancestry data lives in the SQLite `ancestry` table (X1_* / X2_*
          # admixture proportions, super-population labels, etc.). The
          # previous version grepped X1_/X2_ columns out of
          # numerical_donor_info.csv, but those columns were moved to the
          # ancestry table — anc.cols ended up character(0), temp dropped to
          # a 1-column data.frame, and the subsequent row-filter dropped
          # dimensions entirely. `dim(temp)[1]` then threw "argument is of
          # length zero" and crashed the whole download (frontend hung on
          # dwnld_loading forever, since Java returned null).
          if(dbExistsTable(con, "ancestry")){
            temp <- dbReadTable(con, "ancestry")
            temp <- temp[temp$record_id %in% donors, , drop = FALSE]
            if(nrow(temp) > 0) {
              sqlite.tables[["ancestry"]] <- temp
            } else {
              print("No ancestry rows for these donors!")
            }
          }

       }else {
        temp.nm <- outcome.tables[i]
        temp <- dbReadTable(con, outcome.tables[i])
        temp <- temp[temp$record_id %in% donors, ]
        
        if(dim(temp)[1] > 0) {
          sqlite.tables[[temp.nm]] <- temp
        } else {
          print("No table entries for these donors!")
          outcome.tables <- outcome.tables[-i]
        }
        
        # check if perifusion data. if so, add summary values
        if(temp.nm %in% c('peri_gluc', 'peri_leu', 'peri_olp')){
          temp <- dbReadTable(con, "function_summary")
          temp <- temp[temp$record_id %in% donors, ]
          
          if(dim(temp)[1] > 0) {
            sqlite.tables[["function_summary"]] <- temp
          }
        }
       

      }
    }
    dbDisconnect(con)
  }


  # get omics tables
  if(length(omics.tables) > 0){
    con <- dbConnect(RSQLite::SQLite(), omics.path)
    for(i in c(1:length(omics.tables))){
      table.temp <- omics.tables[i]

      if(table.temp == "proc_prot"){

        if(version == "v2"){

          # v2 proteomics is acquired in TWO batches, so the download ships all
          # four matrices instead of the single batch-1 table v1 exposed:
          #   proc_prot_b1      batch 1 donors only          (n = 134)
          #   proc_prot_b2      batch 2 donors only          (n = 100)
          #   proc_prot_combine both batches side by side, restricted to the
          #                     proteins quantified in BOTH batches, NOT
          #                     batch-corrected (batch is meant to be modelled
          #                     as a covariate)              (n = 234)
          #   proc_prot_combat  both batches, ComBat-corrected (n = 234)
          # batch_info_prot.csv (B1/B2 per donor) is added to the zip further
          # down so users can fit the covariate themselves.
          #
          # VERSION GATE: HI_omics.sqlite (v1) contains only raw_prot /
          # proc_prot / unproc_prot -- none of the four tables above exist
          # there, so a v1 request must keep the original two-table export in
          # the else branch or dbReadTable would error out and kill the whole
          # download.
          #
          # These tables carry 5 leading columns: id, gene_id, symbol, name,
          # Protein_Group. `id` is an internal DEA key (Protein_Group with
          # ';' -> '.', needed because the Java feature-id validator rejects
          # ';') and is dropped on export -- Protein_Group already holds the
          # same accession group in its original form. Remaining annotation is
          # re-ordered to put the gene symbol first.
          prot.anno <- c("symbol", "gene_id", "name", "Protein_Group")
          for(pt in c("proc_prot_b1", "proc_prot_b2", "proc_prot_combine", "proc_prot_combat")){
            if(!dbExistsTable(con, pt)) next
            temp <- dbReadTable(con, pt)
            donor.inds <- which(colnames(temp) %in% donors)
            if(length(donor.inds) > 0){
              anno.inds <- match(prot.anno, colnames(temp))
              anno.inds <- anno.inds[!is.na(anno.inds)]
              temp <- temp[, c(anno.inds, donor.inds[!donor.inds %in% anno.inds])]
              sqlite.tables[[pt]] <- temp
            }
          }

          # RAW (unprocessed) proteomics, one matrix per batch: the search-engine
          # output before any normalization. Layout DIFFERS from the proc_* tables
          # above -- 4 annotation columns (Protein_Group, Protein_Names, Genes,
          # First_Protein_Description) then donors -- so it needs its own column
          # list. A `temp[, c(1, donor.inds)]` slice would keep Protein_Group and
          # silently drop the other three.
          unproc.prot.anno <- c("Protein_Group", "Protein_Names", "Genes", "First_Protein_Description")
          for(pt in c("unproc_prot_b1", "unproc_prot_b2")){
            if(!dbExistsTable(con, pt)) next
            temp <- dbReadTable(con, pt)
            donor.inds <- which(colnames(temp) %in% donors)
            if(length(donor.inds) > 0){
              anno.inds <- match(unproc.prot.anno, colnames(temp))
              anno.inds <- anno.inds[!is.na(anno.inds)]
              temp <- temp[, c(anno.inds, donor.inds[!donor.inds %in% anno.inds])]
              sqlite.tables[[pt]] <- temp
            }
          }

        } else {

          temp <- dbReadTable(con, "unproc_prot")
          donor.inds = which(colnames(temp) %in% donors)
          if(length(donor.inds) > 0){
            temp <- temp[, c(1, donor.inds)]
            sqlite.tables[["unproc_prot"]] <- temp
          }

          temp <- dbReadTable(con, "proc_prot")
          donor.inds = which(colnames(temp) %in% donors)
          if(length(donor.inds) > 0){
            temp <- temp[, c(1, donor.inds)]
            sqlite.tables[["proc_prot"]] <- temp
          }

        }

      } else if (table.temp == "proc_nanostring_merge"){

        # v2: the per-codeset unproc_nanostring_C6555 / C8898 tables were replaced by a
        # single merged table covering all three codesets (C6555, C8898, and the R511+
        # batch) -- 259 donors vs the 190 the two per-codeset tables held.
        temp <- dbReadTable(con, "unproc_nanostring_merge")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_nanostring_merge"]] <- temp
        }

        temp <- dbReadTable(con, "proc_nanostring_merge")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["proc_nanostring_merge"]] <- temp
        }

      } else if (table.temp == "proc_rnaseq"){

        # Keep EVERY non-donor (annotation) column, not just column 1.
        # v2 proc_rnaseq leads with accession, gene_id, symbol, name, biotype, so the
        # old `temp[, c(1, donor.inds)]` slice exported the bare ENSG accession and
        # silently dropped the gene symbol, name, Entrez id and biotype -- leaving the
        # user no way to map a row to a gene. unproc_rnaseq has a single feature_id
        # column, so it is unchanged by this. Works for v1 too (3 annotation columns).
        for(rt in c("unproc_rnaseq", "proc_rnaseq")){
          if(!dbExistsTable(con, rt)) next
          temp <- dbReadTable(con, rt)
          donor.inds <- which(colnames(temp) %in% donors)
          if(length(donor.inds) > 0){
            anno.inds <- setdiff(seq_len(ncol(temp)), donor.inds)
            temp <- temp[, c(anno.inds, donor.inds)]
            sqlite.tables[[rt]] <- temp
          }
        }

      } else if (table.temp == "proc_pbrna_Alpha"){

        temp <- dbReadTable(con, "unproc_pbrna_Alpha")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_pbrna_Alpha"]] <- temp
        }

        temp <- dbReadTable(con, "proc_pbrna_Alpha")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["proc_pbrna_Alpha"]] <- temp
        }

      } else if (table.temp == "proc_pbrna_Beta"){

        temp <- dbReadTable(con, "unproc_pbrna_Beta")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_pbrna_Beta"]] <- temp
        }

        temp <- dbReadTable(con, "proc_pbrna_Beta")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["proc_pbrna_Beta"]] <- temp
        }
      }else if (table.temp == "proc_contaminants"){

        temp <- dbReadTable(con, "unproc_contaminants")
        print(head(temp))
        donor.inds = which(temp$DONOR %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[donor.inds, ]
          sqlite.tables[["unproc_contaminants"]] <- temp
        }
        
        temp <- dbReadTable(con, "proc_contaminants")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          # proc_contaminants has 4 leading info columns in the DB:
          #   1 = LOG_corrected, 2 = Tissue, 3 = Compound (analyte), 4 = Class.
          # NOTE: the DB column "LOG_corrected" is a misnomer for LOD-correction
          # (limit of detection) confirmed with the data contributor; it is NOT a
          # log transform (values are raw pg/g wet weight). Rename it to
          # LOD_corrected on export so the downloaded file is accurate.
          # The old c(1, donor.inds) kept ONLY the flag (true/false) and dropped
          # Tissue/Compound/Class, so the download had no row identifier.
          # Keep all four and reorder to the contributor README order:
          #   Tissue, Compound, Class, LOD_corrected, then the donor columns.
          temp <- temp[, c(2, 3, 4, 1, donor.inds[!donor.inds %in% 1:4])]
          colnames(temp)[colnames(temp) == "LOG_corrected"] <- "LOD_corrected"
          # Insert an explicit Units column after LOD_corrected (contributor asked
          # for the unit to be stated in the file). All concentrations are pg/g
          # (wet weight basis) -- see unproc_contaminants.UNIT and the source data.
          temp <- data.frame(temp[, 1:4, drop = FALSE],
                             Units = "pg/g (wet weight basis)",
                             temp[, 5:ncol(temp), drop = FALSE],
                             check.names = FALSE, stringsAsFactors = FALSE)
          sqlite.tables[["proc_contaminants"]] <- temp
        }

      } else if (table.temp %in% c("proc_metabolite", "proc_metabo")){

        # v2 split the metabolite matrix into 6 per-glucose-condition tables.
        # Bundle every variant that exists in the omics DB so the downloaded
        # zip carries the full set (LG / HG / ratio + ComBat-corrected
        # versions). Each variant shares the same 9 info columns
        # (compound … sub_class), so the donor-column filter logic is
        # identical across them.
        metabo.variants <- c("proc_metabolite_LG", "proc_metabolite_HG", "proc_metabolite_ratio",
                             "proc_metabolite_combat_LG", "proc_metabolite_combat_HG", "proc_metabolite_combat_ratio")
        for(mv in metabo.variants){
          if(!dbExistsTable(con, mv)) next
          mt <- dbReadTable(con, mv)
          donor.inds <- which(colnames(mt) %in% donors)
          if(length(donor.inds) > 0){
            # Keep the 9 info columns plus the matching donor columns
            mt <- mt[, c(1:9, donor.inds[!donor.inds %in% 1:9])]
            sqlite.tables[[mv]] <- mt
          }
        }

      }
    }
    dbDisconnect(con)
  }

  # get methylation table (beta values) - stored in its own database, HI_methylation.sqlite
  if(length(methylation.tables) > 0){
    meth.path <- paste0(sqlite.path, "HI_methylation.sqlite")
    if(file.exists(meth.path)){
      con <- dbConnect(RSQLite::SQLite(), meth.path)
      if(dbExistsTable(con, "proc_methylation_beta")){
        temp <- dbReadTable(con, "proc_methylation_beta")
        # first 8 columns are CpG annotation (feature_id, chr, pos_hg38, strand,
        # gene, gene_region, cgi_relation, cgi_name); the rest are donor record IDs
        donor.inds <- which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1:8, donor.inds[!donor.inds %in% 1:8])]
          sqlite.tables[["proc_methylation_beta"]] <- temp
        }
      }
      dbDisconnect(con)
    }
  }

  # create folder for downloading tables
  dir.nm <- paste0("download_", format(Sys.time(), '%Y%m%d%H%M%S'))
  dir.create(dir.nm)
  setwd(paste0("./", dir.nm))
  
  # Load the canonical column-label dictionary once. Used below to write a
  # <table>_README.txt next to each exported file so users can decode the
  # otherwise-cryptic short column names (e.g. avglgcp → "Average LG C-peptide
  # concentration (pM)"). If the dictionary file is missing the READMEs are
  # skipped silently — never block the download itself.
  var.summary <- NULL
  var.summary.path <- paste0(other.tables.path, "display_interface/proc_variable_summary_v2.csv")
  if(file.exists(var.summary.path)){
    var.summary <- tryCatch(
      read.csv(var.summary.path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  # write files into zip directory
  table.nms <- names(sqlite.tables)
  for(i in c(1:length(table.nms))){
    # write file
    # fwrite instead of write.csv / write.table. Both base writers are
    # single-threaded and dominate the download on wide omics matrices:
    # RNA-seq costs 30.3s with write.csv vs 0.6s with fwrite, the four
    # proteomics tables 22.7s vs 0.8s.
    #
    # na = "NA" IS MANDATORY, DO NOT DROP IT. fwrite defaults to na = "" while
    # write.csv/write.table default to na = "NA". Without it every missing value
    # is written as an empty field, and on CHARACTER columns that does not read
    # back as missing -- read.csv returns "" and the user's is.na() silently
    # returns FALSE. Real exposure in this export: donor.hlaa is NULL for
    # 280 of 651 donors, donor.hba1c for 91, proc_prot_b1.gene_id for 86.
    # (Numeric columns happen to survive, because an empty field parses to NA
    # for a numeric type -- which is exactly what makes the bug easy to miss.)
    #
    # NUMERIC FIDELITY: values are the same doubles, but the printed digits can
    # differ in the last place. write.csv/write.table round to 15 significant
    # digits; fwrite emits the shortest string that reads back as the exact
    # stored double -- so fwrite is the more faithful of the two, not the less.
    # Measured spread between the two: 0 for proc_prot_*, <= 1e-15 for
    # proc_metabolite_ratio (i.e. floating-point epsilon, not a data change).
    if(filetype == "csv") {
      fwrite(sqlite.tables[[i]], paste0(table.nms[i], ".csv"), quote = TRUE, na = "NA")
    } else if(filetype == "txt") {
      fwrite(sqlite.tables[[i]], paste0(table.nms[i], ".txt"), quote = FALSE, sep = "\t", na = "NA")
    } else if(filetype == "rds") {
      saveRDS(sqlite.tables[[i]], paste0(table.nms[i], ".rds"))
    } else {
      # DNA Methylation (~729k rows) cannot be written via write.xlsx2 (rJava
      # OOM/hang). Fall back to CSV for that one table so an xlsx request still
      # completes instead of stalling the whole download. The frontend also
      # blocks methylation+xlsx; this is the server-side safety net.
      if(table.nms[i] == "proc_methylation_beta"){
        write.csv(sqlite.tables[[i]], paste0(table.nms[i], ".csv"), quote = TRUE, row.names = FALSE)
      } else {
        write.xlsx2(sqlite.tables[[i]], paste0(table.nms[i], ".xlsx"), sheetName = "Details",
                    col.names = TRUE, row.names = FALSE, append = FALSE, quote = TRUE)
      }
    }

    # Write a per-table column-description README. For metadata/phenotype
    # tables (prohormone, gsis, ephys_donor, …) columns map 1:1 to entries in
    # proc_variable_summary_v2.csv. For omics tables, columns 2+ are donor
    # record IDs, so write a short structural note instead of a per-column map.
    readme.path <- paste0(table.nms[i], "_README.txt")
    cols <- colnames(sqlite.tables[[i]])
    readme.lines <- c(
      paste0("Table: ", table.nms[i]),
      paste0("Rows: ", nrow(sqlite.tables[[i]]), "    Columns: ", length(cols)),
      "",
      "Column descriptions:",
      ""
    )
    is.omics <- table.nms[i] %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq',
                                    'proc_pbrna_Alpha', 'proc_pbrna_Beta',
                                    'proc_contaminants', 'proc_metabolite',
                                    'proc_metabolite_LG', 'proc_metabolite_HG', 'proc_metabolite_ratio',
                                    'proc_metabolite_combat_LG', 'proc_metabolite_combat_HG', 'proc_metabolite_combat_ratio',
                                    'unproc_nanostring_merge',
                                    'unproc_prot', 'unproc_rnaseq',
                                    'unproc_pbrna_Alpha', 'unproc_pbrna_Beta',
                                    'unproc_contaminants')
    if(table.nms[i] == "unproc_contaminants"){
      # Column descriptions provided by the data contributor (Bruin lab).
      readme.lines <- c(readme.lines,
        "All data is reported in pg/g tissue.",
        "",
        "  - Column A includes the DONOR ID.",
        "  - Column B includes the tissue in which contaminant concentrations were measured (included adipose and pancreas).",
        "  - Column C includes the analyte measured.",
        "  - Column D includes SGS-AXYS lab flags:",
        "      - ND flag  = not detected",
        "      - NDR flag = a peak was detected but did not meet quantification criteria",
        "      - G flag   = lock mass interference present",
        "      - D flag   = data is reported from a dilution",
        "      - C flags  = coelution of compounds (e.g. C156 flag for compound for PCB 193 means that PCB 156 and PCB 193 coeluted together. This flag is reflected in the compound name)",
        "      - NA       = not applicable (ie no flag reported)",
        "  - Column E includes the sample limit of detection (LOD) for a specific donor, tissue and analyte. These values were used for LOD-adjusted data.",
        "  - Column F includes analyte concentrations that were blank corrected.",
        "  - Column G includes analyte concentrations that were blank and LOD corrected.",
        "  - Column H identifies whether analytes are \"dioxin-like\", which means they activate the aryl hydrocarbon receptor pathway.",
        "  - Column I identifies the contaminant class each analyte belongs to (i.e. dioxin/furan, PCB, or OCP).",
        "  - Column J identifies the analyte concentration units.",
        "",
        "  - 'unproc_*' = raw; 'proc_*' = pre-processed matrix used by the platform.")
    } else if(table.nms[i] == "proc_contaminants"){
      # Column descriptions provided by the data contributor (Bruin lab).
      # Order matches the export reorder above: Tissue, Compound, Class,
      # LOD_corrected, then donor IDs. The flag is LOD-correction (limit of
      # detection), NOT a log transform; values are raw pg/g (wet weight basis).
      readme.lines <- c(readme.lines,
        "  - Column A includes the tissue in which contaminant concentrations were measured (includes adipose and pancreas).",
        "  - Column B includes the analyte measured.",
        "  - Column C includes the contaminant class (i.e. dioxin/furan, PCB, or OCP).",
        "  - Column D indicates whether the data was LOD-corrected (LOD = limit of detection). Each analyte appears twice: FALSE = blank-adjusted (non-detects = 0); TRUE = blank- and LOD-adjusted (non-detects assigned half the batch-specific LOD).",
        "  - Column E gives the concentration units for all donor values: pg/g (wet weight basis).",
        "  - Columns F onward are donor record IDs, with each cell containing the analyte concentration for that donor.", 
        "  - 'unproc_*' = raw matrix; 'proc_*' = pre-processed matrix used by the platform.")
    } else if(table.nms[i] %in% c("proc_prot_b1", "proc_prot_b2", "proc_prot_combine", "proc_prot_combat")){
      # v2 proteomics carries 4 annotation columns rather than the single
      # feature column the generic omics note below describes.
      readme.lines <- c(readme.lines,
        "  - symbol: gene symbol(s) for the protein group. ';'-separated when the group maps to more than one gene (shared peptides).",
        "  - gene_id: Entrez gene ID of the protein group's leading gene (used for pathway analysis).",
        "  - name: protein name.",
        "  - Protein_Group: UniProt accession group as reported by the search engine (';'-separated when peptides are shared). Unique row identifier.",
        "  - All remaining columns are donor record IDs; each cell is the log2 normalized (Normics median) protein abundance for that donor.",
        "",
        "Proteomics was acquired in TWO batches. batch_info_prot.csv, included in this download,",
        "gives the acquisition batch (B1 / B2) for every donor. The four matrices are:",
        "",
        "  - proc_prot_b1      : batch 1 donors only.",
        "  - proc_prot_b2      : batch 2 donors only.",
        "  - proc_prot_combine : both batches side by side, restricted to proteins quantified in BOTH",
        "                        batches. NOT batch-corrected - include batch as a covariate in your model.",
        "  - proc_prot_combat  : both batches, ComBat batch-corrected. Use when batch cannot be modelled explicitly.",
        "",
        "Row counts differ between the four matrices, for two different reasons:",
        "  - b1 (8373) vs b2 (11033): each batch quantifies a different set of protein groups.",
        "  - combine (8223): the protein groups quantified in BOTH batches, i.e. the b1/b2 intersection.",
        "  - combat (7822): supplied already ComBat-corrected, and covers a subset of combine - every",
        "    protein group in combat is also in combine. The 401 present in combine but not in combat",
        "    are fully observed in both batches, so they were excluded upstream during the ComBat",
        "    processing rather than by any filter applied when building this download.")
    } else if(table.nms[i] %in% c("unproc_prot_b1", "unproc_prot_b2")){
      # Raw proteomics: search-engine output, 4 annotation columns then donors.
      readme.lines <- c(readme.lines,
        "  - Protein_Group: UniProt accession group as reported by the search engine (';'-separated when peptides are shared). Unique row identifier.",
        "  - Protein_Names: UniProt protein name(s) for the group.",
        "  - Genes: gene symbol(s) for the group (';'-separated for shared peptides).",
        "  - First_Protein_Description: description of the leading protein in the group.",
        "  - All remaining columns are donor record IDs; each cell is the RAW (un-normalized) protein intensity for that donor.",
        "",
        "These are the unprocessed per-batch matrices: unproc_prot_b1 (batch 1) and unproc_prot_b2",
        "(batch 2), before normalization, filtering or batch correction. Use the proc_prot_* matrices",
        "for analysis; use these to start from the raw intensities. Row counts differ from the",
        "corresponding proc_prot_* matrix because processing filters protein groups.")
    } else if(table.nms[i] %in% c("unproc_rnaseq", "proc_rnaseq")){
      # RNA-seq: v2 proc_rnaseq carries 5 annotation columns; unproc_rnaseq has 1.
      readme.lines <- c(readme.lines,
        "  - unproc_rnaseq: column A is feature_id (Ensembl gene accession); all remaining columns are donor record IDs holding RAW counts.",
        "  - proc_rnaseq: columns A-E are accession (Ensembl), gene_id (Entrez), symbol (gene symbol), name (gene name) and biotype; all remaining columns are donor record IDs.",
        "  - 'unproc_*' = raw count matrix; 'proc_*' = pre-processed (filtered + normalized) matrix used by the platform.",
        "",
        "RNA-seq was sequenced in several batches. batch_rnaseq.csv, included in this download,",
        "gives the acquisition batch for every donor (columns: record_id, batch).",
        "Include batch as a covariate in your model where appropriate.")
    } else if(table.nms[i] == "proc_methylation_beta"){
      readme.lines <- c(readme.lines,
        "  - feature_id: CpG probe ID (Illumina EPICv2, e.g. cg########).",
        "  - chr, pos_hg38, strand: genomic location of the CpG (hg38).",
        "  - gene, gene_region: UCSC RefGene gene symbol and gene region (blank for intergenic CpGs).",
        "  - cgi_relation, cgi_name: CpG-island relation (Island / Shore / Shelf / OpenSea) and the island name.",
        "  - All remaining columns are donor record IDs; each cell is the CpG's beta value (0-1 methylation proportion) for that donor.",
        "  - 'unproc_*' = raw/un-normalized matrix; 'proc_*' = pre-processed (filtered + normalized) matrix used by the platform.")
    } else if(is.omics){
      readme.lines <- c(readme.lines,
        "  - The first column identifies the feature (gene_id / symbol / compound / rxn / DONOR).",
        "  - All remaining columns are donor record IDs; each cell is the feature's value for that donor.",
        "  - 'unproc_*' = raw/un-normalized matrix; 'proc_*' = pre-processed (filtered + normalized) matrix used by the platform.")
    } else if(!is.null(var.summary)){
      for(cn in cols){
        if(cn == "record_id"){
          readme.lines <- c(readme.lines, "  record_id: Donor record ID")
          next
        }
        match.row <- var.summary[var.summary$column == cn, , drop = FALSE]
        if(nrow(match.row) >= 1){
          desc <- match.row$display[1]
          if(is.null(desc) || is.na(desc) || nchar(desc) == 0) desc <- match.row$axis_title[1]
          readme.lines <- c(readme.lines, paste0("  ", cn, ": ", desc))
        } else {
          readme.lines <- c(readme.lines, paste0("  ", cn, ": (no description available)"))
        }
      }
    } else {
      readme.lines <- c(readme.lines, "  (column dictionary unavailable on the server; please contact the platform team.)")
    }
    writeLines(readme.lines, readme.path)
  }

  # add batch info for RNAseq. Guarded with file.exists() like the proteomics copy
  # below -- file.copy() on a missing source only warns, so an absent file would
  # otherwise pass unnoticed (which is exactly what happened before).
  if(any(c("proc_rnaseq", "unproc_rnaseq") %in% table.nms)){
    if(file.exists(batch.path)){
      file.copy(from = batch.path, to = "batch_rnaseq.csv", overwrite = TRUE)
    }
  }

  # add batch info for proteomics (v2 only -- the v1 export has a single batch,
  # so there is nothing to disambiguate). Fires whenever any of the four v2
  # proteomics matrices made it into the zip.
  if(any(c("proc_prot_b1", "proc_prot_b2", "proc_prot_combine", "proc_prot_combat") %in% table.nms)){
    if(file.exists(prot.batch.path)){
      file.copy(from = prot.batch.path, to = "batch_info_prot.csv", overwrite = TRUE)
    }
  }

  # create zip folder
  zip.nm <- paste0(dir.nm, ".zip")
  setwd("..")
  files2zip <- dir(dir.nm, full.names = TRUE)
  print(files2zip)
  # R's zip() defaults to flags = "-r9X", i.e. MAXIMUM deflate. On the ~93 MB
  # proteomics export that alone costs ~38s. Level 1 does the same job in ~8s
  # for 2.8 MB more (41.5 -> 44.3 MB, +6.7%) -- a good trade when the file is
  # streamed straight to the browser. Only the compression digit changed; -r
  # (recurse) and -X (drop extra file attributes) are kept as before.
  zip(zipfile = zip.nm, files = files2zip, flags = "-r1X")
  
  return(paste0(zip.nm, "RES-OK"))

}


################################################################################
####generate the files for patch SC data download

getPatchSCDownload_fun <- function(filename, version="v2"){
  require(RSQLite)
  require(dplyr)
  require(data.table)
  require(stringr)
  require(jsonlite)
    
  isFilter <- FALSE
  donor.info <- readRDS(paste0(other.tables.path, "analysis_input/patchseq_donors.rds"));
  donor.table <- as.data.table(donor.info);
  donors <- donor.info$record_id
  # check if filter donor
  if(file.exists("filter.json")){
    filter <- read_json("filter.json");
    filter <- filter[names(filter) %in% c("sex", "diagnosis","age","bmi","hba1c")]
    
    filter <- lapply(filter, function(x){str_split(x, "\\|")[[1]]}) %>% as.data.frame() %>% t() %>% as.data.frame();
    filter <- filter[filter$V1 == "true", ]
    filter.vars <- rownames(filter)
    filter <- as.list(filter$V2)
    names(filter) <- filter.vars

    # for loop with switch statement inside to progressively filter donor list
    if (length(filter.vars) > 0) {
      isFilter <- TRUE
      for (i in c(1:length(filter.vars))) {
        if (length(donors) > 0) {
          filter.var <- filter.vars[i];
          switch (filter.var,
                  sex = {
                    filt <- str_split(filter$sex, ",")[[1]];
                    filt.list <- donor.table[(donorsex %in% filt), record_id];
                    donors <- intersect(donors, filt.list);
                  },
                  
                  diagnosis = {
                    filt <- str_split(filter$diagnosis, ",")[[1]];
                    if (sum(filt %in% c("Type1", "Type2", "None")) != 3) {
                      
                      filt.list <- c();
                      
                      if ("Type1" %in% filt) {
                        temp <- donor.table[diagnosis == "Type1", record_id];
                        filt.list <- unique(c(filt.list, temp));
                      }
                      
                      if ("Type2" %in% filt) {
                        temp <- donor.table[diagnosis == "Type2", record_id];
                        filt.list <- unique(c(filt.list, temp));
                      }
                      
                      if ("None" %in% filt) {
                        temp <- donor.table[diagnosis == "None", record_id];
                        filt.list <- unique(c(filt.list, temp));
                      }
                    }
                    donors <- intersect(donors, filt.list);
                  },
                  
                  age = {
                    filt <- str_split(filter$age, ",")[[1]] %>% as.numeric();
                    filt.list <- donor.table[(donorage > filt[1] & donorage < filt[2]), record_id];
                    donors <- intersect(donors, filt.list);
                  },
                  
                  hba1c = {
                    filt <- str_split(filter$hba1c, ",")[[1]] %>% as.numeric();
                    filt.list <- donor.table[(hba1c > filt[1] & hba1c < filt[2]), record_id];
                    donors <- intersect(donors, filt.list);
                  },
                  
                  bmi = {
                    filt <- str_split(filter$bmi, ",")[[1]] %>% as.numeric();
                    filt.list <- donor.table[(bodymassindex > filt[1] & bodymassindex < filt[2]), record_id];
                    donors <- intersect(donors, filt.list);
                  },
                  { print('filter name not found'); })
        } else {
          break
        }
      }
    }
    donors <- unique(donors);
  }
  if (length(donors) > 0) {
    
    dir.nm <- paste0("download_", format(Sys.time(), '%Y%m%d%H%M%S'))
    dir.create(dir.nm)
    setwd(paste0("./", dir.nm))
    if( filename =="ephys_cell.csv"){
      # ephys DB by version from API call: v1 -> HI_tables_v1.sqlite; default (v2, main project) -> HI_tables.sqlite
      sqlite.path <- ifelse(version=="v1", paste0(sqlite.path, "HI_tables_v1.sqlite"), paste0(sqlite.path, "HI_tables.sqlite"));
      con <- dbConnect(RSQLite::SQLite(), sqlite.path)
      ephys_cell <- dbReadTable(con, "ephys_cell")
      dbDisconnect(con)
      ephys_cell <- ephys_cell[ephys_cell$record_id %in% donors,]
      write.csv(ephys_cell,"ephys_cell.csv",row.names = F)
     
    }else{
      sc.h5.path <- ifelse(version=="v2", h5.v2.path, h5.path)
      table.path <- paste0(sc.h5.path, filename)
      if(isFilter){
          library(rhdf5)
        cells <- h5read(table.path, "meta/cells/cellid")
        feature_table <- h5read(table.path, "data/norm_expression") 
        cell_scores <- h5read(table.path, "/meta/cells/cellscore")
        target_cell_ids = donor.info$cell_id[donor.info$record_id %in% donors]
        matching_indices <- which(cells %in% target_cell_ids)
        feature_table <- feature_table[, matching_indices]
        cells <- cells[matching_indices]
        cell_scores <- cell_scores[matching_indices] 
        genes <- h5read(table.path, "meta/genes") %>% as.data.frame()
        
        ##creat customized h5 file
        if (file.exists(filename)) {
          file.remove(filename)
        }
        
        h5createFile(filename)
        h5createGroup(filename, "meta")
        h5createGroup(filename, "data")
        h5createDataset(filename, "/data/norm_expression", dims = dim(feature_table), storage.mode = "double")
        h5write(feature_table, filename, "data/norm_expression")   # write the actual expression values (was missing -> all-zero output)

        h5createGroup(filename, "meta/cells")
        h5write(cells, filename, "meta/cells/cellid")
        h5write(cell_scores, filename, "meta/cells/cellscore")
        
        
        h5createGroup(filename, "meta/genes")
        if(version=="v2"){
          # v2: symbol is the primary identifier
          h5write(genes[,1], filename, "meta/genes/symbol")
          if(ncol(genes) >= 2) h5write(genes[,2], filename, "meta/genes/name")
        } else {
          h5write(genes$entrez, filename, "meta/genes/entrez")
          h5write(genes$symbol, filename, "meta/genes/symbol")
          h5write(genes$name, filename, "meta/genes/name")
        }
        H5close()
      
      }else{
        file.copy(table.path, filename, overwrite = TRUE)
      }
      
    } 
  
  }
   
  zip.nm <- paste0(dir.nm, ".zip")
  setwd("..")
  files2zip <- dir(dir.nm, full.names = TRUE)
  print(files2zip)
  zip(zipfile = zip.nm, files = files2zip)
  
  return(paste0(zip.nm, "RES-OK"))
}

################################################################################


uploadDonors_fun <- function(ids, idType = 'redcap', donorListKey = ''){

    require(stringr)
    require(RSQLite)
    require(dplyr)

    donors <- str_split(ids, ",")[[1]]

    if(idType == 'rrid'){
        # set sqlite file path
        sqlite.path <- paste0(sqlite.path, "HI_tables.sqlite");
        con <- dbConnect(RSQLite::SQLite(), sqlite.path)
        donor.info <- dbReadTable(con, "donor")
        dbDisconnect(con)
        donors <- donor.info$record_id[donor.info$rrid %in% donors]
    }

    if (length(donors) > 0) {
      donor_filename <- if (nzchar(donorListKey)) paste0("donors_", donorListKey, ".rds") else "donors.rds"
      saveRDS(donors, donor_filename);
    }

    num.donors <- as.character(length(donors));
    res <- paste0(num.donors, "RES-OK");

    return(res)
}


################################################################################

getDonorInfo <- function(record_id, dataType,version="v2"){
    print(c(record_id,dataType,version))
    file.nm <- paste0(record_id, "_", dataType, ".json")
    
    if(file.exists(file.nm)) { # return name; do not recreate file
        return(paste0("RES-OK;", file.nm))

    } else { # need to create file
        require(RSQLite)
        require(dplyr)
        require(rjson)

        # set file paths
        omics.path <- ifelse(version=="v1",paste0(sqlite.path, "HI_omics.sqlite"),paste0(sqlite.path, "HI_omics_v2.sqlite"));
        outcomes.path <- ifelse(version=="v1",paste0(sqlite.path, "HI_tables_v1.sqlite"),paste0(sqlite.path, "HI_tables.sqlite"));
        print(  omics.path)
        # Convert to RedCap ID if RRID
        if(grepl("SAMN", record_id)){
          mydb <- dbConnect(SQLite(), outcomes.path)
            record_id <- dbGetQuery(mydb, "SELECT record_id FROM donor WHERE rrid = ?", params = c(record_id))[1,1]
          dbDisconnect(mydb)
          if(is.na(record_id)){ return("RES-NO") }
        } else {
          # Check if donor exists in the database
          mydb <- dbConnect(SQLite(), outcomes.path)
            check.donor <- dbGetQuery(mydb, "SELECT record_id FROM donor WHERE record_id = ?", params = c(record_id))
          dbDisconnect(mydb)
          if(dim(check.donor)[1] == 0){ return("RES-NO") }
        }

    if(dataType == "donor"){

        # get donor info
        mydb <- dbConnect(SQLite(), outcomes.path)
          if(version == "v2"){
            other.info <- dbGetQuery(mydb, "SELECT rrid, donorsex, diagnosis, diagnosis_computed, hla_a2, yearsdiabetic, donationtype FROM donor WHERE record_id = ?", params = c(record_id))
          } else {
            other.info <- dbGetQuery(mydb, "SELECT rrid, donorsex, diagnosis, diagnosis_computed, hla_a2, donationtype FROM donor WHERE record_id = ?", params = c(record_id))
            other.info$yearsdiabetic <- NA
          }
          comp.info <- dbGetQuery(mydb, "SELECT * FROM computed WHERE record_id = ?", params = c(record_id))
        dbDisconnect(mydb)

        all.info <- read.csv(paste0(other.tables.path, "display_data/numerical_donor_info.csv"))
        all.iso <- read.csv(paste0(other.tables.path, "display_data/numerical_isolation_info.csv"))
        info <- all.info[all.info$record_id == record_id,]
        iso <- all.iso[all.iso$record_id == record_id,]
     
        # create JSON structure
        donor.res <- list();
        donor.res$record_id = record_id
        donor.res$rrid = other.info$rrid
        donor.res$sex = other.info$donorsex
        donor.res$diagnosis = other.info$diagnosis
        donor.res$corStatus = other.info$diagnosis_computed
        donor.res$HLA = other.info$hla_a2
        donor.res$yearsdiabetic = ifelse(is.na(other.info$yearsdiabetic), "NA", other.info$yearsdiabetic)
        donor.res$donationtype = ifelse(is.na(other.info$donationtype) || other.info$donationtype == "", "NA", other.info$donationtype)
        donor.res$hba1c = list(value = signif(info$hba1c,2), perc = getPercentile(info$hba1c, all.info$hba1c))
        donor.res$age = list(value = info$donorage, perc = getPercentile(info$donorage, all.info$donorage))
        donor.res$bmi = list(value = round(info$bodymassindex, 1), perc = getPercentile(info$bodymassindex, all.info$bodymassindex))
        donor.res$cryotubesremaining = list(value = createIsolationItem("cryotubesremaining", "Cryopreserved Tubes", "Number remaining", iso, all.iso)[["value"]])
        donor.res$sftubesremaining = list(value = createIsolationItem("sftubesremaining", "Snap-frozen Tubes", "Number remaining", iso, all.iso)[["value"]])
        

         # write out results to json file
        donor.json <- toJSON(donor.res)
        write(donor.json, file.nm)
    }else if(dataType == "grs") {
        mydb <- dbConnect(SQLite(), outcomes.path)
        all.info <- dbReadTable(mydb, "grs")
        dbDisconnect(mydb)
        all.info[,-1] <- lapply(all.info[,-1], as.numeric)
        info <- all.info[all.info$record_id == record_id,]
         donor.res <- list();
         
           donor.res$t1d_grs = list(value = round(info$t1d_grs,2), perc = getPercentile(info$t1d_grs, all.info$t1d_grs))
           
           donor.res$t2d_grs_hard = list(value = round(info$t2d_grs_hard, 2), perc = getPercentile(info$t2d_grs_hard, all.info$t2d_grs_hard))
     
        parsc = info$t2d_beta_cell_pineg + info$t2d_beta_cell_pipos +info$t2d_body_fat+info$t2d_lipodys+info$t2d_liver_lipid+info$t2d_metabo_syndrome+info$t2d_obesity+info$t2d_glyc;

           donor.res$t2d_beta_cell_pineg = list(value = round(info$t2d_beta_cell_pineg,2), perc = getPercentile(info$t2d_beta_cell_pineg, all.info$t2d_beta_cell_pineg),par_risk= round(info$t2d_beta_cell_pineg/parsc,4))
           donor.res$t2d_beta_cell_pipos = list(value = round(info$t2d_beta_cell_pipos, 2), perc = getPercentile(info$t2d_beta_cell_pipos, all.info$t2d_beta_cell_pipos),par_risk= round(info$t2d_beta_cell_pipos/parsc,4))
           donor.res$t2d_body_fat = list(value = round(info$t2d_body_fat,2), perc = getPercentile(info$t2d_body_fat, all.info$t2d_body_fat),par_risk= round(info$t2d_body_fat/parsc,4))
           donor.res$t2d_lipodys = list(value = round(info$t2d_lipodys, 2), perc = getPercentile(info$t2d_lipodys, all.info$t2d_lipodys),par_risk= round(info$t2d_lipodys/parsc,4))
           donor.res$t2d_liver_lipid = list(value = round(info$t2d_liver_lipid,2), perc = getPercentile(info$t2d_liver_lipid, all.info$t2d_liver_lipid),par_risk= round(info$t2d_liver_lipid/parsc,4))
            donor.res$t2d_metabo_syndrome = list(value = round(info$t2d_metabo_syndrome, 2), perc = getPercentile(info$t2d_metabo_syndrome, all.info$t2d_metabo_syndrome),par_risk= round(info$t2d_metabo_syndrome/parsc,4))
           donor.res$t2d_obesity = list(value = round(info$t2d_obesity,2), perc = getPercentile(info$t2d_obesity, all.info$t2d_obesity),par_risk= round(info$t2d_obesity/parsc,4))
           donor.res$t2d_glyc = list(value = round(info$t2d_glyc,2), perc = getPercentile(info$t2d_glyc, all.info$t2d_glyc),par_risk= round(info$t2d_glyc/parsc,4))

       parsc = info$t2d_alp_neg + info$t2d_beta_cell1 +info$t2d_beta_cell2+info$t2d_bilirubin+info$t2d_cholesterol+info$t2d_hyperins+info$t2d_lipodys_1+info$t2d_lipodys_2+info$t2d_liver_lipid_soft+info$t2d_obesity_soft + info$t2d_proins+info$t2d_shbg_lpa;

           donor.res$t2d_alp_neg = list(value = round(info$t2d_alp_neg,2), perc = getPercentile(info$t2d_alp_neg, all.info$t2d_alp_neg),par_risk= round(info$t2d_alp_neg/parsc,4))
           donor.res$t2d_beta_cell1 = list(value = round(info$t2d_beta_cell1, 2), perc = getPercentile(info$t2d_beta_cell1, all.info$t2d_beta_cell1),par_risk= round(info$t2d_beta_cell1/parsc,4))
           donor.res$t2d_beta_cell2 = list(value = round(info$t2d_beta_cell2,2), perc = getPercentile(info$t2d_beta_cell2, all.info$t2d_beta_cell2),par_risk= round(info$t2d_beta_cell2/parsc,4))
           donor.res$t2d_bilirubin = list(value = round(info$t2d_bilirubin, 2), perc = getPercentile(info$t2d_bilirubin, all.info$t2d_bilirubin),par_risk= round(info$t2d_bilirubin/parsc,4))
           donor.res$t2d_cholesterol = list(value = round(info$t2d_cholesterol,2), perc = getPercentile(info$t2d_cholesterol, all.info$t2d_cholesterol),par_risk= round(info$t2d_cholesterol/parsc,4))
           donor.res$t2d_hyperins = list(value = round(info$t2d_hyperins, 2), perc = getPercentile(info$t2d_hyperins, all.info$t2d_hyperins),par_risk= round(info$t2d_hyperins/parsc,4))
           donor.res$t2d_lipodys_1 = list(value = round(info$t2d_lipodys_1,2), perc = getPercentile(info$t2d_lipodys_1, all.info$t2d_lipodys_1),par_risk= round(info$t2d_lipodys_1/parsc,4))
           donor.res$t2d_lipodys_2 = list(value = round(info$t2d_lipodys_2,2), perc = getPercentile(info$t2d_lipodys_2, all.info$t2d_lipodys_2),par_risk= round(info$t2d_lipodys_2/parsc,4))
           donor.res$t2d_liver_lipid_soft = list(value = round(info$t2d_liver_lipid_soft,2), perc = getPercentile(info$t2d_liver_lipid_soft, all.info$t2d_liver_lipid_soft),par_risk= round(info$t2d_liver_lipid_soft/parsc,4))
          donor.res$t2d_obesity_soft = list(value = round(info$t2d_obesity_soft,2), perc = getPercentile(info$t2d_obesity_soft, all.info$t2d_obesity_soft),par_risk= round(info$t2d_obesity_soft/parsc,4))
          donor.res$t2d_proins = list(value = round(info$t2d_proins,2), perc = getPercentile(info$t2d_proins, all.info$t2d_proins),par_risk= round(info$t2d_proins/parsc,2))
          donor.res$t2d_shbg_lpa = list(value = round(info$t2d_shbg_lpa,2), perc = getPercentile(info$t2d_shbg_lpa, all.info$t2d_shbg_lpa),par_risk= round(info$t2d_shbg_lpa/parsc,2))
          
     

 
     parsc =info$t1d_drdq + info$t1d_class1 + info$t1d_class2 + info$t1d_nonhla+4.71055452;
           donor.res$t1d_drdq = list(value = round(info$t1d_drdq,2), perc = getPercentile(info$t1d_drdq, all.info$t1d_drdq), par_risk= round((info$t1d_drdq+4.71)/parsc,2))
           donor.res$t1d_class1 = list(value = round(info$t1d_class1, 2), perc = getPercentile(info$t1d_class1, all.info$t1d_class1),par_risk= round(info$t1d_class1/parsc,2))
           donor.res$t1d_class2 = list(value = round(info$t1d_class2,2), perc = getPercentile(info$t1d_class2, all.info$t1d_class2),par_risk= round(info$t1d_class2/parsc,2))
           donor.res$t1d_nonhla = list(value = round(info$t1d_nonhla, 2), perc = getPercentile(info$t1d_nonhla, all.info$t1d_nonhla),par_risk= round(info$t1d_nonhla/parsc,2))
 

 
 


        donor.json <- toJSON(donor.res)
        write(donor.json, file.nm)
 
    }else if(dataType == "ancs") {

          mydb <- dbConnect(SQLite(), outcomes.path)
          all.info <- dbReadTable(mydb, "ancestry")
          dbDisconnect(mydb)
          all.info[,-1] <- lapply(all.info[,-1], as.numeric)
          group.info <- read.csv(paste0(other.tables.path, "display_interface/disc_groups.csv"))
          info <- all.info[all.info$record_id == record_id,]
          donor.res <- list();
          super = info[,grepl("X1_",names(info))]

     if(all(is.na(super[1,]))){
         donor.res$measured = FALSE
      }else{

        donor.res$measured = TRUE
        names(super) = gsub("X1_","",names(super))
        names(super) = group.info$display[match(names(super) ,group.info$group)]
        kp = apply(super,2,function(x) x>0)
        super = super[,kp,drop=F]
        pop = info[,grepl("X2_",names(info))]
        names(pop) = gsub("X2_","",names(pop))
        kp = apply(pop,2,function(x) x>0)
        pop = pop[,kp,drop=F]
        names(pop) = group.info$display[match(names(pop) ,group.info$group)]
        donor.res$super = super
        donor.res$pop =pop;
       

      }
        donor.json <- toJSON(donor.res)
        write(donor.json, file.nm)
 
    } else if(dataType == "isolation") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
       comp.info <- dbGetQuery(mydb, "SELECT * FROM computed WHERE record_id = ?", params = c(record_id))
      dbDisconnect(mydb)

      all.info <- read.csv(paste0(other.tables.path, "display_data/numerical_isolation_info.csv"))
      info <- all.info[all.info$record_id == record_id,]

      donor.res <- list();

      donor.res$record_id = record_id
      donor.res$embeddedbiopsy = info$embeddedbiopsy

      donor.res$coldischemiatime = createIsolationItem("coldischemiatime", "Cold Ischemic Time", "Hours", info, all.info)
      donor.res$puritypercentage = createIsolationItem("puritypercentage", "Purity", "Percentage", info, all.info)
      donor.res$trappedpercentage = createIsolationItem("trappedpercentage", "Trapped", "Percentage", info, all.info)
      donor.res$pancreasweight = createIsolationItem("pancreasweight", "Pancreas Weight", "Grams", info, all.info)
      donor.res$digesttime = createIsolationItem("digesttime", "Digestion Time", "Minutes", info, all.info)
      donor.res$totalieq = createIsolationItem("totalieq", "Total Islet Equivalents", "IEQ", info, all.info)
      donor.res$isletparticleindex = createIsolationItem("isletparticleindex", "Islet Particle Index", "Index (after isolation)", info, all.info)
      donor.res$ieqperpancreasweight = createIsolationItem("ieqperpancreasweight", "IEQ per Pancreas Weight", "IEQ/gram", info, all.info)
      donor.res$insulincontent = createIsolationItem("insulincontent", "Insulin Content", "mg", info, all.info)
      donor.res$insulinperieq = createIsolationItem("insulinperieq", "Insulin per IEQ", "Nanograms/IEQ after isolation", info, all.info)
      donor.res$predistributionculturetime = createIsolationItem("predistributionculturetime", "Culture Time", "Hours", info, all.info)
      donor.res$percentieqrecoverypostculture = createIsolationItem("percentieqrecoverypostculture", "IEQ Recovery", "Percentage", info, all.info)
      donor.res$pdisletparticleindex = createIsolationItem("pdisletparticleindex", "Islet Particle Index", "Index (after culture)", info, all.info)
      donor.res$pdinsulinperieq = createIsolationItem("pdinsulinperieq", "Insulin per IEQ", "Nanograms/IEQ after culture", info, all.info)
      donor.res$cryotubesremaining = createIsolationItem("cryotubesremaining", "Cryopreserved Tubes", "Number remaining", info, all.info)
      donor.res$sftubesremaining = createIsolationItem("sftubesremaining", "Snap-frozen Tubes", "Number remaining", info, all.info)

      if(dim(comp.info)[1] > 0){
        donor.res$cellType = list(exo = round(comp.info$exo_per, 2)*100, beta = round(comp.info$beta_end, 2)*100, alpha = round(comp.info$alpha_end, 2)*100, 
                                  delta = round(comp.info$delta_end, 2)*100, gamma = round(comp.info$gamma_end, 2)*100)
      } else {
        donor.res$cellType = list(exo = NA, beta = NA, alpha = NA, delta = NA, gamma = NA)
      }

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "gluc_peri") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      all.info <- dbGetQuery(mydb, "SELECT record_id, auc_baseline_3mmgluc, auc_gluc_15mmgluc, auc_gluc_6mmgluc, auc_gluc_30mmkcl FROM proc_metadata")
      dbDisconnect(mydb)

      # Convert to numeric (handles text "NA" from database)
      for (col in setdiff(names(all.info), "record_id")) {
        all.info[[col]] <- suppressWarnings(as.numeric(all.info[[col]]))
      }
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = safeSignif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$gluc1= list(value = safeSignif(info$auc_gluc_15mmgluc,2), perc = getPercentile(info$auc_gluc_15mmgluc, all.info$auc_gluc_15mmgluc))
      donor.res$gluc2 = list(value = safeSignif(info$auc_gluc_6mmgluc,2), perc = getPercentile(info$auc_gluc_6mmgluc, all.info$auc_gluc_6mmgluc))
      donor.res$kcl = list(value = safeSignif(info$auc_gluc_30mmkcl,2), perc = getPercentile(info$auc_gluc_30mmkcl, all.info$auc_gluc_30mmkcl))

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "olp_peri") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      all.info <- dbGetQuery(mydb, "SELECT record_id, auc_baseline_3mmgluc, auc_olp_1p5mmolp, auc_olp_1p5mmolp_6mmgluc, auc_olp_30mmkcl FROM proc_metadata")
      dbDisconnect(mydb)

      # Convert to numeric (handles text "NA" from database)
      for (col in setdiff(names(all.info), "record_id")) {
        all.info[[col]] <- suppressWarnings(as.numeric(all.info[[col]]))
      }
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = safeSignif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$olp1= list(value = safeSignif(info$auc_olp_1p5mmolp,2), perc = getPercentile(info$auc_olp_1p5mmolp, all.info$auc_olp_1p5mmolp))
      donor.res$olp2 = list(value = safeSignif(info$auc_olp_1p5mmolp_6mmgluc,2), perc = getPercentile(info$auc_olp_1p5mmolp_6mmgluc, all.info$auc_olp_1p5mmolp_6mmgluc))
      donor.res$kcl = list(value = safeSignif(info$auc_olp_30mmkcl,2), perc = getPercentile(info$auc_olp_30mmkcl, all.info$auc_olp_30mmkcl))


      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "leu_peri") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      all.info <- dbGetQuery(mydb, "SELECT record_id, auc_baseline_3mmgluc, auc_leu_5mmleu, auc_leu_5mmleu_6mmgluc, auc_leu_30mmkcl FROM proc_metadata")
      dbDisconnect(mydb)

      # Convert to numeric (handles text "NA" from database)
      for (col in setdiff(names(all.info), "record_id")) {
        all.info[[col]] <- suppressWarnings(as.numeric(all.info[[col]]))
      }
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = safeSignif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$leu1= list(value = safeSignif(info$auc_leu_5mmleu,2), perc = getPercentile(info$auc_leu_5mmleu, all.info$auc_leu_5mmleu))
      donor.res$leu2 = list(value = safeSignif(info$auc_leu_5mmleu_6mmgluc,2), perc = getPercentile(info$auc_leu_5mmleu_6mmgluc, all.info$auc_leu_5mmleu_6mmgluc))
      donor.res$kcl = list(value = safeSignif(info$auc_leu_30mmkcl,2), perc = getPercentile(info$auc_leu_30mmkcl, all.info$auc_leu_30mmkcl))
      
      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "seahorse") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      if(version == "v2"){
        all.info <- dbGetQuery(mydb, "SELECT record_id, calc_nonmito_oc, calc_basal_resp, calc_atp_resp, calc_proton_leak, calc_max_gluc_resp, calc_stim_gluc_resp, calc_max_resp, calc_spare_cap, calc_gluc_stim_oci, calc_gluc_stim_sparecap, calc_atp_resp_gluc, calc_bhi, calc_bhi_gluc FROM seahorse_outcome")
      } else {
        all.info <- dbGetQuery(mydb, "SELECT record_id, calc_nonmito_oc, calc_basal_resp, calc_atp_resp, calc_proton_leak, calc_max_gluc_resp, calc_stim_gluc_resp, calc_max_resp, calc_spare_cap, calc_gluc_stim_oci, calc_gluc_stim_sparecap FROM seahorse")
        all.info[,2:11] <- lapply(all.info[,2:11], as.numeric)

      # average replicates
      all.info <- na.omit(all.info)
      all.info <- aggregate(all.info[,2:11], by = list(all.info$record_id), FUN = mean)
      colnames(all.info)[1] <- "record_id"
     
    }
      dbDisconnect(mydb)

  info <- all.info[all.info$record_id == record_id, ] 
      # convert to numeric (columns may be character from t() transpose)
      

      # set units based on version
      dna_unit <- ifelse(version == "v2", "log10(pmol/(min*ug DNA))", "pmol/(min*ug DNA*bl OCR)")

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$calc_nonmito_oc = list(value = signif(info$calc_nonmito_oc, 2), perc = getPercentile(info$calc_nonmito_oc, all.info$calc_nonmito_oc), units = dna_unit)
      donor.res$calc_basal_resp = list(value = signif(info$calc_basal_resp, 2), perc = getPercentile(info$calc_basal_resp, all.info$calc_basal_resp), units = "pmol/(min*ug DNA)")
      donor.res$calc_atp_resp = list(value = signif(info$calc_atp_resp, 2), perc = getPercentile(info$calc_atp_resp, all.info$calc_atp_resp), units = dna_unit)
      donor.res$calc_proton_leak = list(value = signif(info$calc_proton_leak, 2), perc = getPercentile(info$calc_proton_leak, all.info$calc_proton_leak), units = dna_unit)
      donor.res$calc_max_gluc_resp = list(value = signif(info$calc_max_gluc_resp, 2), perc = getPercentile(info$calc_max_gluc_resp, all.info$calc_max_gluc_resp), units = dna_unit)
      donor.res$calc_stim_gluc_resp = list(value = signif(info$calc_stim_gluc_resp, 2), perc = getPercentile(info$calc_stim_gluc_resp, all.info$calc_stim_gluc_resp), units = dna_unit)
      donor.res$calc_max_resp = list(value = signif(info$calc_max_resp, 2), perc = getPercentile(info$calc_max_resp, all.info$calc_max_resp), units = dna_unit)
      donor.res$calc_spare_cap = list(value = signif(info$calc_spare_cap, 2), perc = getPercentile(info$calc_spare_cap, all.info$calc_spare_cap), units = dna_unit)
      donor.res$calc_gluc_stim_oci = list(value = signif(info$calc_gluc_stim_oci, 2), perc = getPercentile(info$calc_gluc_stim_oci, all.info$calc_gluc_stim_oci), units = "Index (unitless)")
      donor.res$calc_gluc_stim_sparecap = list(value = signif(info$calc_gluc_stim_sparecap, 2), perc = getPercentile(info$calc_gluc_stim_sparecap, all.info$calc_gluc_stim_sparecap), units = "Percent (%)")
      donor.res$calc_atp_resp_gluc = list(value = signif(info$calc_atp_resp_gluc, 2), perc = getPercentile(info$calc_atp_resp_gluc, all.info$calc_atp_resp_gluc), units = dna_unit)
      donor.res$calc_bhi = list(value = signif(info$calc_bhi, 2), perc = getPercentile(info$calc_bhi, all.info$calc_bhi), units = "Index (unitless)")
      donor.res$calc_bhi_gluc = list(value = signif(info$calc_bhi_gluc, 2), perc = getPercentile(info$calc_bhi_gluc, all.info$calc_bhi_gluc), units = "Index (unitless)")

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "gsis") {

      # get data
      all.info <- read.csv(paste0(other.tables.path, "display_data/numerical_gsis.csv"))
      info <- all.info[all.info$record_id == record_id, ]

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$culturetime2 = list(value = signif(info$culturetime2, 2), perc = getPercentile(info$culturetime2, all.info$culturetime2))
      donor.res$total_insulin_content = list(value = signif(info$total_insulin_content, 2), perc = getPercentile(info$total_insulin_content, all.info$total_insulin_content))
      donor.res$insulin_secretion_1 = list(value = signif(info$insulin_secretion_1, 2), perc = getPercentile(info$insulin_secretion_1, all.info$insulin_secretion_1))
      donor.res$insulin_secretion_2p8 = list(value = signif(info$insulin_secretion_2p8, 2), perc = getPercentile(info$insulin_secretion_2p8, all.info$insulin_secretion_2p8))
      donor.res$insulin_secretion_10 = list(value = signif(info$insulin_secretion_10, 2), perc = getPercentile(info$insulin_secretion_10, all.info$insulin_secretion_10))
      donor.res$insulin_secretion_16p7 = list(value = signif(info$insulin_secretion_16p7, 2), perc = getPercentile(info$insulin_secretion_16p7, all.info$insulin_secretion_16p7))
      donor.res$gsis_index_1_10 = list(value = signif(info$gsis_index_1_10, 2), perc = getPercentile(info$gsis_index_1_10, all.info$gsis_index_1_10))
      donor.res$gsis_index_1_16p7 = list(value = signif(info$gsis_index_1_16p7, 2), perc = getPercentile(info$gsis_index_1_16p7, all.info$gsis_index_1_16p7))
      donor.res$gsis_index_2p8_16p7 = list(value = signif(info$gsis_index_2p8_16p7, 2), perc = getPercentile(info$gsis_index_2p8_16p7, all.info$gsis_index_2p8_16p7))
      
      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "omics") {

        mydb <- dbConnect(SQLite(), omics.path)
            rnaseq <- if( record_id %in% dbListFields(mydb, "proc_rnaseq") ){ "Yes" } else { "No" }
            nano <- if( record_id %in% dbListFields(mydb, "proc_nanostring_merge") ){ "Yes" } else { "No" }
            prot <- if( record_id %in% dbListFields(mydb, "proc_prot") ){ "Yes" } else { "No" }
            if(version=="1"){
              
             met <- "No"
             contam <- "No"
            }else{
            met_cols <- dbListFields(mydb, "proc_metabolite_ratio")
            met <-  if( record_id %in% met_cols ){ "Yes" } else { "No" }
              contam <-  if( record_id %in% dbListFields(mydb, "proc_contaminants") ){ "Yes" } else { "No" }

            }
        dbDisconnect(mydb)

        mydb <- dbConnect(SQLite(), outcomes.path)
            sc <- if( record_id %in% dbReadTable(mydb, "ephys_cell")[,"record_id"]){ "Yes" } else { "No" }
        dbDisconnect(mydb)

        avail <- c(rnaseq, nano, sc, prot, met, contam)
        omics.types <- c("Bulk gene expression (RNA-seq)", "Bulk gene expression (Nanostring)", "Single-cell gene expression (RNA-seq)",
                         "Bulk protein expression")
    
       if(version=="v2"){

          omics.types <- c(    omics.types ,"Metabolomics","Environmental contaminants data")
      }

        donor.res <- list()
        for(i in c(1:length(omics.types))){
          donor.res[[i]] <- list(type = omics.types[i], avail = avail[i])
        }

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "exocytosis") {

      # get data
      all.info <- read.csv(paste0(other.tables.path, "display_data/numerical_ephys.csv"))
      info <- all.info[all.info$record_id == record_id,]

      donor.res <- list();

      donor.res$record_id = record_id

      donor.res$Beta = list(conc1 = createEphysSet("Beta", "1", info, all.info), 
                            conc5 = createEphysSet("Beta", "5", info, all.info),
                            conc10 = createEphysSet("Beta", "10", info, all.info))

      donor.res$Alpha = list(conc1 = createEphysSet("Alpha", "1", info, all.info), 
                             conc5 = createEphysSet("Alpha", "5", info, all.info),
                             conc10 = createEphysSet("Alpha", "10", info, all.info))

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    }else if(dataType == "prohormone") {
    
   all.info <-  read.csv(paste0(other.tables.path, "display_data/numerical_prohormone.csv"))  
   
   donor.res <- list();
   donor.res$record_id = record_id;
   if(!record_id %in% all.info$record_id){
     for(measure in c("lg","hg","lysate")){
       donor.res[[measure]]$cp <- list(value="NA")
       donor.res[[measure]]$pi <- list(value="NA")
       donor.res[[measure]]$cppiratio <- list(value="NA")
       donor.res[[measure]]$picpratio <- list(value="NA")
      }
  
   }else{
      info <- all.info[all.info$record_id == record_id, ]
      print(info)
      for(measure in c("lg","hg","lysate")){
        donor.res[[measure]]$cp <-  list(value=round(info[[paste0("avg",measure,"cp")]],3),perc = getPercentile(info[[paste0("avg",measure,"cp")]], all.info[[paste0("avg",measure,"cp")]][!is.na(all.info[[paste0("avg",measure,"cp")]])]))
        donor.res[[measure]]$pi <- list(value=round(info[[paste0("avg",measure,"pi")]],3),perc = getPercentile(info[[paste0("avg",measure,"pi")]], all.info[[paste0("avg",measure,"pi")]][!is.na(all.info[[paste0("avg",measure,"pi")]])]))
        donor.res[[measure]]$cppiratio <- list(value=round(info[[paste0(measure,"cppiratio")]],3),perc = getPercentile(info[[paste0(measure,"cppiratio")]], all.info[[paste0(measure,"cppiratio")]][!is.na(all.info[[paste0(measure,"cppiratio")]])]))
        donor.res[[measure]]$picpratio <- list(value=round(info[[paste0(measure,"picpratio")]],3),perc = getPercentile(info[[paste0(measure,"picpratio")]], all.info[[paste0(measure,"picpratio")]][!is.na(all.info[[paste0(measure,"picpratio")]])]))
     }
      
   }
     
   # write out results to json file
   donor.json <- toJSON(donor.res)
   write(donor.json, file.nm)
   
 }else if(dataType == "lip_extract") {

   all.info <-  read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))

   donor.res <- list();
   donor.res$record_id = record_id;
   if(!record_id %in% all.info$record_id){
       donor.res$tc <- list(value="NA")
       donor.res$tg <- list(value="NA")
       donor.res$fc <- list(value="NA")
       donor.res$ce <- list(value="NA")
   }else{
      info <- all.info[all.info$record_id == record_id, ]
      print(info)
      donor.res$tc <- list(value=round(info[["tc_weight_recovery"]],3), perc = getPercentile(info[["tc_weight_recovery"]], all.info[["tc_weight_recovery"]][!is.na(all.info[["tc_weight_recovery"]])]))
      donor.res$tg <- list(value=round(info[["tg_weight_recovery"]],3), perc = getPercentile(info[["tg_weight_recovery"]], all.info[["tg_weight_recovery"]][!is.na(all.info[["tg_weight_recovery"]])]))
      donor.res$fc <- list(value=round(info[["fc_weight_recovery"]],3), perc = getPercentile(info[["fc_weight_recovery"]], all.info[["fc_weight_recovery"]][!is.na(all.info[["fc_weight_recovery"]])]))
      donor.res$ce <- list(value=round(info[["ce"]],3), perc = getPercentile(info[["ce"]], all.info[["ce"]][!is.na(all.info[["ce"]])]))
   }

   # write out results to json file
   donor.json <- toJSON(donor.res)
   write(donor.json, file.nm)

 }
        return(paste0("RES-OK;", file.nm))
    }
}

getPercentile <- function(x, arr.x) {
  x <- suppressWarnings(as.numeric(x))
  arr.x <- suppressWarnings(as.numeric(arr.x))
  if (length(x) == 0 || all(is.na(x))) return(NA)
  arr.x <- arr.x[!is.na(arr.x)]
  if (length(arr.x) == 0) return(NA)
  cdf <- ecdf(arr.x)
  perc <- round(cdf(x), 2)
  return(perc)
}

safeSignif <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA)
  signif(x, digits)
}

createIsolationItem <- function(var, name, units, info, all.info){
  res <- list(name = name, 
              units = units, 
              value = signif(info[, var], 2), 
              perc = getPercentile(info[, var], all.info[,var]))

  return(res)
}

createEphysSet <- function(cell, glucose, info, all.info){
    info <- info[info$cell_type == cell & info$glucose_mM == glucose, ]
    all.info <- all.info[all.info$cell_type == cell & all.info$glucose_mM == glucose, ]

    ephys.set <- list()
    vars <- c("cell_size_pF_donor", "total_exocytosis_fF_pF_donor", "early_exocytosis_fF_pF_donor", "late_exocytosis_fF_pF_donor",
              "calcium_entry_pC_pF_donor", "early_ca_current_pA_pF_donor", "late_ca_current_pA_pF_donor", 
              "na_current_amp_pA_pF_donor", "na_half_inactivation_mV_donor")

    names <- c("Cell size", "Total exocytosis", "Early (RRP) exocytosis", "Late exocytosis", "Calcium entry", "Early calcium current amplitude",
               "Late calcium current amplitude", "Sodium current amplitude", "Half inactivation sodium current")

    units <- c("pF", "fF/pF", "fF/pF", "fF/pF", "pC/pF", "pA/pF", "pA/pF", "pA/pF", "mV")

    for(i in c(1:length(vars))){
        ephys.set[[i]] <- createIsolationItem(vars[i], names[i], units[i], info, all.info)
    }

    names(ephys.set) <- vars
    return(ephys.set)
}


################################################################################

getSearchResults <- function(variable_ID, variable_type, contrast, version = "v1"){

  library(RSQLite)

  # select precomputed file + column names based on version
  if(version == "v2"){
    db.file      <- "HI_precomputed_v2.sqlite"
    col.id       <- "ID"              # v2 renamed Gene_ID -> ID
    col.metaID   <- "Metadata_ID_col" # v2 renamed Metadata_ID -> Metadata_ID_col
  } else {
    db.file      <- "HI_precomputed.sqlite"
    col.id       <- "Gene_ID"
    col.metaID   <- "Metadata_ID"
  }
  sqlite.path <- paste0(sqlite.path, db.file)

  # retrieve results (column name is controlled, not user-supplied — safe to interpolate)
  mydb <- dbConnect(SQLite(), sqlite.path)
  if(variable_type == "omics"){
    # variable_ID may be a comma-separated list of underlying IDs for the same
    # gene/feature (Ensembl + Entrez + ...) — the frontend deduplicates the
    # search-name dropdown by symbol and joins all source IDs into one string
    # so that a single search returns matches across every omics that uses any
    # of those IDs (RNA-seq via Ensembl, proteomics via Entrez, etc.).
    ids <- trimws(strsplit(variable_ID, ",", fixed = TRUE)[[1]])
    ids <- ids[nzchar(ids)]
    if (length(ids) == 0) ids <- variable_ID  # safety: fall back to raw input
    placeholders <- paste(rep("?", length(ids)), collapse = ",")
    qry <- sprintf("SELECT * FROM omics_outcomes WHERE %s IN (%s)", col.id, placeholders)
    results <- dbGetQuery(mydb, qry, params = as.list(ids))
  } else {
    qry <- sprintf("SELECT * FROM omics_outcomes WHERE %s = ?", col.metaID)
    results <- dbGetQuery(mydb, qry, params = c(variable_ID))
    if(contrast != "NA"){
      results <- results[results$Contrast == contrast, ]
    }
  }
  dbDisconnect(mydb)

  # re-order by p-value & format for display
  results <- results[order(results$Adjusted_p__value), ]
  results <- results[, c("Feature", "Metadata", "Omics_type", "Coefficient", "P__value", "Adjusted_p__value",
                         col.id, "Description", "Metadata_group", "Metadata_Group_ID", col.metaID, "Omics_ID")]

  colnames(results) <- c("Feature", "Metadata", "Omics type", "Coefficient", "P_value", "Adjusted p_value",
                         "Gene ID", "Description", "Metadata group", "Metadata Group ID", "Metadata ID", "Omics ID")

  # write out results to CSV
  write.csv(results, "search_results.csv", row.names = FALSE)

  return("RES-OK")

}


################################################################################
# Knowledge Search: flexible feature <-> phenotype correlation lookup
#
# Inputs:
#   var_id       - query string (gene symbol | ENSG | entrez | metabolite name |
#                  contaminant | rxn | phenotype id)
#   var_type     - "auto" | "gene" | "ensembl" | "entrez" | "metabolite" |
#                  "contaminant" | "reaction" | "phenotype"
#   omics_filter - "all" or pipe-separated Omics_type names to restrict to
#   contrast     - for phenotype searches with categorical contrasts; "NA" for none
#   fdr_thresh   - keep Adjusted_p__value < fdr_thresh   (default 0.05)
#   rho_thresh   - keep abs(Coefficient) >= rho_thresh   (default 0.15)
#
# Behavior:
#   - Reads knowledge_search_index_v2.json to resolve var_id -> (omics_type, ID)
#     tuples. For metabolite queries, contaminants are included.
#   - Pulls matching rows from omics_outcomes, joins donor_counts for N,
#     applies thresholds, writes kg_correlations.csv, returns "RES-OK".
getKgCorrelations <- function(var_id, var_type = "auto",
                              omics_filter = "all", contrast = "NA",
                              fdr_thresh = 0.05, rho_thresh = 0.15,
                              use_raw = FALSE){

  library(RSQLite)
  library(jsonlite)

  json.path  <- paste0(other.tables.path, "display_kg/knowledge_search_index_v2.json")
  db.path    <- paste0(sqlite.path, "HI_precomputed_v2.sqlite")

  if(!file.exists(json.path)) return("RES-NO-INDEX")
  if(!file.exists(db.path))   return("RES-NO-DB")

  # cache the JSON index in the R session: first call parses (~3s),
  # subsequent calls reuse the in-memory object (~0s). The cache is
  # invalidated when the underlying file's mtime changes.
  mt <- file.info(json.path)$mtime
  if(is.null(.GlobalEnv$.kg_index_cache) ||
     !identical(.GlobalEnv$.kg_index_mtime, mt)){
    .GlobalEnv$.kg_index_cache <- fromJSON(json.path, simplifyDataFrame = FALSE)
    .GlobalEnv$.kg_index_mtime <- mt
  }
  idx <- .GlobalEnv$.kg_index_cache

  # ---- helpers ---------------------------------------------------------------
  norm <- function(x) tolower(trimws(as.character(x)))
  q    <- norm(var_id)

  # ---- resolve var_type=auto by scanning the index --------------------------
  if(var_type == "auto"){
    if(grepl("^ENSG[0-9]+", var_id, ignore.case = TRUE)){
      var_type <- "ensembl"
    } else if(grepl("^[0-9]+$", var_id)){
      var_type <- "entrez"
    } else if(grepl("^MAR[0-9]+$", var_id, ignore.case = TRUE)){
      var_type <- "reaction"
    } else if(grepl("^C[0-9]{4,5}$", var_id, ignore.case = TRUE) ||
              grepl("^HMDB[0-9]+$",  var_id, ignore.case = TRUE) ||
              grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", var_id, ignore.case = TRUE) ||
              grepl("^[A-Z]{14}$",                 var_id, ignore.case = TRUE)){
      var_type <- "metabolite"
    } else {
      # try phenotype first (short tokens often are phenotypes)
      phen_hit <- FALSE
      for(e in idx){
        if(!is.null(e$type) && e$type == "phenotype" &&
           (norm(e$id) == q || norm(e$display) == q)){
          phen_hit <- TRUE; break
        }
      }
      if(phen_hit){
        var_type <- "phenotype"
      } else {
        # then gene symbol
        sym_hit <- FALSE
        for(e in idx){
          if(!is.null(e$type) && e$type == "gene" && norm(e$symbol) == q){
            sym_hit <- TRUE; break
          }
        }
        if(sym_hit) var_type <- "gene"
        else if(any(sapply(idx, function(e) !is.null(e$type) && e$type == "reaction" && norm(e$rxn) == q))) var_type <- "reaction"
        else if(any(sapply(idx, function(e) !is.null(e$type) && e$type == "contaminant" && norm(e$compound) == q))) var_type <- "contaminant"
        else var_type <- "metabolite"
      }
    }
  }

  # ---- per-search output filename (unique, timestamped) --------------------
  safe_varID <- gsub("[^A-Za-z0-9]", "_", var_id)
  ts         <- format(Sys.time(), "%Y%m%d_%H%M%S")
  csv_name   <- paste0("kg_correlations_", var_type, "_", safe_varID, "_", ts, ".csv")

  # ---- resolve (omics_type, ID-or-Feature) query tuples ---------------------
  # each tuple: list(match_col = "ID"|"Feature", value = "...", omics_type = "...")
  tuples <- list()
  pheno_meta_id  <- NULL
  pheno_contrast <- NULL

  if(var_type == "gene"){
    for(e in idx){
      if(is.null(e$type) || e$type != "gene") next
      if(norm(e$symbol) == q){
        if(!is.null(e$omics_queries)){
          for(ot in names(e$omics_queries)){
            tuples[[length(tuples)+1]] <- list(match_col="ID", value=as.character(e$omics_queries[[ot]]), omics_type=ot)
          }
        }
        break
      }
    }
  } else if(var_type == "ensembl"){
    for(e in idx){
      if(is.null(e$type) || e$type != "gene") next
      ens <- e$ensembl
      if(!is.null(ens) && length(ens) > 0 && any(!is.na(ens) & norm(ens) == q)){
        # Ensembl maps to RNA-seq only
        if(!is.null(e$omics_queries) && !is.null(e$omics_queries[["Bulk gene expression (RNA-seq)"]])){
          tuples[[length(tuples)+1]] <- list(match_col="ID",
                                             value=as.character(e$omics_queries[["Bulk gene expression (RNA-seq)"]]),
                                             omics_type="Bulk gene expression (RNA-seq)")
        }
        break
      }
    }
  } else if(var_type == "entrez"){
    for(e in idx){
      if(is.null(e$type) || e$type != "gene") next
      ent <- e$entrez
      if(!is.null(ent) && length(ent) > 0 && any(!is.na(ent) & norm(ent) == q)){
        if(!is.null(e$omics_queries)){
          for(ot in names(e$omics_queries)){
            # Entrez does not map to Bulk gene expression (RNA-seq) (which uses Ensembl)
            if(ot == "Bulk gene expression (RNA-seq)") next
            tuples[[length(tuples)+1]] <- list(match_col="ID", value=as.character(e$omics_queries[[ot]]), omics_type=ot)
          }
        }
        break
      }
    }
  } else if(var_type == "metabolite"){
    # match against compound name, KEGG, HMDB, InChIKey (full or 14-char connectivity), GEM ID
    met_match <- function(e){
      if(is.null(e$compound)) return(FALSE)
      if(norm(e$compound) == q) return(TRUE)
      if(!is.null(e$kegg_id)  && !is.na(e$kegg_id)  && norm(e$kegg_id)  == q) return(TRUE)
      if(!is.null(e$hmdb_id)  && !is.na(e$hmdb_id)  && norm(e$hmdb_id)  == q) return(TRUE)
      if(!is.null(e$inchikey) && !is.na(e$inchikey) && norm(e$inchikey) == q) return(TRUE)
      if(!is.null(e$ik_conn)  && !is.na(e$ik_conn)  && norm(e$ik_conn)  == q) return(TRUE)
      if(!is.null(e$gem_id)   && !is.na(e$gem_id)   && norm(e$gem_id)   == q) return(TRUE)
      FALSE
    }
    for(e in idx){
      if(is.null(e$type)) next
      if(e$type == "metabolite" && met_match(e)){
        for(ot in e$omics_types){
          tuples[[length(tuples)+1]] <- list(match_col="Feature", value=as.character(e$compound), omics_type=ot)
        }
      } else if(e$type == "contaminant" && !is.null(e$compound) && norm(e$compound) == q){
        tuples[[length(tuples)+1]] <- list(match_col="Feature", value=as.character(e$compound), omics_type=as.character(e$omics_type))
      }
    }
  } else if(var_type == "contaminant"){
    for(e in idx){
      if(is.null(e$type) || e$type != "contaminant") next
      if(!is.null(e$compound) && norm(e$compound) == q){
        tuples[[length(tuples)+1]] <- list(match_col="Feature", value=as.character(e$compound), omics_type=as.character(e$omics_type))
      }
    }
  } else if(var_type == "reaction"){
    # match by rxn ID (e.g. MAR07889) or display name
    for(e in idx){
      if(is.null(e$type) || e$type != "reaction") next
      hit <- (!is.null(e$rxn)  && norm(e$rxn)  == q) ||
             (!is.null(e$name) && norm(e$name) == q)
      if(hit){
        for(ot in e$omics_types){
          tuples[[length(tuples)+1]] <- list(match_col="Feature", value=as.character(e$rxn), omics_type=ot)
        }
        break
      }
    }
  } else if(var_type == "phenotype"){
    for(e in idx){
      if(is.null(e$type) || e$type != "phenotype") next
      if(norm(e$id) == q || (!is.null(e$display) && norm(e$display) == q)){
        pheno_meta_id  <- e$id
        pheno_contrast <- e$contrast
        break
      }
    }
  } else if(var_type == "group"){
    # var_id format: "pheno1:contrast1|pheno2:contrast2|..."
    # each pair becomes a separate query; results are unioned.
    # No index lookup needed — caller already resolved the group.
  }

  # ---- optional omics_filter restriction ------------------------------------
  if(!is.null(omics_filter) && nzchar(omics_filter) && omics_filter != "all"){
    keep <- strsplit(omics_filter, "|", fixed = TRUE)[[1]]
    keep <- trimws(keep)
    if(var_type != "phenotype"){
      tuples <- Filter(function(t) t$omics_type %in% keep, tuples)
    }
  }

  # ---- query omics_outcomes -------------------------------------------------
  mydb <- dbConnect(SQLite(), db.path)
  on.exit(dbDisconnect(mydb), add = TRUE)

  if(var_type == "group"){
    # var_id format: "pheno1:contrast1|pheno2:contrast2|..."
    pairs_raw <- strsplit(var_id, "|", fixed = TRUE)[[1]]
    pairs_raw <- trimws(pairs_raw)
    pairs_raw <- pairs_raw[nzchar(pairs_raw)]
    if(length(pairs_raw) == 0){
      write.csv(data.frame(), csv_name, row.names = FALSE)
      return(paste0("RES-EMPTY:", csv_name))
    }
    parts <- list()
    qry <- "SELECT Feature, ID, Omics_type, Metadata, Contrast, Coefficient, P__value, Adjusted_p__value, Description, Metadata_group, Metadata_ID_col FROM omics_outcomes WHERE Metadata_ID_col = ? AND Contrast = ?"
    for(pr in pairs_raw){
      kv <- strsplit(pr, ":", fixed = TRUE)[[1]]
      if(length(kv) < 2) next
      ph <- trimws(kv[1]); ct <- trimws(paste(kv[-1], collapse = ":"))
      if(!nzchar(ph) || !nzchar(ct)) next
      parts[[length(parts)+1]] <- dbGetQuery(mydb, qry, params = list(ph, ct))
    }
    if(length(parts) == 0){
      write.csv(data.frame(), csv_name, row.names = FALSE)
      return(paste0("RES-EMPTY:", csv_name))
    }
    res <- do.call(rbind, parts)
    if(!is.null(omics_filter) && nzchar(omics_filter) && omics_filter != "all"){
      keep <- trimws(strsplit(omics_filter, "|", fixed = TRUE)[[1]])
      res  <- res[res$Omics_type %in% keep, , drop = FALSE]
    }
  } else if(var_type == "phenotype"){
    if(is.null(pheno_meta_id)){
      write.csv(data.frame(), csv_name, row.names = FALSE)
      return(paste0("RES-EMPTY:", csv_name))
    }
    qry <- "SELECT Feature, ID, Omics_type, Metadata, Contrast, Coefficient, P__value, Adjusted_p__value, Description, Metadata_group, Metadata_ID_col FROM omics_outcomes WHERE Metadata_ID_col = ?"
    res <- dbGetQuery(mydb, qry, params = list(pheno_meta_id))
    # resolve contrast: function arg "contrast" ("NA" means none) or the index's native contrast
    use_contrast <- contrast
    if(use_contrast == "NA" && !is.null(pheno_contrast) && !is.na(pheno_contrast) && nzchar(pheno_contrast)){
      use_contrast <- pheno_contrast
    }
    if(!is.na(use_contrast) && use_contrast != "NA" && nzchar(use_contrast)){
      res <- res[!is.na(res$Contrast) & res$Contrast == use_contrast, , drop = FALSE]
    }
    if(!is.null(omics_filter) && nzchar(omics_filter) && omics_filter != "all"){
      keep <- trimws(strsplit(omics_filter, "|", fixed = TRUE)[[1]])
      res  <- res[res$Omics_type %in% keep, , drop = FALSE]
    }
  } else {
    if(length(tuples) == 0){
      write.csv(data.frame(), csv_name, row.names = FALSE)
      return(paste0("RES-EMPTY:", csv_name))
    }
    parts <- list()
    for(t in tuples){
      col <- if(t$match_col == "Feature") "Feature" else "ID"
      qry <- sprintf("SELECT Feature, ID, Omics_type, Metadata, Contrast, Coefficient, P__value, Adjusted_p__value, Description, Metadata_group, Metadata_ID_col FROM omics_outcomes WHERE %s = ? AND Omics_type = ?", col)
      parts[[length(parts)+1]] <- dbGetQuery(mydb, qry, params = list(t$value, t$omics_type))
    }
    res <- do.call(rbind, parts)
  }

  if(is.null(res) || nrow(res) == 0){
    write.csv(data.frame(), csv_name, row.names = FALSE)
    return(paste0("RES-EMPTY:", csv_name))
  }

  # ---- apply threshold ----
  # Default: filter on Adjusted_p__value (BH-FDR). When use_raw=TRUE, fall back
  # to raw P__value < fdr_thresh — surfaced via UI button when FDR returns nothing.
  if(isTRUE(use_raw)){
    res <- res[!is.na(res$P__value) & res$P__value < fdr_thresh, , drop = FALSE]
  } else {
    res <- res[!is.na(res$Adjusted_p__value) & res$Adjusted_p__value < fdr_thresh, , drop = FALSE]
  }
  res <- res[!is.na(res$Coefficient), , drop = FALSE]
  if(nrow(res) == 0){
    write.csv(data.frame(), csv_name, row.names = FALSE)
    return(paste0("RES-EMPTY:", csv_name))
  }

  # ---- join donor_counts to get N per (Omics_type, Metadata_ID_col) ---------
  dc <- dbGetQuery(mydb, "SELECT omics_type, phenotype, N FROM donor_counts")
  res$N <- dc$N[match(paste(res$Omics_type, res$Metadata_ID_col),
                      paste(dc$omics_type, dc$phenotype))]

  # ---- format + write CSV to user folder (for user download) ---------------
  res <- res[order(res$Adjusted_p__value), , drop = FALSE]
  out <- data.frame(
    Feature         = res$Feature,
    ID              = res$ID,
    Omics_type      = res$Omics_type,
    Metadata        = res$Metadata,
    Metadata_group  = res$Metadata_group,
    Metadata_ID     = res$Metadata_ID_col,
    Contrast        = res$Contrast,
    Coefficient     = res$Coefficient,
    P_value         = res$P__value,
    Adjusted_p_value= res$Adjusted_p__value,
    N               = res$N,
    Description     = res$Description,
    stringsAsFactors = FALSE
  )
  write.csv(out, csv_name, row.names = FALSE)
  return(paste0("RES-OK:", csv_name))
}
