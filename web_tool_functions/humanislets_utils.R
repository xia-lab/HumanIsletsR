# R Functions for HumanIslets web tool
# Author: Jessica Ewald

################################################################################

getDonors_fun <- function(){
  require(RSQLite)
  require(dplyr)
  require(data.table)
  require(stringr)
  require(jsonlite)

  # set sqlite file path
  omics.path <- paste0(sqlite.path, "HI_omics.sqlite");
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite");
  
  # parse filters
  filter <- read_json("filter.json");
  
  if(filter$omics == "bulk"){
    # get donor information table
    con <- dbConnect(RSQLite::SQLite(), outcomes.path)
    donor.info <- dbReadTable(con, "donor")
    iso.info <- dbReadTable(con, "isolation")
    dbDisconnect(con)
    donor.table <- as.data.table(donor.info)
    iso.table <- as.data.table(iso.info)
  } else if(filter$omics == "singlecell"){
    donor.info <- readRDS(paste0(other.tables.path, "analysis_input/patchseq_donors.rds"));
    donor.table <- as.data.table(donor.info);
  }

  filter <- filter[!(names(filter) %in% c("usrDir", "omics"))]
  filter <- lapply(filter, function(x){str_split(x, "\\|")[[1]]}) %>% as.data.frame() %>% t() %>% as.data.frame();
  filter <- filter[filter$V1 == "true", ]

  filter.vars <- rownames(filter)
  filter <- as.list(filter$V2)
  names(filter) <- filter.vars

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
                               { print("table name not found"); });
                     } else {
                       break;
                     }
                   }
                     }
                   dbDisconnect(con);

                   # Omics tables
                   filt.omics <- filt[filt %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta')]
                   if(length(filt.omics) > 0){
                     con <- dbConnect(RSQLite::SQLite(), omics.path)
                     filt.list <- c()
                     for (j in c(1:length(filt.omics))) {
                       if (length(donors) > 1){
                         table <- dbReadTable(con, name = filt.omics[j])
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
                                 { print("table name not found"); });
                       } else {
                         break;
                       }
                     }
                     dbDisconnect(con);
                   };
     
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
               
               { print('filter name not found'); })
      } else {
        break
      }
    }
  }
  
  donors <- unique(donors);
  if (length(donors) > 0) {
    saveRDS(donors, "donors.rds");
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

createHITables_fun <- function(tables, filetype = 'csv'){

  require(RSQLite)
  require(stringr)
  require(xlsx)
  
  # set sqlite file path
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite");
  omics.path <- paste0(sqlite.path, "HI_omics.sqlite");
  batch.path <- paste0(other.tables.path, "omics_processing_input/raw/rnaseq_batch_info.txt");
  
  donors <- readRDS("donors.rds")
  tables <- str_split(tables, ",")[[1]]

  # split into different table types
  outcome.tables <- tables[tables %in% c('donor', 'isolation', 'gsis', 'peri_gluc', 'peri_leu', 'peri_olp', 'ephys_donor', 'seahorse')]
  omics.tables <- tables[tables %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta')]

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
        
      } else {
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

      } else if (table.temp == "proc_nanostring_merge"){

        temp <- dbReadTable(con, "unproc_nanostring_C6555")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_nanostring_C6555"]] <- temp
        }

        temp <- dbReadTable(con, "unproc_nanostring_C8898")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_nanostring_C8898"]] <- temp
        }

        temp <- dbReadTable(con, "proc_nanostring_merge")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["proc_nanostring_merge"]] <- temp
        }

      } else if (table.temp == "proc_rnaseq"){

        temp <- dbReadTable(con, "unproc_rnaseq")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["unproc_rnaseq"]] <- temp
        }

        temp <- dbReadTable(con, "proc_rnaseq")
        donor.inds = which(colnames(temp) %in% donors)
        if(length(donor.inds) > 0){
          temp <- temp[, c(1, donor.inds)]
          sqlite.tables[["proc_rnaseq"]] <- temp
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
      }
    }
    dbDisconnect(con)
  }
  
  # create folder for downloading tables
  dir.nm <- paste0("download_", format(Sys.time(), '%Y%m%d%H%M%S'))
  dir.create(dir.nm)
  setwd(paste0("./", dir.nm))
  
  # write files into zip directory
  table.nms <- names(sqlite.tables)
  for(i in c(1:length(table.nms))){
    # write file
    if(filetype == "csv") {
      write.csv(sqlite.tables[[i]], paste0(table.nms[i], ".csv"), quote = TRUE, row.names = FALSE)
    } else if(filetype == "txt") {
      write.table(sqlite.tables[[i]], paste0(table.nms[i], ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
    } else if(filetype == "rds") {
      saveRDS(sqlite.tables[[i]], paste0(table.nms[i], ".rds"))
    } else {
      write.xlsx2(sqlite.tables[[i]], paste0(table.nms[i], ".xlsx"), sheetName = "Details",
                  col.names = TRUE, row.names = FALSE, append = FALSE, quote = TRUE)
    }
  }

  # add batch info for RNAseq
  if("proc_rnaseq" %in% table.nms){
    file.copy(from = batch.path, to = "rnaseq_batch_info.txt", overwrite = TRUE)
  }
  
  # create zip folder
  zip.nm <- paste0(dir.nm, ".zip")
  setwd("..")
  files2zip <- dir(dir.nm, full.names = TRUE)
  print(files2zip)
  zip(zipfile = zip.nm, files = files2zip)
  
  return(paste0(zip.nm, "RES-OK"))

}


################################################################################


uploadDonors_fun <- function(ids, idType = 'redcap'){

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
      saveRDS(donors, "donors.rds");
    }
    
    num.donors <- as.character(length(donors));
    res <- paste0(num.donors, "RES-OK");

    return(res)
}


################################################################################

getDonorInfo <- function(record_id, dataType){

    file.nm <- paste0(record_id, "_", dataType, ".json")
    
    if(file.exists(file.nm)) { # return name; do not recreate file
        return(paste0("RES-OK;", file.nm))

    } else { # need to create file
        require(RSQLite)
        require(dplyr)
        require(rjson)

        # set file paths
        omics.path <- paste0(sqlite.path, "HI_omics.sqlite");
        outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite");
        
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
          other.info <- dbGetQuery(mydb, "SELECT rrid, donorsex, diagnosis, diagnosis_computed, hla_a2 FROM donor WHERE record_id = ?", params = c(record_id))
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
        donor.res$hba1c = list(value = signif(info$hba1c,2), perc = getPercentile(info$hba1c, all.info$hba1c))
        donor.res$age = list(value = info$donorage, perc = getPercentile(info$donorage, all.info$donorage))
        donor.res$bmi = list(value = round(info$bodymassindex, 1), perc = getPercentile(info$bodymassindex, all.info$bodymassindex))
        donor.res$cryotubesremaining = list(value = createIsolationItem("cryotubesremaining", "Cryopreserved Tubes", "Number remaining", iso, all.iso)[["value"]])
        donor.res$sftubesremaining = list(value = createIsolationItem("sftubesremaining", "Snap-frozen Tubes", "Number remaining", iso, all.iso)[["value"]])
        # write out results to json file
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
      
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]
      
      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = signif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$gluc1= list(value = signif(info$auc_gluc_15mmgluc,2), perc = getPercentile(info$auc_gluc_15mmgluc, all.info$auc_gluc_15mmgluc))
      donor.res$gluc2 = list(value = signif(info$auc_gluc_6mmgluc,2), perc = getPercentile(info$auc_gluc_6mmgluc, all.info$auc_gluc_6mmgluc))
      donor.res$kcl = list(value = signif(info$auc_gluc_30mmkcl,2), perc = getPercentile(info$auc_gluc_30mmkcl, all.info$auc_gluc_30mmkcl))

      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "olp_peri") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      all.info <- dbGetQuery(mydb, "SELECT record_id, auc_baseline_3mmgluc, auc_olp_1p5mmolp, auc_olp_1p5mmolp_6mmgluc, auc_olp_30mmkcl FROM proc_metadata")
      dbDisconnect(mydb)
      
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]
      
      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = signif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$olp1= list(value = signif(info$auc_olp_1p5mmolp,2), perc = getPercentile(info$auc_olp_1p5mmolp, all.info$auc_olp_1p5mmolp))
      donor.res$olp2 = list(value = signif(info$auc_olp_1p5mmolp_6mmgluc,2), perc = getPercentile(info$auc_olp_1p5mmolp_6mmgluc, all.info$auc_olp_1p5mmolp_6mmgluc))
      donor.res$kcl = list(value = signif(info$auc_olp_30mmkcl,2), perc = getPercentile(info$auc_olp_30mmkcl, all.info$auc_olp_30mmkcl))
      
      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "leu_peri") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
      all.info <- dbGetQuery(mydb, "SELECT record_id, auc_baseline_3mmgluc, auc_leu_5mmleu, auc_leu_5mmleu_6mmgluc, auc_leu_30mmkcl FROM proc_metadata")
      dbDisconnect(mydb)
      
      all.info <- na.omit(all.info)
      info <- all.info[all.info$record_id == record_id, ]
      
      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$baseline = list(value = signif(info$auc_baseline_3mmgluc,2), perc = getPercentile(info$auc_baseline_3mmgluc, all.info$auc_baseline_3mmgluc))
      donor.res$leu1= list(value = signif(info$auc_leu_5mmleu,2), perc = getPercentile(info$auc_leu_5mmleu, all.info$auc_leu_5mmleu))
      donor.res$leu2 = list(value = signif(info$auc_leu_5mmleu_6mmgluc,2), perc = getPercentile(info$auc_leu_5mmleu_6mmgluc, all.info$auc_leu_5mmleu_6mmgluc))
      donor.res$kcl = list(value = signif(info$auc_leu_30mmkcl,2), perc = getPercentile(info$auc_leu_30mmkcl, all.info$auc_leu_30mmkcl))
      
      # write out results to json file
      donor.json <- toJSON(donor.res)
      write(donor.json, file.nm)

    } else if(dataType == "seahorse") {

      # get data
      mydb <- dbConnect(SQLite(), outcomes.path)
        all.info <- dbGetQuery(mydb, "SELECT record_id, calc_nonmito_oc, calc_basal_resp, calc_atp_resp, calc_proton_leak, calc_max_gluc_resp, calc_stim_gluc_resp, calc_max_resp, calc_spare_cap, calc_gluc_stim_oci, calc_gluc_stim_sparecap FROM seahorse")
      dbDisconnect(mydb)

      # average replicates
      all.info <- na.omit(all.info)
      all.info <- aggregate(all.info[,2:11], by = list(all.info$record_id), FUN = mean)
      colnames(all.info)[1] <- "record_id"
      info <- all.info[all.info$record_id == record_id, ]

      donor.res <- list();
      donor.res$record_id = record_id
      donor.res$calc_nonmito_oc = list(value = signif(info$calc_nonmito_oc, 2), perc = getPercentile(info$calc_nonmito_oc, all.info$calc_nonmito_oc), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_basal_resp = list(value = signif(info$calc_basal_resp, 2), perc = getPercentile(info$calc_basal_resp, all.info$calc_basal_resp), units = "pmol/(min*ug DNA)")
      donor.res$calc_atp_resp = list(value = signif(info$calc_atp_resp, 2), perc = getPercentile(info$calc_atp_resp, all.info$calc_atp_resp), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_proton_leak = list(value = signif(info$calc_proton_leak, 2), perc = getPercentile(info$calc_proton_leak, all.info$calc_proton_leak), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_max_gluc_resp = list(value = signif(info$calc_max_gluc_resp, 2), perc = getPercentile(info$calc_max_gluc_resp, all.info$calc_max_gluc_resp), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_stim_gluc_resp = list(value = signif(info$calc_stim_gluc_resp, 2), perc = getPercentile(info$calc_stim_gluc_resp, all.info$calc_stim_gluc_resp), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_max_resp = list(value = signif(info$calc_max_resp, 2), perc = getPercentile(info$calc_max_resp, all.info$calc_max_resp), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_spare_cap = list(value = signif(info$calc_spare_cap, 2), perc = getPercentile(info$calc_spare_cap, all.info$calc_spare_cap), units = "pmol/(min*ug DNA*bl OCR)")
      donor.res$calc_gluc_stim_oci = list(value = signif(info$calc_gluc_stim_oci, 2), perc = getPercentile(info$calc_gluc_stim_oci, all.info$calc_gluc_stim_oci), units = "Index (unitless)")
      donor.res$calc_gluc_stim_sparecap = list(value = signif(info$calc_gluc_stim_sparecap, 2), perc = getPercentile(info$calc_gluc_stim_sparecap, all.info$calc_gluc_stim_sparecap), units = "Percent (%)")

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
            adipose <- "No"
            prot <- if( record_id %in% dbListFields(mydb, "proc_prot") ){ "Yes" } else { "No" }
            met <- "No"
            lip <- "No"
        dbDisconnect(mydb)

        mydb <- dbConnect(SQLite(), outcomes.path)
            sc <- if( record_id %in% dbReadTable(mydb, "ephys_cell")[,"record_id"]){ "Yes" } else { "No" }
        dbDisconnect(mydb)

        avail <- c(rnaseq, nano, sc, adipose, prot, met, lip)
        omics.types <- c("Bulk gene expression (RNA-seq)", "Bulk gene expression (Nanostring)", "Single-cell gene expression (RNA-seq)", 
                         "Bulk gene expression (RNA-seq: adipose tissue)", "Bulk protein expression")

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

    }

        return(paste0("RES-OK;", file.nm))
    }
}

getPercentile <- function(x, arr.x) {
  cdf <- ecdf(arr.x)
  perc <- round(cdf(x), 2)
  return(perc)
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

getSearchResults <- function(variable_ID, variable_type, contrast){

  library(RSQLite)


  # set sqlite file path
  sqlite.path <- paste0(sqlite.path, "HI_precomputed.sqlite");

  # retrieve results
  if(variable_type == "omics"){
    
    # get results
    mydb <- dbConnect(SQLite(), sqlite.path)
    results <- dbGetQuery(mydb, "SELECT * FROM omics_outcomes WHERE Gene_ID = ?", params = c(variable_ID))
    dbDisconnect(mydb)

  } else {

    # get results
    mydb <- dbConnect(SQLite(), sqlite.path)
    results <- dbGetQuery(mydb, "SELECT * FROM omics_outcomes WHERE Metadata_ID = ?", params = c(variable_ID))
    dbDisconnect(mydb)
    
    if(contrast != "NA"){
      results <- results[results$Contrast == contrast, ]
    }
  }

  # re-order by p-value & format for display
  results <- results[order(results$Adjusted_p__value), ]
  results <- results[,c("Feature", "Metadata", "Omics_type", "Coefficient", "P__value", "Adjusted_p__value", "Gene_ID", "Description", 
                        "Metadata_group", "Metadata_Group_ID", "Metadata_ID", "Omics_ID")]

  colnames(results) <- c("Feature", "Metadata", "Omics type", "Coefficient", "P_value", "Adjusted p_value", "Gene ID", "Description", 
                        "Metadata group", "Metadata Group ID", "Metadata ID", "Omics ID")

  # write out results to CSV
  write.csv(results, "search_results.csv", row.names = FALSE)
  
  return("RES-OK")
  
}

