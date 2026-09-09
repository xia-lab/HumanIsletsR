# R Functions for HumanIslets web tool
# Author: Jessica Ewald, Yao Lu

################################################################################

# visually summarize the sex, diagnosis, age, bmi, and HbA1c of donors in the downloaded data.

datasetSummary_fun <- function(tables, version="v2"){

  require(RSQLite)
  require(stringr)
  require(ggplot2)
  require(ggpubr)
  require(RColorBrewer)

  # set sqlite file path
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite")
  omics.db <- ifelse(version=="v2", "HI_omics_v2.sqlite", "HI_omics.sqlite")
  omics.path <- paste0(sqlite.path, omics.db)
  
  donors <- readRDS("donors.rds")
  tables <- str_split(tables, ",")[[1]]

  # Keep these lists in sync with createHITables_fun in humanislets_utils.R.
  # If a selectable table is missing here, both `outcome.tables` and
  # `omics.tables` can end up empty, and `for(i in c(1:length(x)))` then
  # iterates over `c(1, 0)` (R's classic empty-vector quirk) and crashes on
  # `dbReadTable(con, NA)` — this used to break the data summary plot for
  # prohormone / ephys_donor / grs / cell_pro / lip_extract / ancestry /
  # proc_contaminants / proc_metabolite.
  outcome.tables <- tables[tables %in% c('donor', 'isolation', 'gsis', 'peri_gluc', 'peri_leu', 'peri_olp', 'ephys_donor', 'seahorse', 'grs', 'cell_pro', 'prohormone', 'lip_extract', 'ancestry')]
  # The download UI sends `proc_metabo` (the dropdown value) but the v2
  # SQLite database stores the metabolite matrices as per-glucose-condition
  # tables (proc_metabolite_LG / _HG / _ratio, plus combat variants). Accept
  # both spellings here so the summary narrows donors correctly either way.
  omics.tables <- tables[tables %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta', 'proc_contaminants', 'proc_metabolite', 'proc_metabo')]

  # Narrow `donors` to those that actually have a row in EVERY selected
  # outcome table. Previously this loop read each table into `temp` but
  # discarded it, so the summary plot looked identical regardless of which
  # outcome table the user picked (prohormone vs gsis vs ephys_donor, etc.).
  # Now mirrors the omics loop's intersect step. Tables not backed by SQLite
  # (lip_extract / ancestry — CSVs) are skipped, so they neither narrow the
  # set nor crash the read.
  con <- dbConnect(RSQLite::SQLite(), outcomes.path)
  dat <- dbReadTable(con, "donor")
  if(length(outcome.tables) > 0){
    for(i in c(1:length(outcome.tables))){
      temp.nm <- outcome.tables[i]
      if(temp.nm == "cell_pro") temp.nm <- "computed"
      if(!dbExistsTable(con, temp.nm)) next
      temp <- dbReadTable(con, temp.nm)
      donors <- intersect(donors, temp$record_id)
    }
  }
  dbDisconnect(con)

  # get donors from omics tables
  if(length(omics.tables) > 0){
    con <- dbConnect(RSQLite::SQLite(), omics.path)
    for(i in c(1:length(omics.tables))){
      temp.nm <- omics.tables[i]
      # v2 stores metabolites in per-glucose-condition tables; both the
      # frontend's `proc_metabo` value and the legacy `proc_metabolite` name
      # need to resolve to a real table. For the summary plot use the ratio
      # table (proc_metabolite_ratio) as the canonical donor source — it is
      # the per-donor HG/LG ratio and so represents the same donor cohort
      # as the LG/HG split tables. Fall back through the other variants if
      # ratio is absent.
      if(temp.nm %in% c("proc_metabo", "proc_metabolite")){
        for(alt in c("proc_metabolite_ratio", "proc_metabolite_LG", "proc_metabolite_HG",
                     "proc_metabolite_combat_ratio", "proc_metabolite_combat_LG", "proc_metabolite_combat_HG")){
          if(dbExistsTable(con, alt)){ temp.nm <- alt; break }
        }
      }
      # v2 proteomics spans two acquisition batches: the donor cohort is the
      # merged proc_prot_combine matrix (234 donors), not the v1 batch-1-only
      # proc_prot (134). HI_omics.sqlite (v1) has no proc_prot_combine, so this
      # is gated on version to keep the v1 summary plot unchanged.
      if(temp.nm == "proc_prot" && version == "v2" && dbExistsTable(con, "proc_prot_combine")){
        temp.nm <- "proc_prot_combine"
      }
      if(!dbExistsTable(con, temp.nm)) next
      temp <- dbReadTable(con, temp.nm)
      # proc_contaminants stores donor IDs in a DONOR row column rather than
      # as column headers; handle both layouts.
      if("DONOR" %in% colnames(temp)){
        donors <- intersect(donors, temp$DONOR)
      } else {
        donors <- intersect(donors, colnames(temp))
      }
    }
    dbDisconnect(con)
  }
  
  # get donors in filtered list
  dat$selected <- dat$record_id %in% donors
  
  # sex in dataset
  # na.rm is required: a single NA donorsex in the selected set propagates
  # through sum() and the resulting NA Percent makes ggplot silently skip the
  # entire "In dataset" bar (this is why donor / isolation selections used to
  # render only the "All" column — the donor table has 1 NA donorsex).
  df <- data.frame(
    Sex = c("Male","Female"),
    Percent = c(sum(dat[dat$selected,c("donorsex")] == "Male", na.rm = TRUE),
                sum(dat[dat$selected,c("donorsex")] == "Female", na.rm = TRUE))
  )
  df$Percent <- round((df$Percent/length(donors))*100)
  df$Donor <- "In dataset"
  
  # sex all
  df2 <- data.frame(
    Sex = c("Male","Female"),
    Percent = c(sum(dat$donorsex == "Male", na.rm = TRUE), sum(dat$donorsex == "Female", na.rm = TRUE))
  )
  df2$Percent <- round((df2$Percent/sum(df2$Percent))*100)
  df2$Donor <- "All"
  
  df <- rbind(df, df2)
  
  
  bp.sex <- ggplot(df, aes(x=Donor, y=Percent, fill=Sex)) +
    geom_bar(width=0.9, stat="identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle("Sex") + 
    theme(legend.position = c(0.5, 0.8), legend.title = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) 
  
  # diabetes type in dataset
  # na.rm = TRUE on every sum so a single NA T1/T2_diabetes value cannot
  # turn the In-dataset Percent into NA and hide the bar (same failure mode
  # as the Sex section above).
  df <- data.frame(
    Type = c("Type 1", "Type 2", "None"),
    Percent = c(sum(dat[dat$selected,c("T1_diabetes")] == "1", na.rm = TRUE),
                sum(dat[dat$selected,c("T2_diabetes")] == "1", na.rm = TRUE),
                sum(dat[dat$selected,c("T1_diabetes")] == "0" & dat[dat$selected,c("T2_diabetes")] == "0", na.rm = TRUE))
  )
  df$Percent <- round((df$Percent/length(donors))*100)
  df$Donor <- "In dataset"
  
  # diabetes type all
  df2 <- data.frame(
    Type = c("Type 1", "Type 2", "None"),
    Percent = c(sum(dat$T1_diabetes == "1", na.rm = TRUE), sum(dat$T2_diabetes == "1", na.rm = TRUE), 
                sum(dat$T1_diabetes == "0" & dat$T2_diabetes == "0", na.rm = TRUE))
  )
  df2$Percent <- round((df2$Percent/sum(df2$Percent))*100)
  df2$Donor <- "All"
  df <- rbind(df, df2)
  
  bp.type <- ggplot(df, aes(x=Donor, y=Percent, fill=Type)) +
    geom_bar(width=0.9, stat="identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle("Diabetes Type") +
    theme(legend.position = c(0.5, 0.8), legend.title = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black"))
  
  # hist
  df <- data.frame(
    Age = dat$donorage,
    BMI = dat$bodymassindex,
    HbA1c = dat$hba1c,
    Donor = "All"
  )
  df <- rbind(df, data.frame(
    Age = dat$donorage[dat$selected],
    BMI = dat$bodymassindex[dat$selected],
    HbA1c = dat$hba1c[dat$selected],
    Donor = "In dataset"
  ))
  
  hist.age <- ggplot(df, aes(x = Age, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.2, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  hist.bmi <- ggplot(df, aes(x = BMI, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  hist.hba1c <- ggplot(df, aes(x = HbA1c, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  # arrange plots
  col1 <- ggarrange(bp.sex, bp.type, ncol = 1)
  col2 <- ggarrange(hist.age, hist.bmi, hist.hba1c, ncol = 1)
  fig <- ggarrange(col1, col2, nrow = 1, widths = c(1,2))
  
  # plot data 
  Cairo::Cairo(file = "data_summary.png", unit="in", res=300, width=5.5, 
               height= 9.5, type="PNG", bg="white");
  print(fig)
  dev.off()
  
 ggsave("data_summary.svg", plot = fig, width = 5.5, height = 9.5, dpi = 300)
  ggsave("data_summary.pdf", plot = fig, width = 5.5, height = 9.5, dpi = 300)

  # string to return
  donor.num <- length(donors)
  
  per.missing <- NA
  
  res <- paste(c(donor.num, per.missing, "RES-OK"), collapse = ";")
  return(res)

}


################################################################################

metaplot_fun = function(meta, x.label, imgNm) {
  print(c(meta,x.label))
  require(RColorBrewer)
  require(RSQLite)
  require(dplyr)
  require(ggplot2)
  require(Cairo)
  
  # set file paths
  sqlite.path <- paste0(sqlite.path, "HI_tables.sqlite");
  
  donor.meta <- c("donorage", "donorsex", "donationtype", "bodymassindex", "hba1c", "hla_a2", 
                  "diagnosis", "other_condition", "percentieqrecoverypostculture", "pdisletparticleindex",
                  "pdinsulinperieq", "pdinsulindnaratio", "predistributionculturetime");

   grs.meta <- c("t1d_grs","t1d_drdq","t1d_class1","t1d_class2","t1d_nonhla",
                   "t2d_grs","t2d_beta_cell","t2d_proins","t2d_obesity","t2d_lipodys", "t2d_liver_lipid")

 prohormone.meta <- c("avglgcp","avglgpi","lgcppiratio" ,"lgpicpratio","avghgcp","avghgpi","hgcppiratio","avghgcp","avglysatepi","lysatecppiratio", "lysatepicpratio")

  lipid.meta <- c("tc_weight_recovery", "tg_weight_recovery", "fc_weight_recovery", "ce")

  if (meta %in% lipid.meta) {
    # lipid data is in lip_extract.csv (smaller/faster than numerical_donor_info.csv)
    all.info <- read.csv(paste0(other.tables.path, "display_data/lip_extract.csv"))
    dat <- data.frame(meta = all.info[[meta]], stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat$meta), , drop = FALSE]
  } else {
    if (meta %in% donor.meta) {
      table.name <- "donor"
    }else if (meta %in%  grs.meta){
      table.name <- "grs"
    }else if (meta %in%  prohormone.meta){
      table.name <- "prohormone"
    } else {
      table.name <- "isolation"
    }

    query <- paste0('SELECT ', meta, ' FROM ', table.name)

    # get data
    mydb <- dbConnect(RSQLite::SQLite(), sqlite.path)
    dat <- dbGetQuery(mydb, query)
    dbDisconnect(mydb)
  }
  if(colnames(dat) %in% c("cryotubesremaining","sftubesremaining")){
      dat =data.frame(meta=as.numeric(dat[!is.na(dat[,1]),]),stringsAsFactors = F) 
   }else{
     colnames(dat) <- "meta"
   }
  x.label <- gsub("\\.", " ", x.label)
  
  # create ggplot plot object
  if( class(dat$meta) == "numeric" ) {
    a <- ggplot(dat, aes(x = meta)) +
      geom_histogram(fill="black") +
      theme_bw() +
      xlab(x.label) +
      ylab("Count") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    a <- ggplot(dat, aes(x = meta)) +
      geom_bar(fill="black") +
      theme_bw() +
      xlab(x.label) +
      ylab("Count") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # plot data 
  imgNm <- paste0(imgNm, ".png");
  Cairo::Cairo(file = imgNm, unit="in", res=300, width=6, 
               height= 5, type="PNG", bg="white");
  print(a)
  dev.off()
  
  # add status code
  res <- paste0(imgNm, "RES-OK");
  return(res);

}

################################################################################

## Region-level methylation profile data for the DMR click-plot. Given a region
## ("chr:start-end") and the (continuous or discrete) phenotype, pulls the beta value of
## every CpG in the region for every donor + each donor's phenotype value, and writes
## region_profile.csv (long: record_id, cg, pos, beta, meta). The frontend draws one
## beta-vs-position line per donor, coloured by the phenotype. Read-only.
.runMethylationRegionProfile <- function(region.id, meta.var){
  suppressMessages({ library(RSQLite); library(DBI) })
  rid   <- gsub("[^0-9A-Za-z:_-]", "", region.id)
  parts <- strsplit(rid, "[:-]")[[1]]                 # accepts chr:start-end OR chr-start-end
  if(length(parts) < 3) return("RES-NO; bad region id")
  chr <- parts[1]
  s   <- suppressWarnings(as.integer(parts[2]))
  e   <- suppressWarnings(as.integer(parts[3]))
  if(is.na(s) || is.na(e) || !nzchar(chr)) return("RES-NO; bad region id")

  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_methylation.sqlite"))
  on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)
  b <- dbGetQuery(con, sprintf(
        "SELECT * FROM proc_methylation_beta WHERE chr = '%s' AND pos_hg38 >= %d AND pos_hg38 <= %d",
        chr, s, e))
  if(nrow(b) == 0) return("RES-NO; no CpGs found in region")

  donorcol <- colnames(b)[-(1:8)]                                   # cols 1:8 are annotation
  beta_vec <- suppressWarnings(as.numeric(as.matrix(b[, donorcol, drop = FALSE])))  # column-major

  con2 <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  meta <- dbGetQuery(con2, "SELECT * FROM proc_metadata"); try(dbDisconnect(con2), silent = TRUE)
  mv <- meta[[meta.var]][match(donorcol, meta$record_id)]
  mv[mv == "NA" | mv == ""] <- NA
  mv <- suppressWarnings(as.numeric(mv))

  n <- nrow(b)
  long <- data.frame(
    record_id = rep(donorcol, each = n),
    cg        = rep(b$feature_id, times = length(donorcol)),
    pos       = rep(suppressWarnings(as.integer(b$pos_hg38)), times = length(donorcol)),
    beta      = beta_vec,
    meta      = rep(mv, each = n),
    stringsAsFactors = FALSE
  )
  long <- long[!is.na(long$beta) & !is.na(long$meta), ]
  if(nrow(long) == 0) return("RES-NO; no usable data for this region/phenotype")
  data.table::fwrite(long, "region_profile.csv")
  return(paste0("RES-OK;", n, ";", length(unique(long$record_id))))
}

plotOmicsFeature <- function(gene.id, meta.var, omics.type, donors = "all", cell, glucose, batch = NULL, peakmode = NULL, metabo.gluc = NULL, version, mode="tool"){

  library(RSQLite)
  library(dplyr)
  library(ggplot2)
  library(see)
  library(Cairo)
  library(rhdf5)

  # Region-level (DMR) profile plot: gene.id = "chr:start-end". Writes region_profile.csv
  # (per-CpG beta for every donor + each donor's phenotype) and returns early.
  if(identical(omics.type, "proc_methylation_region")){
    return(.runMethylationRegionProfile(gene.id, meta.var))
  }

  # Escape single quotes for SQL — contaminant compound IDs (e.g. "4,4'-DDE") contain them
  gene.id <- gsub("'", "''", gene.id, fixed = TRUE)
print(c(gene.id, meta.var, omics.type,version))
  # set omics category
  if(omics.type == "proc_scrna"){
    sc.h5.path <- ifelse(version=="v2", h5.v2.path, h5.path)
    table.path <- paste0(sc.h5.path, "sc_", cell, "_", glucose, ".h5")
    meta.table <- "ephys_cell"
  } else if(omics.type == "proc_pbrna"){
    omics.category <- "bulk"
    omics.type <- paste0(omics.type, "_", cell)
    meta.table <- "proc_metadata"
  } else if(omics.type == "proc_flux"){
    # Handle flux simulation with glucose treatment
    omics.category <- "bulk"
    if(!is.null(metabo.gluc) && metabo.gluc != ""){
      if(metabo.gluc == "HG_LG_ratio"){
        omics.type <- paste0(omics.type, "_ratio")
      } else {
        omics.type <- paste0(omics.type, "_", metabo.gluc)
      }
    }
    meta.table <- "proc_metadata"
  } else if(omics.type == "proc_metabolite"){
    omics.category <- "bulk"
    # DB stores pre-split tables: proc_metabolite_{HG|LG|ratio}
    #                         and proc_metabolite_combat_{HG|LG|ratio}
    # Resolve to the exact table name here so the SELECT query is correct.
    base_nm <- if(!is.null(batch) && batch == "combat") "proc_metabolite_combat" else "proc_metabolite"
    gluc_sfx <- switch(if(!is.null(metabo.gluc) && metabo.gluc != "") metabo.gluc else "LG",
                       "HG"          = "HG",
                       "LG"          = "LG",
                       "HG_LG_ratio" = "ratio",
                       "LG")
    omics.type <- paste0(base_nm, "_", gluc_sfx)
    meta.table <- "proc_metadata"
  } else if(omics.type == "proc_methylation"){
    # Boxplot shows beta values (interpretable 0-1 methylation fraction);
    # resolve to the beta table (lives in HI_methylation.sqlite, opened below).
    omics.category <- "bulk"
    omics.type <- "proc_methylation_beta"
    meta.table <- "proc_metadata"
  } else if(grepl("^proc_contaminants", omics.type)){
    # proc_contaminants_{pancreas|adipose} -> single raw table proc_contaminants,
    # filtered by Tissue + LOG_corrected (the association used LOG_corrected = 'true').
    omics.category <- "bulk"
    cont.tissue <- sub("^proc_contaminants_?", "", omics.type)               # "pancreas"/"adipose"
    cont.tissue <- paste0(toupper(substr(cont.tissue, 1, 1)), substr(cont.tissue, 2, nchar(cont.tissue)))  # -> "Pancreas"/"Adipose"
    omics.type <- "proc_contaminants"
    meta.table <- "proc_metadata"
  } else {
    omics.category <- "bulk"
    if(grepl("_donor", meta.var)){
        meta.table <- "ephys_donor"
    } else {
        meta.table <- "proc_metadata"
    }
  }
 
  if(mode == "tool"){
    n <-3
    # get data
    if(omics.type == "proc_scrna"){
      # v2: gene IDs are symbols directly; v1: gene IDs are entrez
      if(version=="v2"){
        genes <- h5read(table.path, "meta/genes/symbol")
      } else {
        genes <- h5read(table.path, "meta/genes/entrez")
      }
      gene.ind <- which(genes == gene.id)

      feat.info <- h5read(table.path, "meta/genes") %>% as.data.frame()
      feat.info <- feat.info[gene.ind, ]
      if(version=="v2"){
        feat.info$symbol <- genes[gene.ind]
      }
      cells <- h5read(table.path, "meta/cells/cellid")

      feat.dat <- h5read(table.path, "data/norm_expression", index = list(gene.ind, NULL)) %>% t() %>% as.data.frame()
      rownames(feat.dat) <- cells

    } else {
      if(version=="v2"){
        # Methylation lives in its own DB; everything else is in HI_omics_v2.
        if(grepl("proc_methylation", omics.type)){
          mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_methylation.sqlite"))
        } else {
          mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
        }
        if( omics.type == "proc_rnaseq"){
        colnm <- "accession"   # v2: omics_outcomes.ID = ENSG accession (raw gene_id column is entrez)
        n <- 5
        }else if(omics.type == "proc_contaminants"){
              colnm <- "Compound"   # match by contaminant name; ID = name
              n <- 4                # LOG_corrected, Tissue, Compound, Class
        }else if(grepl("proc_methylation", omics.type)){
              colnm <- "feature_id"
              n <- 8
        }else if(grepl("proc_flux", omics.type)){
              colnm <- "rxn"
               n<-4
        }else if(grepl("proc_metabolite", omics.type)){
              colnm <- "compound"   # v2: match by compound NAME (InChIKey inconsistent for some stereoisomers); frontend passes Feature (name) as gene.id
              n <- 9
        }else if(grepl("proc_metabo", omics.type)){
              colnm <- "peak"
              n <- 6
        }else if(grepl("proc_pbrna", omics.type)){
              colnm <- "symbol"
              n <- 4
        }else if(omics.type == "proc_nanostring_merge"){
              colnm <- "symbol"   # v2: omics_outcomes.ID = symbol (raw gene_id column is entrez)
              n <- 4              # symbol, gene_id, ensembl, name
        }else if(omics.type %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")){
              colnm <- "symbol"   # v2: omics_outcomes.ID = symbol (raw id column is Protein_Group)
              n <- 5          # id, gene_id, symbol, name, Protein_Group
        }else{
               colnm <- "gene_id"
               n <- 3
       }


     }else{
      mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
      colnm <- "gene_id"
      n <- 3
     }

        # proc_prot_v2 is a label; the real table is the combined or combat one.
        qtab <- if(omics.type == "proc_prot_v2"){ if(!is.null(batch) && batch == "combat") "proc_prot_combat" else "proc_prot_combine" } else omics.type
        if(omics.type == "proc_contaminants"){
          # single table holds both tissues x {LOG_corrected true/false}; pick the exact row
          query <- paste0("SELECT * FROM ", qtab, " WHERE ", colnm, " = '", gene.id,
                          "' AND Tissue = '", cont.tissue, "' AND LOG_corrected = 'true'")
        } else {
          query <- paste0(  "SELECT * FROM ", qtab,  " WHERE ", colnm, " = '", gene.id, "'")
        }

      print(query)

        feat.dat <- dbGetQuery(mydb, query)

      # v2 ONLY -- v1 keeps its single gene_id lookup above, untouched.
      # Callers disagree on which identifier they send (the precomputed Feature page sends
      # omics_outcomes.ID, Omics View sends the dea_results Feature). If the declared key
      # misses, retry against the OTHER identifier columns of the SAME omics table, so the
      # symbol <-> Ensembl mapping is resolved inside the table itself -- no mapping file.
      # `gene_id` (Entrez) is deliberately EXCLUDED: Entrez is many-to-one against Ensembl
      # (e.g. 51463 -> GPR89A + 2x GPR89B), so it cannot identify a single feature.
      # Lookup only: whichever column matches, the row returned is the same row, so the
      # plotted values are unchanged.
      if(version == "v2" && nrow(feat.dat) == 0){
        alt.cols <- setdiff(intersect(c("symbol", "ensembl", "accession", "id", "Protein_Group",
                                        "compound", "Compound", "rxn", "peak", "feature_id"),
                                      dbListFields(mydb, qtab)),
                            c(colnm, "gene_id"))
        for(alt in alt.cols){
          if(omics.type == "proc_contaminants"){
            q2 <- paste0("SELECT * FROM ", qtab, " WHERE ", alt, " = '", gene.id,
                         "' AND Tissue = '", cont.tissue, "' AND LOG_corrected = 'true'")
          } else {
            q2 <- paste0("SELECT * FROM ", qtab, " WHERE ", alt, " = '", gene.id, "'")
          }
          print(q2)
          feat.dat <- tryCatch(dbGetQuery(mydb, q2), error = function(e) feat.dat[0, , drop = FALSE])
          if(nrow(feat.dat) > 0) break
        }
      }
      dbDisconnect(mydb)

      # v2 ONLY. Guard the two shapes that used to crash further down with the opaque
      # "'names' attribute [3] must be the same length as the vector [2]": no match at all
      # (t() of a 0-row frame -> a 0-column matrix -> the merge loses the feature column),
      # and a multi-row match (-> several feature columns).
      if(version == "v2"){
        if(nrow(feat.dat) == 0){
          return(paste0("RES-NO; feature '", gene.id, "' not found in ", qtab,
                        " (expects an Ensembl ID or a gene symbol)"))
        }
        if(nrow(feat.dat) > 1) feat.dat <- feat.dat[1, , drop = FALSE]
      }

      # For legacy proc_metabo peak tables only, filter by peakmode before processing
      if(grepl("^proc_metabo(_|$)", omics.type) && !grepl("proc_metabolite", omics.type) &&
         !is.null(peakmode) && peakmode != "" && peakmode != "mix"){
        feat.dat <- feat.dat[feat.dat$mode == peakmode,, drop = FALSE]
      }

       feat.info <- feat.dat[,c(1:n)]
      feat.dat <- feat.dat[,-c(1:n)]

      # NOTE: proc_metabolite data is now stored in pre-split tables (proc_metabolite_HG,
      # proc_metabolite_LG, proc_metabolite_ratio, and combat variants).  The table name
      # was already resolved above, so column names are clean record_ids — no splitting needed.
      # The block below is kept only as a safety guard and will never execute for the current DB.

      feat.dat <- t(feat.dat)
    }
 
    # get metadata
    mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    query <- paste0("SELECT * FROM ", meta.table)
    metadata <- dbGetQuery(mydb, query)

    # filter by cell and glucose if relevant
    if(grepl("_donor", meta.var)){
        metadata <- metadata[metadata$cell_type == cell, ]
    } else if(omics.type == "proc_scrna"){
        metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    }

    query <- paste0("SELECT * FROM proc_variable_summary WHERE column='", meta.var, "'")
    meta.info <- dbGetQuery(mydb, query)
    # Handle variables not in proc_variable_summary (e.g., lipid extraction)
    lipid.vars <- c('tc_weight_recovery', 'tg_weight_recovery', 'fc_weight_recovery', 'ce')
    if(nrow(meta.info) == 0 && meta.var %in% lipid.vars){
      meta.info <- data.frame(column = meta.var, type = "cont", stringsAsFactors = FALSE)
    }
    dbDisconnect(mydb)
     if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }

    saveRDS(feat.info,"savedAnalysis/feat.info.rds")
    saveRDS(feat.dat,"savedAnalysis/feat.dat.rds")
    saveRDS(meta.info,"savedAnalysis/meta.info.rds")
    saveRDS(metadata,"savedAnalysis/metadata.rds")
    rcmd <- gsub("tool","local",rcmd)
    write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);

    if(!file.exists("savedAnalysis/plotOmicsFeature.R")){
       dump("plotOmicsFeature", file = "savedAnalysis/plotOmicsFeature.R",append=T)
    }  
     
  } else {
      feat.info <- readRDS("feat.info.rds")
      feat.dat <- readRDS("feat.dat.rds")
      meta.info <- readRDS("meta.info.rds")
      metadata <- readRDS("metadata.rds")
  }
 
  meta.type <- meta.info$type
  if(omics.type != "proc_scrna"){

   if(grepl("proc_prot", omics.type)){
     ylab = "Log2(normalized protein abundance)"

   }else if(omics.type=="proc_contaminants"){
       ylab = "Normalized contaminant abundance"

   }else if(grepl("proc_methylation", omics.type)){
       ylab = "Methylation (beta value)"

   }else if(grepl("proc_flux", omics.type)){
       ylab = "Normalized flux rate"

  }else if(grepl("proc_metabolite", omics.type)){
      ylab = "Normalized metabolite abundance"

  }else if(grepl("proc_metabo", omics.type)){
      ylab = "Normalized peak abundance"

  }else{
      ylab = "Normalized transcript abundance"

  }

    # Set plot title based on omics type
    if(grepl("proc_metabolite", omics.type)){
      plot.title <- feat.info$compound
    } else if(grepl("proc_metabo", omics.type)){
      plot.title <- ifelse(!is.null(feat.info$Compound) && feat.info$Compound != "" && feat.info$Compound != "NA",
                           feat.info$Compound, feat.info$peak)
    } else if(grepl("proc_methylation", omics.type)){
      plot.title <- if(!is.null(feat.info$gene) && !is.na(feat.info$gene) && feat.info$gene != "")
                      paste0(gene.id, " (", feat.info$gene, ")") else gene.id
    } else if(omics.type == "proc_contaminants"){
      plot.title <- feat.info$Compound
    } else {
      plot.title <- feat.info$symbol
    }

    df <- merge(metadata[,c("record_id", meta.var)], feat.dat, by.x = "record_id", by.y = "row.names", all = FALSE)
   df <- na.omit(df)
    colnames(df) <- c("record_id", "meta", "feature")

    # filter by donors
    if(donors == "subset"){
      # donors.rds is written by the donor filter (humanislets_utils.R ~L466) and ONLY when the
      # filter matched at least one donor. Selecting the "Subset" radio sets donors="subset"
      # immediately, so cancelling the dialog -- or a filter that matched nobody -- leaves this
      # file absent and readRDS died with an unreadable "cannot open compressed file" error.
      # Report it instead. Deliberately NOT falling back to all donors: that would silently
      # analyse a different donor set than the one requested.
      if(!file.exists("donors.rds")){
        return("RES-NO; no donor subset is saved - choose donors in the Subset dialog, or set Donors back to All")
      }
      donor.list <- readRDS("donors.rds")
      df <- df[df$record_id %in% donor.list,]
    }
    df$feature <- as.numeric(df$feature)
    if(meta.type == "disc"){
       p <- ggplot2::ggplot(df, aes(x = meta, y = feature, fill = meta)) +
       # width caps the violin at 0.6 category units (ggplot default is 0.9): with
       # scale="area" a tight distribution such as a large "None" group was pushed out
       # to the full width and dominated the panel. Relative widths are unchanged, so
       # the violin still reflects the density peak.
       geom_violin(trim = FALSE, width = 0.6, aes(color = meta), show.legend = FALSE) +
       geom_jitter(height = 0, width = 0.05, show.legend = FALSE) +
       theme(legend.position = "none") +  xlab(meta.info$axis_title) +
       ylab(ylab) +
        ggtitle(plot.title)+
       stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE) +
       scale_fill_okabeito() + 
       scale_color_okabeito() + 
       theme(axis.text.x = element_text(angle=90, hjust=1),
          plot.title = element_text(hjust = 0.5)) +
       theme_bw() +
       theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
    } else {
      p <- ggplot2::ggplot(df, aes(x=meta, y=feature))+
      geom_point(size=2) + theme_bw()  + geom_smooth(method=lm,se=T)+
      xlab(meta.info$axis_title) + ylab(ylab) +
        ggtitle(plot.title)+
      theme(axis.text.x = element_text(angle=90, hjust=1),  plot.title = element_text(hjust = 0.5)) + 
      guides(size="none") + theme_bw() +
      theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
    }

  } else {
    df <- merge(metadata[,c("cell_id","record_id", meta.var)], feat.dat, by.x = "cell_id", by.y = "row.names", all = FALSE)
    df <- na.omit(df)
    colnames(df) <- c("cell_id", "record_id","meta", "feature")
    df$meta <- as.numeric(df$meta)
    df$feature <- as.numeric(df$feature)

    # filter by donors  (single-cell path; same missing-file case as the bulk branch above)
    if(donors == "subset"){
      if(!file.exists("donors.rds")){
        return("RES-NO; no donor subset is saved - choose donors in the Subset dialog, or set Donors back to All")
      }
      donor.list <- readRDS("donors.rds")
      cells.keep <- metadata$cell_id[metadata$record_id %in% donor.list]
      df <- df[df$cell_id %in% cells.keep,]
    }
    
    p <- ggplot2::ggplot(df, aes(x=meta, y=feature))+
        geom_point(size=2, alpha = 0.25) + theme_bw() + geom_smooth(method=lm, se=T) +
        xlab(meta.info$display) + ylab("Normalized transcript abundance") +
        ggtitle(feat.info$symbol)+
        theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
         guides(size="none") + theme_bw() +
        theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
  }
  
  write.csv(df,"df.csv",row.names=F)
  # Build ONE sanitized filename stem and reuse it for every output and the
  # returned imgNm so the frontend's download URL always matches the file on
  # disk. Previous version used three different gsub patterns across .png /
  # .svg / .pdf / .csv (one even contained the literal string "_|-"), which
  # silently broke downloads for any feature with both "," and "-" — i.e.
  # most environmental contaminants (e.g. "4,4'-DDE", "1,2,3,4,6,7,8-HPCDD").
  fnm.stem <- gsub("[^A-Za-z0-9]+", "_", gene.id)   # sanitize ALL non-alphanumerics (metabolite names have / : ( ) [ ] + ? space) so file writes never break
  imgNm <- paste0(fnm.stem, ".png")
  Cairo(file = imgNm, width=4, height=4, type="png", bg="white", unit="in", dpi=150)
    print(p)
  dev.off()
  ggplot2::ggsave(paste0(fnm.stem, ".svg"), plot = p, width = 4, height = 4, dpi = 300)
  ggplot2::ggsave(paste0(fnm.stem, ".pdf"), plot = p, width = 4, height = 4, dpi = 300)
  write.csv(df, paste0(fnm.stem, ".csv"), row.names = F)
  return(paste0("RES-OK", imgNm))
}


################################################################################

plotExpressionByCell <- function(gene.id, display = FALSE){
  library(rhdf5)
  library(ggplot2)
  library(ggbeeswarm)
  library(dplyr)
  library(data.table)
  library(Cairo)

  # set hdf5 file path (ma_singlecell.h5 is unchanged between v1/v2, stays in hdf5/)
  hdf5.path <- paste0(h5.path, "ma_singlecell.h5");

  symbols <- h5read(hdf5.path, "meta/genes/symbol")
  gene.ind <- which(symbols == gene.id)
  if(length(gene.ind) == 0){
    entrez <- h5read(hdf5.path, "meta/genes/entrez")
    gene.ind <- which(entrez == gene.id)
    if(length(gene.ind) == 0){return("RES-NO")}
  }

  symbol <- symbols[gene.ind]
  
  df <- data.frame(celltypes = h5read(hdf5.path, "meta/cells/celltype"),
                   value = c(h5read(hdf5.path, "data/norm_expression", index = list(gene.ind, 1:143181))))
  H5close()
  
  dt <- as.data.table(df)
  summary <- dt[,.(mean = mean(value), num_obs = .N, num_zero = sum(value==0),
                   pro_zero = sum(value==0)/.N, pro_zero = sum(value==0)/.N), by = celltypes]
  
  df <- df[df$value > 0, ]
  
  df <- merge(df, summary, by = "celltypes")
  df$pro_zero <- round(df$pro_zero*100,2)
  
  # plotting takes a long time if so many cells
  max.num <- 1500
  cell.freq <- table(df$celltypes) %>% as.data.frame()
  abun.cells <- cell.freq$Var1[cell.freq$Freq > max.num] %>% as.character()
  num.cells <- cell.freq$Freq[cell.freq$Freq > max.num] %>% as.numeric()
  if(length(abun.cells) > 0){
    for(i in c(1:length(abun.cells))){
      inds <- which(df$celltypes == abun.cells[i])
      rem.inds <- sample(inds, num.cells[i]-max.num, replace = FALSE)
      df <- df[-rem.inds, ]
    }
  }
  
  df$celltypes <- factor(df$celltypes, 
                         levels = c("beta", "alpha", "delta", "PP", "epsilon",
                                    "acinar", "ductal", "endothelial",
                                    "Schwann",
                                    "activated stellate", "quiescent stellate", "mast", "macrophage", "Cytotoxic T"))
  p <- ggplot(df, aes(x = celltypes, y = value, fill = pro_zero, color = pro_zero)) +
    geom_quasirandom(varwidth = TRUE, alpha = 0.5) +
    theme_bw() +
    coord_flip() +
    ggtitle(paste0("Non-zero expression of ", symbol, " in non-diabetic donors"))+
    ylab("Normalized expression values") +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y=element_text(size=12, color = "black"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA)) +
    labs(fill = "% Zero Counts") +
    guides(color = "none") +
    scale_x_discrete(drop=FALSE)
 
  if(display){
    return(p)
  } else {
    # print out plot
    imgNm <- paste0(gene.id, "_singlecell.png")
    Cairo(file = imgNm, width=7, height=9, type="png", bg="white", unit="in", dpi=150)
    print(p)
    dev.off()
    ggplot2::ggsave(paste0(gene.id, "_singlecell.svg"), plot = p, width = 7, height = 9, dpi = 300)  
    ggplot2::ggsave(paste0(gene.id, "_singlecell.pdf"), plot = p, width = 7, height = 9, dpi = 300)  
    write.csv(df,paste0(gene.id, "_singlecell.csv"),row.names=F)
    return(paste0("RES-OK",imgNm))
  }
}


################################################################################

## Methylation pathway heatmap (self-contained; no HDF5Array/minfi). The gene-level
## path below can't serve methylation (features are CpGs in a separate DB, not gene
## rows in HI_omics). Here we map the pathway's Entrez genes -> CpGs (precomputed
## map), load their beta from HI_methylation.sqlite, aggregate to per-gene mean beta
## (gene x donor), then reuse the same phenotype-ordered, smoothed heatmap. Genes
## that contain a significant CpG are marked '***'.
.plotMethylationPathwayHeatmap <- function(pathName, funcLib, analysisVar, donors = "all", mode = "tool"){
  suppressMessages({ library(RSQLite); library(pheatmap); library(smoother); library(data.table) })
  lib.path <- paste0(other.tables.path, "libraries/")

  # pathway -> Entrez genes
  libraryRDS <- readRDS(paste0(lib.path, funcLib, ".rds"))
  wi <- which(libraryRDS$term == pathName)
  if(length(wi) == 0) return("RES-NO; pathway not found in library")
  pathGenes <- as.character(libraryRDS$sets[[wi[1]]])

  # Entrez genes -> CpGs (precomputed map)
  map <- qs::qread(paste0(lib.path, "methylation_cpg_anno.qs"))
  map <- map[as.character(map$entrez) %in% pathGenes & !is.na(map$symbol) & map$symbol != "", c("cpg", "symbol")]
  if(nrow(map) == 0) return("RES-NO; no CpGs map to this pathway")

  # load beta for those CpGs from the methylation DB
  con  <- dbConnect(SQLite(), paste0(sqlite.path, "HI_methylation.sqlite"))
  qin  <- paste(sprintf("'%s'", unique(map$cpg)), collapse = ",")
  beta <- dbGetQuery(con, paste0("SELECT * FROM proc_methylation_beta WHERE feature_id IN (", qin, ")"))
  dbDisconnect(con)
  if(nrow(beta) == 0) return("RES-NO; no beta values for pathway CpGs")
  anno.cols  <- c("feature_id","chr","pos_hg38","strand","gene","gene_region","cgi_relation","cgi_name")
  donor.cols <- setdiff(colnames(beta), anno.cols)
  rownames(beta) <- beta$feature_id
  bmat <- as.matrix(beta[, donor.cols, drop = FALSE])

  # aggregate to per-gene mean beta (across the gene's CpGs), per donor
  g    <- map[map$cpg %in% rownames(bmat), ]
  bm   <- bmat[g$cpg, , drop = FALSE]
  sums <- rowsum(bm, group = g$symbol, na.rm = TRUE)
  cnts <- rowsum((!is.na(bm)) + 0, group = g$symbol)
  feature_table <- as.data.frame(sums / cnts)
  feature_table[is.na(feature_table)] <- NA   # 0/0 -> NaN -> NA

  # phenotype (donor metadata)
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  metadata  <- dbGetQuery(con, "SELECT * FROM proc_metadata")
  meta.info <- dbGetQuery(con, paste0("SELECT * FROM proc_variable_summary WHERE column='", analysisVar, "'"))
  dbDisconnect(con)
  if(!(analysisVar %in% colnames(metadata))) return("RES-NO; phenotype not found for methylation")
  meta.type <- if(nrow(meta.info)) meta.info$type else "cont"
  if(donors == "subset" && file.exists("donors.rds")){
    metadata <- metadata[metadata$record_id %in% readRDS("donors.rds"), ]
  }
  metadata <- metadata[, c("record_id", analysisVar)]
  if(meta.type == "cont") metadata[, analysisVar] <- as.numeric(metadata[, analysisVar])
  metadata <- metadata[metadata$record_id %in% colnames(feature_table), ]
  metadata <- metadata[!is.na(metadata[, analysisVar]), ]
  metadata <- metadata[order(metadata[, analysisVar]), ]
  if(nrow(metadata) == 0) return("RES-NO; no donors with this phenotype")
  feature_table <- feature_table[, match(metadata$record_id, colnames(feature_table)), drop = FALSE]

  # mark genes that contain a significant CpG ('***')
  if(file.exists("dea_results.csv")){
    dea <- tryCatch(data.table::fread("dea_results.csv", select = c("Feature", "sig")), error = function(e) NULL)
    if(!is.null(dea)){
      sig.syms <- unique(map$symbol[map$cpg %in% dea$Feature[dea$sig != "NS"]])
      di <- which(rownames(feature_table) %in% sig.syms)
      if(length(di)) rownames(feature_table)[di] <- paste0(rownames(feature_table)[di], "***")
    }
  }

  # drop genes with too many missing donors
  feature_table <- feature_table[apply(feature_table, 1, function(x) sum(!is.na(x))) > 7, , drop = FALSE]
  if(nrow(feature_table) == 0) return("RES-NO; too few non-missing values")

  # smooth across phenotype-ordered donors + heatmap (same look as plotPathwayHeatmap)
  colAnn <- data.frame(Metadata = metadata[, 2]); rownames(colAnn) <- metadata$record_id
  sm <- feature_table
  sm[is.na(feature_table)] <- min(sm, na.rm = TRUE)
  sm <- as.data.frame(t(apply(sm, 1, function(x)
          smoother::smth.gaussian(as.numeric(x), window = 0.03, alpha = 2.5, tails = TRUE, na.rm = TRUE))))
  colnames(sm) <- colnames(feature_table); rownames(sm) <- rownames(feature_table)
  sm[is.na(feature_table)] <- NA

  hm <- pheatmap::pheatmap(sm, show_colnames = ncol(feature_table) < 50, annotation_col = colAnn,
          cluster_cols = FALSE, legend = TRUE, annotation_names_col = FALSE, scale = "row",
          silent = TRUE, border_color = NA, main = pathName, na_col = "black",
          color = colorRampPalette(c("#000080", "#FFF300"))(30))
  hm.height <- nrow(feature_table) * 0.2 + 0.75
  Cairo::Cairo(file = "pathway_heatmap.png", unit = "in", res = 300, width = 10, height = hm.height, type = "PNG", bg = "white")
  print(hm); dev.off()
  ggplot2::ggsave("pathway_heatmap.svg", plot = hm, width = 10, height = hm.height, dpi = 300)
  ggplot2::ggsave("pathway_heatmap.pdf", plot = hm, width = 10, height = hm.height, dpi = 300)
  write.csv(sm, "pathway_heatmap.csv")
  return("RES-OK")
}

plotPathwayHeatmap <- function(pathName, funcLib, analysisVar, omicsType, varGroup, cell, glucose, donors = "all", version="v2", mode="local", batch = "NA", metaboGluc = "NA"){

    # Methylation: features are CpGs in a separate DB (HI_methylation.sqlite), not
    # gene rows in HI_omics -> handle in a self-contained, HDF5Array-free helper.
    if(identical(omicsType, "proc_methylation")){
      return(.plotMethylationPathwayHeatmap(pathName, funcLib, analysisVar, donors, mode))
    }
   
    library(dplyr)
    library(pheatmap)
    library(RSQLite)
    library(rhdf5)
    library(smoother)
    library(RColorBrewer)
    
    # set omics category
    if(omicsType == "proc_scrna"){
      sc.h5.path <- ifelse(version=="v2", h5.v2.path, h5.path)
      table.path <- paste0(sc.h5.path, "sc_", cell, "_", glucose, ".h5")
      meta.table <- "ephys_cell"
    } else if(omicsType == "proc_pbrna"){
      omics.category <- "bulk"
      omics.type <- paste0(omicsType, "_", cell)
      meta.table <- "proc_metadata"
    } else {
      omics.category <- "bulk"
      if(grepl("_donor", analysisVar)){
        meta.table <- "ephys_donor"
      } else {
        meta.table <- "proc_metadata"
      }
    }
    if(mode=="tool"){
     # set library file path
    lib.path <- paste0(other.tables.path, "libraries/");
    
    # get omics data
    if(omicsType == "proc_scrna"){

      feature_info <- h5read(table.path, "meta/genes") %>% as.data.frame()
      H5close()

      if(version=="v2"){
        # v2 H5 has named columns: ensembl, entrez, name, symbol
        feature_info$gene_id <- feature_info$entrez
      } else {
        feature_info <- feature_info[,c(1,3,2)]
        colnames(feature_info) <- c("gene_id", "symbol", "name")
      }

      cells <- h5read(table.path, "meta/cells/cellid")
      H5close()

      feature_table <- h5read(table.path, "data/norm_expression") %>% as.data.frame()
      H5close()

      colnames(feature_table) <- cells

      genes.keep <- !is.na(feature_info$gene_id)
      feature_table <- feature_table[genes.keep, ]
      feature_info <- feature_info[genes.keep, ]
      # Use symbol as rowname for display; gene_id (entrez) is used for pathway matching
      rownames(feature_table) <- if(version=="v2") feature_info$symbol else feature_info$gene_id
      feature_table[feature_table < 0.0001] <- NA
    } else {
      omics.db <- ifelse(version=="v2", "HI_omics_v2.sqlite", "HI_omics.sqlite")
      mydb <- dbConnect(SQLite(), paste0(sqlite.path, omics.db))
      if(omicsType == "proc_pbrna"){
        table.nm <- paste0(omicsType, "_", cell)
        feature_table <- dbReadTable(mydb, table.nm)
      } else if(omicsType == "proc_metabolite"){
        # Metabolite data is stored in pre-split tables: proc_metabolite_{HG|LG|ratio}
        # and proc_metabolite_combat_{HG|LG|ratio}
        base_nm  <- if(!is.null(batch) && batch == "combat") "proc_metabolite_combat" else "proc_metabolite"
        gluc_sfx <- switch(if(!is.null(metaboGluc) && metaboGluc != "" && metaboGluc != "NA") metaboGluc else "LG",
                           "HG" = "HG", "LG" = "LG", "HG_LG_ratio" = "ratio", "LG")
        table.nm <- paste0(base_nm, "_", gluc_sfx)
        feature_table <- dbReadTable(mydb, table.nm)
      } else {
        # proc_prot_v2 is a label -> the combined or combat table.
        qtab <- if(omicsType == "proc_prot_v2"){ if(!is.null(batch) && batch == "combat") "proc_prot_combat" else "proc_prot_combine" } else omicsType
        feature_table <- dbReadTable(mydb, qtab)
      }
      dbDisconnect(mydb)
      if(version=="v2" && omicsType == "proc_pbrna"){
        # v2 pbrna: 4 info columns (symbol, gene_id, ensembl, name)
        # symbol is the primary key; deduplicate by p-value (keep smallest P_value per symbol)
        feature_info  <- feature_table[, c(1:4)]
        feature_table <- feature_table[, -c(1:4)]
        dup_syms <- unique(feature_info$symbol[duplicated(feature_info$symbol)])
        if(length(dup_syms) > 0){
          if(file.exists("dea_results.csv")){
            dea.tmp <- read.csv("dea_results.csv")
            for(sym in dup_syms){
              dup_idx <- which(feature_info$symbol == sym)
              pvals   <- dea.tmp$P_value[match(feature_info$symbol[dup_idx], dea.tmp$Feature)]
              keep    <- dup_idx[which.min(replace(pvals, is.na(pvals), Inf))]
              remove  <- dup_idx[dup_idx != keep]
              feature_table <- feature_table[-remove, ]
              feature_info  <- feature_info[-remove, ]
            }
          } else {
            # fallback: keep first occurrence when dea_results not available
            keep_idx      <- !duplicated(feature_info$symbol)
            feature_info  <- feature_info[keep_idx, ]
            feature_table <- feature_table[keep_idx, ]
          }
        }
        rownames(feature_table) <- as.character(feature_info$symbol)
      } else if(omicsType == "proc_metabolite"){
        # 9 info columns: compound, inchikey, ik_conn, hmdb_id, kegg_id, gem_id,
        #                  super_class, main_class, sub_class
        rownames(feature_table) <- feature_table$compound
        feature_info <- feature_table[,c(1:9)]
        feature_table <- feature_table[,-c(1:9)]
      } else if(version == "v2" && omicsType == "proc_rnaseq"){
        # v2 rnaseq: 5 info columns (accession, gene_id, symbol, name, biotype)
        # accession (Ensembl) is unique; gene_id (Entrez) is NOT — multiple Ensembl IDs can share one Entrez ID
        # Use accession as internal rowname; match pathGenes (Entrez IDs) via gene_id column later
        feature_info  <- feature_table[, c(1:5)]
        feature_table <- feature_table[, -c(1:5)]
        rownames(feature_table) <- feature_info$accession
      } else if(omicsType %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")){
        # proteomics v2: 5 info cols (id, gene_id, symbol, name, Protein_Group); rowname = id.
        # Match pathGenes (Entrez) via the gene_id column (is_v2_gene_expr below); display = symbol.
        prot.cols     <- c("id","gene_id","symbol","name","Protein_Group")
        feature_info  <- feature_table[, prot.cols]
        feature_table <- feature_table[, !(colnames(feature_table) %in% prot.cols)]
        rownames(feature_table) <- feature_info$id
      } else if(version == "v2" && omicsType == "proc_nanostring_merge"){
        # v2 nanostring: 4 info columns (symbol, gene_id, ensembl, name).
        # gene_id (Entrez) stays the rowname so pathGenes match directly -- nanostring
        # is deliberately excluded from is_v2_gene_expr below. v1 keeps 3 info columns
        # (gene_id, symbol, name) and falls through to the branch below.
        feature_info  <- feature_table[, c(1:4)]
        feature_table <- feature_table[, -c(1:4)]
        rownames(feature_table) <- as.character(feature_info$gene_id)
      } else {
        rownames(feature_table) <- feature_table$gene_id
        feature_info <- feature_table[,c(1:3)]
        feature_table <- feature_table[,-c(1:3)]
      }
    }

    # get metadata
    mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    query <- paste0("SELECT * FROM ", meta.table)
    metadata <- dbGetQuery(mydb, query)

    # filter by cell and glucose if relevant
    if(grepl("_donor", analysisVar)){
      metadata <- metadata[metadata$cell_type == cell, ]
    } else if(omicsType == "proc_scrna"){
      metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    }

    query <- paste0("SELECT * FROM proc_variable_summary WHERE column='", analysisVar, "'")
    meta.info <- dbGetQuery(mydb, query)
    # Handle variables not in proc_variable_summary (e.g., lipid extraction)
    lipid.vars <- c('tc_weight_recovery', 'tg_weight_recovery', 'fc_weight_recovery', 'ce')
    if(nrow(meta.info) == 0 && analysisVar %in% lipid.vars){
      meta.info <- data.frame(column = analysisVar, type = "cont", stringsAsFactors = FALSE)
    }
    dbDisconnect(mydb)

    # get pathway features from library
    if(funcLib == "hsa_kegg"){
      # KEGG metabolite library: kegg_hsa_met.qs
      # mset.list: named by KEGG pathway IDs (hsa00010 ...) → vectors of KEGG compound IDs
      # path.ids:  named vector  (human-readable name → KEGG pathway ID)
      raw_lib       <- qs::qread(paste0(lib.path, "kegg_hsa_met.qs"))
      name_lookup   <- setNames(names(raw_lib$path.ids), raw_lib$path.ids)  # hsa00010 → "Glycolysis..."
      readable_names <- name_lookup[names(raw_lib$mset.list)]
      readable_names[is.na(readable_names)] <- names(raw_lib$mset.list)[is.na(readable_names)]
      libraryRDS <- list(term = unname(readable_names), sets = raw_lib$mset.list)
    } else {
      libraryRDS <- readRDS(paste0(lib.path, funcLib, ".rds"))
    }
    pathGenes <- libraryRDS$sets[[which(libraryRDS$term == pathName)]]

    saveRDS(feature_info,"savedAnalysis/feature_info.rds")
    saveRDS(feature_table,"savedAnalysis/feature_table.rds")
    saveRDS(meta.info,"savedAnalysis/meta.info.rds")
    saveRDS(pathGenes,"savedAnalysis/pathGenes.rds")
    rcmd <- gsub("tool","local",rcmd)
    write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
     if(!file.exists("savedAnalysis/plotPathwayHeatmap.R")){
       dump("plotPathwayHeatmap", file = "savedAnalysis/plotPathwayHeatmap.R",append=T)
     }     
    
   }else{
       feature_info <- readRDS("feature_info.rds")
     feature_table <- readRDS("feature_table.rds")
     meta.info <- readRDS("meta.info.rds")
      pathGenes <- readRDS("pathGenes.rds")
   }
    meta.type <- meta.info$type

    # get only relevant metadata
    if(omicsType == "proc_scrna"){
      # filter by donors
      if(donors == "subset"){
        donor.list <- readRDS("donors.rds")
        cells.keep <- metadata$cell_id[metadata$record_id %in% donor.list]
        metadata <- metadata[metadata$cell_id %in% cells.keep,]
        if(dim(metadata)[1] < 10){return("RES-NO")}
      }        
      metadata <- metadata[,c("cell_id", analysisVar)]
      colnames(metadata)[1] <- "record_id"
    } else {
      # filter by donors
      if(donors == "subset"){
        donor.list <- readRDS("donors.rds")
        metadata <- metadata[metadata$record_id %in% donor.list,]
        if(dim(metadata)[1] < 10){return("RES-NO")}
      }

      metadata <- metadata[,c("record_id", analysisVar)]
    }

    if(meta.type == "cont"){metadata[,analysisVar] <- as.numeric(metadata[,analysisVar])}
    metadata <- metadata[metadata$record_id %in% colnames(feature_table), ]
    metadata <- metadata[order(metadata[,analysisVar]), ]
    metadata <- metadata[!is.na(metadata[,analysisVar]), ]
    if(nrow(metadata) == 0){ return("RES-NO") }

    # make feature table match metadata table
    # For v2 gene expression (except nanostring): rownames are Ensembl/symbol (not Entrez),
    # so match pathGenes (Entrez IDs from libraries) via the gene_id column.
    # For nanostring and v1: gene_id IS the rowname, match directly.
    is_v2_gene_expr <- (version == "v2" && omicsType %in% c("proc_rnaseq", "proc_scrna")) ||
                       (version == "v2" && grepl("proc_pbrna", omicsType)) ||
                       omicsType %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")
    if(omicsType == "proc_metabolite"){
      # pathGenes are KEGG compound IDs — match by kegg_id column
      keep_rows <- rownames(feature_table)[feature_info$kegg_id %in% pathGenes]
      feature_table <- feature_table[keep_rows, colnames(feature_table) %in% metadata$record_id, drop = FALSE]
    } else if(is_v2_gene_expr){
      # rownames are Ensembl accession / symbol; match pathGenes (Entrez) via gene_id column
      keep_rows <- rownames(feature_table)[feature_info$gene_id %in% pathGenes]
      feature_table <- feature_table[keep_rows, colnames(feature_table) %in% metadata$record_id, drop = FALSE]
    } else {
      feature_table <- feature_table[rownames(feature_table) %in% pathGenes, colnames(feature_table) %in% metadata$record_id]
    }
    if(nrow(feature_table) == 0){ return("RES-NO") }
    feature_table <- feature_table[, match(metadata$record_id, colnames(feature_table)), drop = FALSE]
    if(ncol(feature_table) == 0){ return("RES-NO") }

    # sync feature_info to subsetted rows, then set display labels
    if(omicsType == "proc_metabolite"){
      feat.vec     <- rownames(feature_table)
      hit.inx      <- match(feat.vec, feature_info$compound)
      feature_info <- feature_info[hit.inx, ]
      # rownames are already compound names — no renaming needed
    } else if(version == "v2" && omicsType == "proc_rnaseq"){
      feat.vec     <- rownames(feature_table)      # Ensembl accessions
      hit.inx      <- match(feat.vec, feature_info$accession)
      feature_info <- feature_info[hit.inx, ]
      display_labels <- feature_info$symbol
      display_labels[is.na(display_labels) | display_labels == ""] <- feature_info$gene_id[is.na(display_labels) | display_labels == ""]
      rownames(feature_table) <- display_labels
    } else if(version == "v2" && (grepl("proc_pbrna", omicsType) || omicsType == "proc_scrna")){
      # rownames are unique symbols (deduplicated by p-value) — sync feature_info directly
      feat.vec     <- rownames(feature_table)
      hit.inx      <- match(feat.vec, as.character(feature_info$symbol))
      feature_info <- feature_info[hit.inx, ]
      rownames(feature_table) <- as.character(feature_info$symbol)
    } else if(omicsType %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")){
      feat.vec       <- rownames(feature_table)   # id
      hit.inx        <- match(feat.vec, feature_info$id)
      feature_info   <- feature_info[hit.inx, ]
      display_labels <- feature_info$symbol
      display_labels[is.na(display_labels) | display_labels == ""] <- feature_info$id[is.na(display_labels) | display_labels == ""]
      rownames(feature_table) <- display_labels
    } else {
      feat.vec <- rownames(feature_table)
      hit.inx  <- match(feat.vec, feature_info$gene_id)
      feature_info <- feature_info[hit.inx, ]
      feature_info$symbol[is.na(feature_info$symbol)] <- feature_info$gene_id[is.na(feature_info$symbol)]
      rownames(feature_table) <- feature_info$symbol
    }

    # highlight DEGs
    dea.res <- read.csv("dea_results.csv")
    deg.ind <- which(rownames(feature_table) %in% dea.res$Feature[dea.res$sig != "NS"])
    if(length(deg.ind) > 0){
      rownames(feature_table)[deg.ind] <- paste0(rownames(feature_table)[deg.ind], "***")
    }

    # remove features with too many NAs
    num.vals <- apply(feature_table, 1, function(x){sum(!is.na(x))})
    feature_table <- feature_table[num.vals > 7, ]
    if(nrow(feature_table) == 0){ return("RES-NO") }

    # create heatmap column annotation
    colAnn <- data.frame(Metadata = metadata[,2])
    rownames(colAnn) <- metadata$record_id

    # smooth values
    smooth.df <- feature_table
    smooth.df[is.na(feature_table)] <- min(smooth.df, na.rm = T) # must replace NA with low value for smoothing
    smooth.df <- apply(smooth.df, 1, function(x){
      smth.gaussian(as.numeric(x), window = 0.03, alpha = 2.5, tails = TRUE, na.rm = TRUE)
    }) %>% t() %>% as.data.frame()
    colnames(smooth.df) <- colnames(feature_table)
    smooth.df[is.na(feature_table)] <- NA # add NAs back

    # create heatmap
    if(dim(feature_table)[2] < 50){ show.donorID = TRUE } else { show.donorID = FALSE }

    hm <- pheatmap(smooth.df, show_colnames = show.donorID, annotation_col = colAnn, cluster_cols = FALSE,
                   legend = TRUE, annotation_names_col = FALSE, scale = "row", silent = TRUE, border_color = NA,
                   main = pathName, na_col = "black", color=colorRampPalette(c("#000080", "#FFF300"))(30))

    # plot heatmap
    hm.height <- dim(feature_table)[1]*0.2 + 0.75
    Cairo::Cairo(file = "pathway_heatmap.png", unit="in", res=300, width=10, 
                 height= hm.height, type="PNG", bg="white");
        print(hm)
    dev.off()

      ggplot2::ggsave("pathway_heatmap.svg", plot = hm, width = 10, height =  hm.height, dpi = 300)
    ggplot2::ggsave("pathway_heatmap.pdf", plot = hm, width = 10, height =  hm.height, dpi = 300)
    
    write.csv(smooth.df,"pathway_heatmap.csv")
    return("RES-OK")
}


################################################################################
plotDataAvailHeatmap <- function(){
  
  ## heatmap of data availability
  mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  donors <- dbReadTable(mydb, "donor")
  donors <- donors$record_id
  dbDisconnect(mydb)
  
  avail_df <- data.frame(donor = rep(1, length(donors)),
                         tech = rep(1, length(donors)),
                         gsis = rep(0, length(donors)),
                         peri = rep(0, length(donors)),
                         seahorse = rep(0, length(donors)),
                         ephys = rep(0, length(donors)),
                         prot = rep(0, length(donors)),
                         nano = rep(0, length(donors)),
                         bulkrna = rep(0, length(donors)),
                         pbrna = rep(0, length(donors)),
                         scrna = rep(0, length(donors)))
  
  # get lists of record_ids in each outcomes table
  mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  gsis <- dbGetQuery(mydb, "SELECT record_id FROM gsis")[,1] %>% unique() 
  peri <- dbGetQuery(mydb, "SELECT record_id FROM function_summary")[,1] %>% unique()
  seahorse <- dbGetQuery(mydb, "SELECT record_id FROM seahorse")[,1] %>% unique()
  ephys <- dbGetQuery(mydb, "SELECT record_id FROM ephys_donor")[,1] %>% unique()
  scrna <- dbGetQuery(mydb, "SELECT record_id FROM ephys_cell")[,1] %>% unique()
  dbDisconnect(mydb)
  
  # get lists of record_ids in each bulk omics table
  mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  nano <- dbListFields(mydb, "proc_nanostring_merge")
  bulkrna <- dbListFields(mydb, "proc_rnaseq")
  pbrna <- c(dbListFields(mydb, "proc_pbrna_Alpha"), dbListFields(mydb, "proc_pbrna_Beta")) %>% unique()
  prot <- dbListFields(mydb, "proc_prot")
  dbDisconnect(mydb)
  
  # update data avail matrix
  avail_df$gsis[donors %in% gsis] <- 1
  avail_df$peri[donors %in% peri] <- 1
  avail_df$seahorse[donors %in% seahorse] <- 1
  avail_df$ephys[donors %in% ephys] <- 1
  avail_df$prot[donors %in% prot] <- 1
  avail_df$nano[donors %in% nano] <- 1
  avail_df$bulkrna[donors %in% bulkrna] <- 1
  avail_df$pbrna[donors %in% pbrna] <- 1
  avail_df$scrna[donors %in% scrna] <- 1
  
  avail_df <- t(avail_df)
  
  disp.rownames <- c("Clinical metadata", "Technical metadata", "Static insulin secretion", "Dynamic insulin secretion",
                     "Oxygen consumption", "Electrophysiology", "Proteomics (bulk)", "Nanostring (bulk)", "RNAseq (bulk)",
                     "RNAseq (pseudobulk)", "RNAseq (patchSeq)")
  
  pheatmap(avail_df, legend = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
           color = c("white", "#43A047"), angle_col=90, labels_row = disp.rownames, border_color = NA,
           filename=paste0(other.tables.path, "/display_data/avail_heatmap.png"), width=10, height=5)
}
