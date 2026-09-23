# ==============================================================================
# kg_functions/kg_property.R  --  kgProperty: ONE COHORT METADATA VARIABLE, summarised.
#
#   variable    a metadata column (any of the 258 in metadata_sum_raw.csv)
#   subset      "all" | "col<op>value" (the donor row filter, op in = > < >= <=)
#   use_raw     "true" (default) | "false"
#
# WHY THIS EXISTS (user ruling 2026-09-21). "Donor property" is a FRAME — one question
# shape over ALL the metadata — not the handful of columns that happen to be stored on a
# neo4j :Donor node. The engine's op15 could only answer for what the graph held (16 donor
# columns + the numeric phenotype edges), so six variables with values in this very file
# answered NOTHING: hla_a2 (371 donors), collagenase (628), collagenasetype (629),
# fatinfiltration (634), pancreasconsistency (636), final_cluster (625). They are STRINGS,
# and the 2026-08-08 donor-edge load wrote numbers only. The frame's reach was being
# decided by storage, which is the bug this closes: the values live in a CSV, the CSV
# belongs to R, so R answers.
#
# ⚠ RAW, NOT NORM, AND THAT IS THE WHOLE POINT OF THE DEFAULT. metadata_sum_norm.csv is
# log10 on 40 columns (hba1c mean 0.755 there vs 5.784 raw) — a mean read off it is not the
# number a person asked for. The project's rule is already written down: raw is the source
# for VALUES, norm is the source for MODELLING. Only a caller that explicitly wants the
# modelling scale passes use_raw="false".
#
# ⚠ THE TYPE IS INFERRED FROM THE VALUES, never from a list of column names, so a column
# added to the file tomorrow is answered the day it appears (.kgInferType, kg_common.R).
#
# ⚠ A VARIABLE NOBODY CARRIES IS AN ANSWER, NOT AN ERROR. It returns RES-OK with one row
# and N=0, so "no donor has a value for this" can be said plainly. An error code there
# would read as "the question failed", which is a different thing entirely.
#
# Writes kg_property.csv (+ a timestamped copy); returns
#   "RES-OK;<n_rows>" | "RES-NO-META" | "RES-NO-VAR" | "RES-NO-SUBSET"
# ==============================================================================

# ---- the donor subset: ALL FOUR SHAPES op15 ACCEPTS --------------------------
# ⚠ THE FRAME ALREADY PROMISES THEM. `generate_reason_v3.OP_SUBSET_KINDS[15]` is
# {"eq","cmp","pheno","argmax"} — op15 is declared able to take every shape — and the graph side
# implements all four. This function existed with only `col<op>value`, so the other two fell back
# to the graph, where cluster membership and hla_a2 do not exist: MEASURED 2026-09-21 on the F37
# build, 23 of 209 records answered "The answer was not computed …" for data this very file holds
# ("Among high-purity preparations, how many are in cluster 3?").
# ⚠ THE SEMANTICS ARE COPIED FROM op15, NOT INVENTED (`agent._cohort_property`):
#     pheno:<phenotype_id><op><n>   the donor's VALUE for a phenotype, op in > < >= <= (NO "="),
#                                   the id with or without its `local_` prefix
#     top:<phenotype_id>            the donor's LARGEST value among its SIBLINGS (an argmax, not a
#                                   threshold) — the family is the id with its last part replaced
#                                   by "_", e.g. local_X1_EUR -> X1_
# ⚠ AN UNAPPLIABLE FILTER RETURNS NULL, which the caller turns into RES-NO-SUBSET. A filter that
# applies and matches NOBODY returns zero rows, which is an ANSWER. The two must not be confused.
.kgPropertySubset <- function(m, subset){
  s <- trimws(as.character(subset))
  if(!nzchar(s) || identical(s, "all")) return(m)

  pm <- regmatches(s, regexec("^pheno:([A-Za-z0-9_.]+)\\s*(>=|<=|>|<)\\s*(-?[0-9.]+)$", s))[[1]]
  if(length(pm) == 4){
    col <- sub("^local_", "", pm[2]); op <- pm[3]; v <- suppressWarnings(as.numeric(pm[4]))
    if(is.na(v) || !(col %in% names(m))) return(NULL)
    lhs  <- suppressWarnings(as.numeric(as.character(m[[col]])))
    if(all(is.na(lhs))) return(NULL)                 # not a numeric column: the filter cannot apply
    keep <- switch(op, ">" = lhs > v, "<" = lhs < v, ">=" = lhs >= v, "<=" = lhs <= v,
                   rep(FALSE, length(lhs)))
    keep[is.na(keep)] <- FALSE
    return(m[keep, , drop = FALSE])
  }

  tm <- regmatches(s, regexec("^top:([A-Za-z0-9_.]+)$", s))[[1]]
  if(length(tm) == 2){
    pid  <- sub("^local_", "", tm[2])
    stem <- sub("_[A-Za-z]+$", "_", pid)             # X1_EUR -> X1_
    fam  <- names(m)[startsWith(names(m), stem)]
    if(!(pid %in% fam) || length(fam) < 2) return(NULL)
    M <- suppressWarnings(vapply(fam, function(cc) as.numeric(as.character(m[[cc]])),
                                 numeric(nrow(m))))
    if(!is.matrix(M) || all(is.na(M))) return(NULL)
    # the donor's OWN largest sibling, NA-safe; a donor measured for none of the family is dropped
    hit <- apply(M, 1, function(row){
      if(all(is.na(row))) return(NA_character_)
      fam[which.max(replace(row, is.na(row), -Inf))]
    })
    keep <- !is.na(hit) & hit == pid
    return(m[keep, , drop = FALSE])
  }

  .kgDonorSubset(m, s)                                # col<op>value — the shared parser
}


.kgPropertyImpl <- function(variable = "", subset = "all", use_raw = "true"){
  .kgSetPaths()
  variable <- trimws(as.character(variable))
  if(!nzchar(variable)) return("RES-NO-VAR")

  m <- .kgLoadMetaFrame(use_raw = !identical(tolower(trimws(use_raw)), "false"))
  if(is.null(m)) return("RES-NO-META")

  # ⚠ R000 IS NOT A DONOR (user rule 2026-08-08, the same rule the graph refresh follows: the
  # cohort is R001-R650). The file carries 651 rows and R000 holds five real-looking values —
  # diagnosis, diagnosis_computed, puritypercentage, trappedpercentage, totalieq — so leaving it
  # in would answer "of 651" and put every one of those five counts one above the graph's.
  if("record_id" %in% names(m)) m <- m[as.character(m[["record_id"]]) != "R000", , drop = FALSE]

  # the column, matched EXACTLY first and case-insensitively second — the caller sends a
  # resolved column id, and a near-miss must not silently answer about another column.
  col <- if(variable %in% names(m)) variable
         else { hit <- names(m)[tolower(names(m)) == tolower(variable)]
                if(length(hit) == 1) hit else NA_character_ }
  if(is.na(col)) return("RES-NO-VAR")

  # ---- donors: all, or the subset the question asked for -----------------------
  # ⚠ THE ONE PARSER (.kgDonorSubset, kg_common.R), so "bodymassindex>30" means here exactly
  # what it means in kg_corr/kg_assoc. An unappliable filter REFUSES rather than answering
  # the unrestricted question — the failure kg_corr documents at its own subset block.
  subset_note <- "all donors"
  md <- m
  if(!identical(subset, "all") && nzchar(subset)){
    md <- .kgPropertySubset(m, subset)
    if(is.null(md)) return("RES-NO-SUBSET")
    subset_note <- subset
  }
  denom <- nrow(md)

  vals     <- as.character(md[[col]])
  nonblank <- !is.na(vals) & !(vals %in% c("", "NA", "NaN"))
  ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 32))

  .write <- function(out){
    utils::write.csv(out, paste0("kg_property_", safe(col), "_", ts, ".csv"), row.names = FALSE)
    utils::write.csv(out, "kg_property.csv", row.names = FALSE)
    paste0("RES-OK;", nrow(out))
  }

  if(!any(nonblank))
    return(.write(data.frame(Variable = col, Type = "none", Level = NA_character_, N = 0L,
                             Mean = NA_real_, Min = NA_real_, Max = NA_real_,
                             Denominator = denom, Subset = subset_note,
                             stringsAsFactors = FALSE)))

  type <- .kgInferType(vals)
  if(identical(type, "cont")){
    num <- suppressWarnings(as.numeric(vals[nonblank]))
    num <- num[!is.na(num)]
    # a column the type call read as continuous but that holds no coercible number is
    # reported as held-but-unsummarisable rather than as an invented mean of nothing.
    if(!length(num))
      return(.write(data.frame(Variable = col, Type = "none", Level = NA_character_, N = 0L,
                               Mean = NA_real_, Min = NA_real_, Max = NA_real_,
                               Denominator = denom, Subset = subset_note,
                               stringsAsFactors = FALSE)))
    return(.write(data.frame(Variable = col, Type = "cont", Level = NA_character_,
                             N = length(num), Mean = mean(num), Min = min(num), Max = max(num),
                             Denominator = denom, Subset = subset_note,
                             stringsAsFactors = FALSE)))
  }

  # ---- discrete: one row per LEVEL, biggest first ------------------------------
  # ⚠ EVERY LEVEL IS REPORTED, not a top-N. A distribution IS the answer here, and the
  # levels of a metadata column are few (2-7 measured across this file's categoricals).
  tb  <- sort(table(vals[nonblank]), decreasing = TRUE)
  .write(data.frame(Variable = col, Type = "disc", Level = names(tb), N = as.integer(tb),
                    Mean = NA_real_, Min = NA_real_, Max = NA_real_,
                    Denominator = denom, Subset = subset_note,
                    stringsAsFactors = FALSE))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; an R error comes back as "RES-ERR;<fn>;<message>"
# instead of a null -> empty body, and every RES-OK / RES-NO* return passes through.
kgProperty <- function(...) .kgGuard("kgProperty", .kgPropertyImpl, ...)
