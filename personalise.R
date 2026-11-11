## ===========================================================================
##  Final Assessment -- personalised dataset generator
##  Introduction to Data Science 11372 / 11516
##
##  Derives a unique but REAL dataset for each student from their student ID:
##    * a seeded subset of countries (stratified by region)
##    * a seeded analysis window
##    * a seeded set of planted data-quality defects
##
##  Values are genuine OWID figures. Nothing is simulated.
## ===========================================================================

seed_from_id <- function(student_id) {
  id <- toupper(gsub("[^0-9A-Za-z]", "", as.character(student_id)))
  if (nchar(id) == 0L) stop("student_id must contain letters or digits.")
  ch <- utf8ToInt(id)
  h  <- 5381                                   # djb2-style rolling hash
  for (k in ch) h <- (h * 33 + k) %% 2147483647
  as.integer(h %% 1000000000) + 1L
}

## Alternative spellings used to break the name-based join.
NAME_VARIANTS <- list(
  "South Korea"          = c("Korea, South", "Korea, Rep.", "Republic of Korea"),
  "Russia"               = c("Russian Federation", "Russia "),
  "Vietnam"              = c("Viet Nam", "Vietnam, Soc. Rep."),
  "Czechia"              = c("Czech Republic", "Czechia "),
  "Turkey"               = c("Turkiye", "Türkiye", "Turkey, Rep. of"),
  "United States"        = c("United States of America", "USA", " United States"),
  "United Kingdom"       = c("UK", "United Kingdom of Great Britain"),
  "Iran"                 = c("Iran, Islamic Rep.", "Iran (Islamic Republic of)"),
  "Egypt"                = c("Egypt, Arab Rep.", "Egypt "),
  "Democratic Republic of Congo" = c("Congo, Dem. Rep.", "DR Congo"),
  "Bolivia"              = c("Bolivia (Plurinational State of)", "Bolivia "),
  "Tanzania"             = c("Tanzania, United Rep. of", "United Republic of Tanzania"),
  "Syria"                = c("Syrian Arab Republic", "Syria "),
  "Laos"                 = c("Lao PDR", "Lao People's Dem. Rep."),
  "Myanmar"              = c("Burma", "Myanmar (Burma)"),
  "Slovakia"             = c("Slovak Republic", "Slovakia "),
  "Cote d'Ivoire"        = c("Ivory Coast", "Côte d'Ivoire"),
  "Brunei"               = c("Brunei Darussalam", "Brunei "),
  "Venezuela"            = c("Venezuela, RB", "Venezuela (Bolivarian Republic of)"),
  "Moldova"              = c("Republic of Moldova", "Moldova, Rep.")
)

make_my_dataset <- function(student_id,
                            raw_dir = "data_raw",
                            out_dir = "data",
                            quiet   = FALSE) {

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")  # stable on R >= 3.6
  set.seed(seed_from_id(student_id))

  rd <- function(f) read.csv(file.path(raw_dir, f), check.names = FALSE,
                             stringsAsFactors = FALSE)
  emis <- rd("Emissions.csv")
  renw <- rd("Renewables.csv")
  cons <- rd("Consumption.csv")
  ctry <- rd("Countries.csv")

  ## ---- 1. Seeded analysis window ----------------------------------------
  y_start <- sample(2000:2004, 1)
  y_end   <- sample(2021:2023, 1)
  yrs     <- y_start:y_end

  ## ---- 2. Seeded country subset, stratified by region -------------------
  ## Every region is kept. Within a region, larger emitters are more likely
  ## to survive, but none is guaranteed -- so top-10 tables genuinely differ.
  tot <- tapply(emis$annual_co2_mt, emis$iso_code, sum)
  keep <- character(0)
  for (rg in unique(emis$region)) {
    iso_rg <- unique(emis$iso_code[emis$region == rg])
    sz     <- tot[iso_rg]
    rnk    <- rank(sz) / length(sz)                  # 0..1, 1 = largest
    p      <- 0.58 + 0.34 * rnk                      # 0.58 .. 0.92
    taken  <- iso_rg[runif(length(iso_rg)) < p]
    if (length(taken) < 4) {                         # each region needs >= 4
      taken <- union(taken, names(sort(sz, decreasing = TRUE))[1:min(4, length(sz))])
    }
    keep <- c(keep, taken)
  }
  keep <- sort(unique(keep))

  emis <- emis[emis$iso_code %in% keep & emis$year %in% yrs, ]
  cons <- cons[cons$`Country Code` %in% keep & cons$Year %in% yrs, ]
  ctry <- ctry[ctry$countryCode %in% keep, ]
  renw <- renw[, c("Country.Region", intersect(as.character(yrs), names(renw)))]
  keep_names <- unique(emis$location)
  renw <- renw[renw$Country.Region %in% c(keep_names, "World",
                                          "European Union (27)"), ]

  defects <- list(window = c(y_start, y_end), n_countries = length(keep))

  ## ---- 3. Planted defect: negative annual emissions ----------------------
  ## Plausible as a downward revision of a previously published estimate.
  n_neg <- sample(2:4, 1)
  cand  <- which(emis$annual_co2_mt > 5 & emis$year > y_start + 2)
  i_neg <- sample(cand, n_neg)
  emis$annual_co2_mt[i_neg] <- -round(emis$annual_co2_mt[i_neg] *
                                      runif(n_neg, 0.05, 0.30), 3)
  defects$negative_rows <- emis[i_neg, c("location", "year", "annual_co2_mt")]

  ## ---- 4. Planted defect: sentinel values in the oil series --------------
  n_sent <- sample(4:8, 1)
  i_sent <- sample(nrow(emis), n_sent)
  emis$annual_oil_co2_mt[i_sent] <- -999
  defects$sentinel_rows <- emis[i_sent, c("location", "year")]

  ## ---- 5. Planted defect: duplicated observations ------------------------
  n_dup <- sample(3:6, 1)
  i_dup <- sample(setdiff(seq_len(nrow(emis)), c(i_neg, i_sent)), n_dup)
  defects$duplicate_rows <- emis[i_dup, c("location", "year")]
  defects$n_duplicates   <- n_dup
  emis  <- rbind(emis, emis[i_dup, ])
  emis  <- emis[sample(nrow(emis)), ]                # shuffle so dups are not adjacent
  emis  <- emis[order(emis$iso_code, emis$year), ]
  rownames(emis) <- NULL

  ## ---- 6. Planted defect: country-name variants in the wide file ---------
  present <- intersect(names(NAME_VARIANTS), renw$Country.Region)
  n_var   <- min(length(present), sample(4:7, 1))
  chosen  <- if (n_var > 0) sample(present, n_var) else character(0)
  mapping <- character(0)
  for (nm in chosen) {
    alt <- sample(NAME_VARIANTS[[nm]], 1)
    renw$Country.Region[renw$Country.Region == nm] <- alt
    mapping[alt] <- nm
  }
  defects$renamed <- mapping

  ## ---- 7. Planted defect: a block of missing renewables ------------------
  if (nrow(renw) > 5) {
    i_blank <- sample(which(renw$Country.Region %in% keep_names), 1)
    ycols   <- setdiff(names(renw), "Country.Region")
    j0      <- sample(seq_len(max(1, length(ycols) - 5)), 1)
    jj      <- ycols[j0:min(length(ycols), j0 + sample(3:6, 1))]
    renw[i_blank, jj] <- NA
    defects$blank_block <- list(country = renw$Country.Region[i_blank], years = jj)
  }

  ## ---- 8. Planted defect: messy keys and formats in the patchy file ------
  cons$Year <- as.character(cons$Year)
  n_fmt <- max(1, round(nrow(cons) * 0.06))
  i_fmt <- sample(nrow(cons), n_fmt)
  style <- sample(c("fy", "dec", "date"), n_fmt, replace = TRUE)
  cons$Year[i_fmt] <- ifelse(style == "fy",   paste0("FY", cons$Year[i_fmt]),
                      ifelse(style == "dec",  paste0(cons$Year[i_fmt], ".0"),
                                              paste0("01/01/", cons$Year[i_fmt])))
  defects$year_format_rows <- n_fmt

  cons$`Consumption CO2` <- as.character(cons$`Consumption CO2`)
  i_sep <- sample(nrow(cons), max(1, round(nrow(cons) * 0.05)))
  big   <- suppressWarnings(as.numeric(cons$`Consumption CO2`[i_sep])) >= 1000
  i_sep <- i_sep[!is.na(big) & big]
  if (length(i_sep))
    cons$`Consumption CO2`[i_sep] <- formatC(
      as.numeric(cons$`Consumption CO2`[i_sep]), format = "f",
      digits = 3, big.mark = ",")
  defects$thousands_sep_rows <- length(i_sep)

  i_key <- sample(nrow(cons), max(1, round(nrow(cons) * 0.04)))
  cons$`Country Code`[i_key] <- paste0(" ", tolower(cons$`Country Code`[i_key]))
  defects$messy_key_rows <- length(i_key)

  ## ---- 9. Planted defect: a population recorded in thousands -------------
  ## Restricted to countries that also have GDP and GDP/capita, so the error
  ## is always detectable by cross-checking the implied population.
  mid <- which(ctry$popData2022 > 2e6 & ctry$popData2022 < 6e7 &
               !is.na(suppressWarnings(as.numeric(ctry$GDP))) &
               !is.na(suppressWarnings(as.numeric(ctry$`GDP/capita`))) &
               ctry$GDP != "" & ctry$`GDP/capita` != "")
  if (length(mid)) {
    i_unit <- sample(mid, 1)
    defects$unit_error <- list(country = ctry$Country[i_unit],
                               true_pop = ctry$popData2022[i_unit])
    ctry$popData2022[i_unit] <- round(ctry$popData2022[i_unit] / 1000)
  }

  ## ---- 10. Write out ------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(emis, file.path(out_dir, "Emissions.csv"),   row.names = FALSE)
  write.csv(renw, file.path(out_dir, "Renewables.csv"),  row.names = FALSE, na = "")
  write.csv(cons, file.path(out_dir, "Consumption.csv"), row.names = FALSE)
  write.csv(ctry, file.path(out_dir, "Countries.csv"),   row.names = FALSE)

  fp <- dataset_fingerprint(emis, renw, cons, ctry)

  if (!quiet) {
    cat("Student ID   :", student_id, "\n")
    cat("Seed         :", seed_from_id(student_id), "\n")
    cat("Window       :", y_start, "-", y_end, "\n")
    cat("Countries    :", length(keep), "\n")
    cat("Files written to", normalizePath(out_dir), "\n")
    cat("  Emissions.csv  ", nrow(emis), "rows\n")
    cat("  Renewables.csv ", nrow(renw), "rows x", ncol(renw), "cols\n")
    cat("  Consumption.csv", nrow(cons), "rows\n")
    cat("  Countries.csv  ", nrow(ctry), "rows\n")
    cat("FINGERPRINT  :", fp, "\n")
  }

  invisible(list(fingerprint = fp, defects = defects,
                 window = c(y_start, y_end), countries = keep))
}

dataset_fingerprint <- function(emis, renw, cons, ctry) {
  a <- as.integer(round(sum(emis$annual_co2_mt, na.rm = TRUE)))        %% 1000000L
  b <- as.integer(sum(utf8ToInt(paste0(renw$Country.Region, collapse = "")))) %% 10000L
  d <- as.integer(nrow(emis) * 31L + nrow(cons) * 7L + nrow(ctry))     %% 10000L
  sprintf("EM-%06d-%04d-%04d", a, b, d)
}
