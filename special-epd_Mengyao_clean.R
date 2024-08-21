library(dplyr)

`%>%` <- magrittr::`%>%`
WDIR <- "/Users/nanqiangbeidiao/Bacon/"
data_raw_path <- "/Users/nanqiangbeidiao/Bacon/files from Mary/"


GPD_METADATA <-
  file.path(data_raw_path,
            "metadata.csv") %>%
  readr::read_csv() %>%
  janitor::clean_names() |>
  dplyr::arrange(entity_name)

aux <- GPD_METADATA |> 
  dplyr::select(entity_name, id_entity) |> 
  dplyr::distinct()

GPD_COUNTS1 <-
  file.path(data_raw_path,
            "samples.csv") %>%
  readr::read_csv() %>%
  janitor::clean_names() |>
  dplyr::inner_join(aux) |>
  dplyr::select(entity_name, everything())

GPD_DATES <-
  file.path(data_raw_path,
            "dates.csv") %>%
  readr::read_csv() %>%
  janitor::clean_names() |>
  dplyr::inner_join(aux) |>
  dplyr::select(entity_name, everything()) |>
  dplyr::mutate(depth_cm = as.double(depth_cm)) |>
  dplyr::arrange(entity_name, depth_cm)  

gpd_subset <- GPD_METADATA %>% 
  dplyr::filter(entity_name%in%c("SOSEDNEE", "STAD1", "STAD2", "STLAWR2", "STLAWR3", "SUMINSK", "	
TANON", "TIANCHI", "TIANCHIYL", "TIKHAN", "TOC11-04", "TOLMACH", "TOM1", "TS"))

gpd_subset_counts <- GPD_COUNTS1 %>%
  dplyr::filter(entity_name %in% gpd_subset$entity_name) %>%
  dplyr::rename(depth = depth_cm)

gpd_subset_dates <- GPD_DATES %>%
  dplyr::filter(entity_name %in% gpd_subset$entity_name)

#epd_subset_dates_coretops <- EPD_DATES_coretops %>%
  #plyr::filter(entity_name %in% epd_subset$entity_name)

# Create vectors of entities----
entities <- unique(gpd_subset$entity_name)

# Extract hiatuses----
gpd_subset_dates_hiatus <- gpd_subset_dates %>%
  dplyr::filter(date_type %>% stringr::str_detect("hiatus|Hiatus"))

# Clean dates ----
charc_cats <- "carbon|AMS|Counting|counting"
marine_cats <- "Marine|Sea|Caspian|marine|sea|caspian"

na_codes <- c(-777777, -888888, -999999)

#if (nrow(epd_subset_dates_coretops) > 0) {
  #epd_subset_dates <- epd_subset_dates %>%
    #dplyr::bind_rows(epd_subset_dates_coretops)
#}

gpd_subset_dates_clean <- gpd_subset_dates %>%
  dplyr::arrange(entity_name, depth_cm) %>%
  dplyr::filter(date_type %>% stringr::str_detect("hiatus|Hiatus", TRUE)) %>%
  dplyr::filter(age_c14 %>% stringr::str_detect("not given", TRUE)) %>%
  dplyr::mutate(age_c14 =  as.double(age_c14)) |>
  dplyr::mutate(age_calib = case_when(is.character(age_calib) ~ NA)) |>
  dplyr::mutate(age_calib =  as.double(age_calib)) |>
  dplyr::left_join(gpd_subset %>% # Append lat and description from the metadata
                     dplyr::select(entity_name, latitude),
                   by = "entity_name") %>%
  dplyr::mutate(age_c14 = dplyr::case_when(is.na(age_c14) ~ NA_real_,
                                           age_c14 %in% na_codes ~ NA_real_,
                                           TRUE ~ age_c14),
                age_calib = dplyr::case_when(is.na(age_calib) ~ NA_real_,
                                           age_calib %in% na_codes ~ NA_real_,
                                           TRUE ~ age_calib),
                age = dplyr::coalesce(age_c14, age_calib),
                error = ifelse(is.na(error) | error <= 0, 1, error),
                cc = dplyr::case_when(stringr::str_detect(date_type, charc_cats) ~ 1,
                                      TRUE ~ 0),
                cc = dplyr::case_when(cc == 1 & latitude <= -15 ~ 3,
                                      cc == 1 &
                                        latitude < 15 &
                                        latitude > -15 ~ 4,
                                      TRUE ~ cc)
  ) %>%
  dplyr::rename(depth = depth_cm) %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::select(entity_name, lab_num, age, error, depth, cc)

# Create mixed calibration curve (cc = 4) ----
#ccdir <- file.path(WDIR, "ccurves")
#ageR::mix_curves(proportion = 0.5, cc1 =1, cc = 1, name = "neotropics.14C", dir = ccdir)

# Create input directories ----
clean_entities <- entities %>%
  purrr::map_chr(function(ent) {
    dates <- gpd_subset_dates_clean %>%
      dplyr::filter(entity_name == ent) %>%
      dplyr::select(lab_num, age, error, depth, cc)
    sample_depths <- gpd_subset_counts %>%
      dplyr::filter(entity_name == ent) %>%
      dplyr::select(depth) %>%
      dplyr::arrange(depth) %>%
      dplyr::mutate(id = seq_along(depth), .before = 1)
    hiatus <- gpd_subset_dates_hiatus %>%
      dplyr::filter(entity_name == ent) %>%
      dplyr::select(depth = depth_cm)
    clean_ent_name <- ageR:::cln_str(ent)
    ageR::create_input(data =
                         list(sample_depths = sample_depths,
                              core = dates,
                              hiatus = hiatus),
                       wdir = file.path(WDIR, "GPD"),
                       entity = clean_ent_name,
                       am = "bacon")
    clean_ent_name
  })

## Create a tibble with the original and clean (without special characters)
## entity names
entities2 <- tibble::tibble(
  original = entities,
  clean = clean_entities)


entities2 <- tibble::tibble(
  original = "TS",
  clean = "TS")
 
#entities2 <- entities2 %>% dplyr::filter(original %in% c("AGHNAGHA"))

# Global parameters ----
MAX_SCENARIOS <- 210
CPUS <- 4 # Change based on your computer specs
## Run the age models (and hope for the best) ----
output <- entities2 %>%
  purrr::pmap(function(original, clean) {
    tryCatch({
      ageR:::msg(original)
      ageR::Bacon2(wdir = file.path(WDIR, "GPD"),
                   entity = clean,
                  # cc = 1,
                   #cc4 = 1
                   #ccdir = ccdir,
                   max_scenarios = MAX_SCENARIOS,
                   #thick_step = 5,
                   #thick_lower = 10,
                   #thick_upper = 100,
                   #acc_lower = 5, 
                  # acc_upper = 100,
                  # acc_step = 10,
                   cpus = CPUS,
                   dry_run = F) %>%
        ageR::pb()
    }, error = function(e) {
      message("ERROR while running `", original, "`: ", e)
      return(NULL)
    })
  })
