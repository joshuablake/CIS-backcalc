suppressMessages(library(dplyr))
library(rlang)
library(tidyr)
source(here::here("utils.R"))

day0 = lubridate::ymd("2020-08-31")

rename_regions = function(x) {
    case_match(
        x,
        "0_Eng" ~ "England",
        "1_NE" ~ "North East",
        "2_NW" ~ "North West",
        "3_YH" ~ "Yorkshire",
        "4_EM" ~ "East Midlands",
        "5_WM" ~ "West Midlands",
        "6_EE" ~ "East of England",
        "7_LD" ~ "London",
        "8_SE" ~ "South East",
        "9_SW" ~ "South West",
    )
}

incidence = load_incidence() |>
    mutate(
        date = daynr + day0,
        incidence = pmax(0, incidence),
    )
poststrat_table = load_poststrat_table()

# Poststratify region-age incidence and prevalence
region_age_incidence = poststratify(incidence, poststrat_table, incidence, region, age_group) |>
    mutate(region = rename_regions(region)) |>
    rename(incidence = val)
region_age_prevalence = poststratify(incidence, poststrat_table, prev_predict, region, age_group) |>
    mutate(region = rename_regions(region)) |>
    rename(prevalence = val)
region_age = inner_join(
    region_age_incidence,
    region_age_prevalence,
    by = c("daynr", "date", "region", "age_group", ".draw"),
    unmatched = "error",
    relationship = "one-to-one"
)
rm(region_age_incidence, region_age_prevalence)

# Poststratify region incidence and prevalence
region_incidence = poststratify(incidence, poststrat_table, incidence, region) |>
    bind_rows(
        poststratify(incidence, poststrat_table, incidence) |>
            mutate(region = "0_Eng")
    )  |>
    mutate(
        region = rename_regions(region),
    ) |> 
    rename(incidence = val)
region_prevalence = poststratify(incidence, poststrat_table, prev_predict, region) |>
    bind_rows(
        poststratify(incidence, poststrat_table, prev_predict) |>
            mutate(region = "0_Eng")
    )  |>
    mutate(
        region = rename_regions(region),
    ) |> 
    rename(prevalence = val)
region = inner_join(
    region_incidence,
    region_prevalence,
    by = c("daynr", "date", "region", ".draw"),
    unmatched = "error",
    relationship = "one-to-one"
)
rm(region_incidence, region_prevalence)

# Poststratify age incidence and prevalence
age_incidence = poststratify(incidence, poststrat_table, incidence, age_group) |>
    rename(incidence = val)
age_prevalence = poststratify(incidence, poststrat_table, prev_predict, age_group) |>
    rename(prevalence = val)
age = inner_join(
    age_incidence,
    age_prevalence,
    by = c("daynr", "date", "age_group", ".draw"),
    unmatched = "error",
    relationship = "one-to-one"
)
rm(age_incidence, age_prevalence)

# function which saves its arguments to rds files named after the arguments
save_vars = function(..., dir_name = here::here("for_thesis")) {
    for (x in enquos(...)) {
        saveRDS(eval_tidy(x), file.path(dir_name, paste0(as_label(x), ".rds")))
    }
}
save_vars(
    age,
    region,
    region_age,
)
