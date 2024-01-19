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

region_age_incidence = poststratify(incidence, poststrat_table, incidence, region, age_group)

region_incidence = poststratify(incidence, poststrat_table, incidence, region) |>
    bind_rows(
        poststratify(incidence, poststrat_table, incidence) |>
            mutate(region = "0_Eng")
    )  |>
    mutate(
        region = rename_regions(region),
    )
age_incidence = poststratify(incidence, poststrat_table, incidence, age_group)

# function which saves its arguments to rds files named after the arguments
save_vars = function(..., dir_name = here::here("for_thesis")) {
    for (x in enquos(...)) {
        saveRDS(eval_tidy(x), file.path(dir_name, paste0(as_label(x), ".rds")))
    }
}
save_vars(
    age_incidence,
    region_incidence,
    region_age_incidence,
)
