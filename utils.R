########################################################################
## GENERIC UTILITY FUNCTIONS
########################################################################
expit = function(x) 1 / (1 + exp(-x))
logit = function(x) log(x) - log1p(-x)

# Create a square matrix of 0s of size x size
square_matrix_0s <- function(size) {
  matrix(0, nrow = size, ncol = size)
}

# Pad a vector to length_out with fill
pad_to_length <- function(vec, length_out, fill = 0) {
  if (length(vec) >= length_out) return(vec[1:length_out])
  c(vec, rep(fill, length_out - length(vec)))
}

# Check if mat is a lower triangular matrix
is_lower_triangular <- function(mat) {
  all(mat == 0 | lower.tri(mat, diag = TRUE))
}

########################################################################
## DATA FUNCTIONS
########################################################################

base_data_dir = here::here("data/STATS18596/")
load_prev = function() {
    readRDS(file.path(base_data_dir, "predict_thin.rds")) |>
    dplyr::left_join(
        readr::read_csv(
            file.path(base_data_dir, "groups.csv"),
            show_col_types = FALSE
        ),
        by = ".group"
    ) |>
    dplyr::mutate(prev_predict = expit(.predict))
}

group_by_strata = function(x, ...) {
  all_cols = rlang::exprs(region, daynr, sex, age_group, ethnicityg)
  exclude = rlang::ensyms(...)
  stopifnot(all(exclude %in% all_cols))
  include = setdiff(all_cols, exclude)
  return(dplyr::group_by(x, !!!include))
}

load_poststrat_table = function() {
    readr::read_csv(file.path(base_data_dir, "poststrat.csv"), show_col_types = FALSE) |>
        dplyr::mutate(
            region = dplyr::case_match(
                Region_Name,
                "North_East_England"  ~ "1_NE",
                "North_West_England"  ~ "2_NW",
                "Yorkshire"  ~ "3_YH",
                "East_Midlands"  ~ "4_EM",
                "West_Midlands"  ~ "5_WM",
                "East_England"  ~ "6_EE",
                "London"  ~ "7_LD",
                "South_East_England"  ~ "8_SE",
                "South_West_England"  ~ "9_SW",
                .default = NA_character_
            )
        ) |>
        assertr::verify(
            !is.na(region)
            | Region_Name %in% c("Wales", "Scotland", "Northern_Ireland")
        ) |>
        dplyr::filter(!is.na(region))
}
