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