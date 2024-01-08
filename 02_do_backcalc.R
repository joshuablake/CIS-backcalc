suppressMessages(library(dplyr))
library(ggplot2)
library(tidybayes)
library(tidyr)
source(here::here("utils.R"))

########################################################################
## DECONVOLUTION FUNCTIONS
########################################################################

#' Construct a matrix S such that E(P) = S %*% Z
#' 
#' Here, Z is the incidence. and P is the prevalence.
#' See Prevalence Process of thesis for interpretation of the quantities.
#' 
#' @param S a giving the survival function for 1, ..., dmax-1
#' @param size the length of the prevalence vector, which determines the size of the output matrix
#' @return a matrix S of dimension size x size such that E(P) = S %*% Z
construct_matrix_S <- function(S, size) {
  # Pre-conditions
  stopifnot(all(is.finite(S)))
  stopifnot(length(size) == 1)
  stopifnot(round(size) == size)

  # Calculate outputs
  out <- square_matrix_0s(size)
  S <- pad_to_length(S, size)
  for (i in 1:size) {
    # Set first i columns of row i of matrix to f reversed
    out[i, 1:i] <- S[i:1]
  }

  # Some post-conditions
  stopifnot(is_lower_triangular(out))
  stopifnot(all(diag(out) == S[1]))
  stopifnot(all(is.finite(out)))
  stopifnot(
    all.equal(
        (out %*% rep(1, size))[,1],
        cumsum(S),
        check.attributes = FALSE
    )
  )

  return(out)
}

#' Given incidence and a survival function, calculate the prevalence
#' 
#' @param incidence a vector of incidence
#' @param S a vector of survival probabilities
#' @return a vector of prevalence
conv <- function(incidence, S) {
  convolve(incidence, rev(S), type = "o")[seq_along(incidence)]
}

#' Given prevalence and a survival function, calculate the incidence
#' 
#' @param prevalence a vector of prevalence
#' @param S a vector of survival probabilities
#' @param constant_prev_inc assume prevalence is constant historically
#' @return a vector of incidence
deconv <- function(prev, S, constant_prev_inc = TRUE) {
  # Some pre-conditions
  stopifnot(S[1] >= 0.3) # for numeric stability
  stopifnot(all(is.finite(prev) & prev >= 0))
  stopifnot(all(is.finite(S) & S >= 0))
  
  # Remove the prevalence explained by historic incidence
  if (constant_prev_inc) {
    assumed_historic_prev <- prev[1] / sum(S)
    fraction_still_pos_by_time <- S[-1] |>
      rev() |>
      cumsum() |>
      rev() |>
      pad_to_length(length(prev))
    prev_new <- pmax(0, prev - assumed_historic_prev * fraction_still_pos_by_time)
    stopifnot(prev_new <= prev)
    prev <- prev_new
  }
  
  # Do the deconvolution
  S_matrix <- construct_matrix_S(S, length(prev))
  incidence <- forwardsolve(S_matrix, prev)
  
  # Some post-conditions
  stopifnot(length(prev) == length(incidence))
  stopifnot(all.equal(conv(incidence, S), prev))
  return(incidence)
}

########################################################################
## READ DATA AND DECONVOLVE BY STRATA
########################################################################

base_dir = here::here("data/STATS18596/")
S = readr::read_csv(
    file.path(base_dir, "meanS.csv"),
    show_col_types = FALSE
) |>
    arrange(time) |>
    pull("S")
prev = load_prev()
postrat_table = load_poststrat_table()

results = prev |>
  group_by(region, sex, age_group, ethnicityg, .draw) |>
  arrange(daynr, .by_group = TRUE) |>
  mutate(
      incidence = deconv(prev_predict, S),
  ) |>
  ungroup()

saveRDS(
  results,
  here::here("outputs/deconv.rds")
)

########################################################################
## PLOTS
########################################################################
incidence_summary = results |>
    group_by_strata() |>
    median_qi(incidence)
    
poststratify(results, postrat_table, incidence, region, age_group) |>
    group_by(daynr, region, age_group) |>
    median_qi(val) |>
    ggplot(aes(daynr, val, ymin = .lower, ymax = .upper)) +
    geom_lineribbon(alpha = 0.3) +
    facet_grid(region~age_group)
