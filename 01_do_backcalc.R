suppressMessages(library(dplyr))
library(tidyr)

########################################################################
## UTILITY FUNCTIONS
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
results = readRDS(file.path(base_dir, "predict_thin.rds")) |>
    group_by(.draw, .group, region) |>
    mutate(
        prevalence = expit(.predict),
        incidence = deconv(prevalence, S),
    )
group_lookup = readr::read_csv(
    file.path(base_dir, "groups.csv"),
    show_col_types = FALSE
)
postrat_table = readr::read_csv(
    file.path(base_dir, "postrat.csv"),
    show_col_types = FALSE
)

results |>
    filter(incidence < 0)
    left_join(group_lookup, by = ".group") |>
    ggplot() +
    geom_line(aes(time, incidence, .))

########################################################################
## OLD FUNCTIONS (PROBABLY REMOVE)
########################################################################
# deconv(c(100construct_matrix_S.5, 1, 1))

plot <- expand_grid(
  forward_model = tbl_prev$source |> unique(),
  backward_model = tbl_survival$source |> unique(),
) |>
  rowwise() |>
  mutate(
    prevalence = list(
      tbl_prev |>
        filter(source == forward_model, date >= truncation_date) |>
        arrange(date) |>
        `$`(prevalence)
    ),
    f = list(
      tbl_survival |>
        filter(source == backward_model, time == round(time)) |>
        arrange(time) |>
        `$`(value)
    ),
    zero = list(deconv(prevalence, f)),
    constant = list(deconv(prevalence, f, TRUE)),
    date = list(filter(tbl_rtm_out, date >= truncation_date)$date),
  ) |>
  unnest(c(zero, constant, date)) |>
  pivot_longer(c(zero, constant), names_to = "assumption", values_to = "incidence") |>
  # bind_rows(
  #   tbl_prev |>
  #     arrange(date) |>
  #     group_by(source) |>
  #     summarise(
  #       forward_model = source,
  #       backward_model = "prev / 8",
  #       incidence = prevalence / 8,
  #       date,
  #       .groups = "drop"
  #     )
# ) |>
ggplot(aes(x = date, y = incidence, colour = backward_model)) +
  geom_line() +
  facet_grid(rows = vars(forward_model), cols = vars(assumption)) +
  geom_line(aes(x = date, y = incidence, colour = "true"), data = tbl_rtm_out) +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  ylim(0, 150000) +
  labs(
     x = "Date",
     y = "Incidence",
     colour = "Backwards model"
  )

if (FALSE) {
  ggsave(
    "/home/joshuab/COVID/ons-incidence/documents/feb-writeup/figures/deterministic-simulation-results.pdf",
    width = 19, height = 20, units = "cm"
  )
}

############################################################################

## Creates a function S(t) that returns the probability of a positive swab t - 1 days after infection
generate_S <- function (spec = 0, sens_curve = read_median()) {
  times_use_sens <- seq_along(sens_curve) - 1
  prob_swab_pos_given_inf_time <- function(t) {
    indices <- t + 1
    not_valid <- which(indices <= 0 | indices > length(sens_curve))
    indices[not_valid] <- NA
    out <- sens_curve[indices]
    out[not_valid] <- spec
    return(out)
  }
}

deconv_deterministic <- pracma::deconv(c(tbl_det_prev$B, rep(0, 30)), f(0:30))
conv(deconv_deterministic$q, f(0:30))
length(tbl_det_prev$B)

tbl_inc |>
  mutate(deconv = deconv_deterministic$q) |>
  pivot_longer(-t) |>
  ggplot(aes(x = t, y = value, colour = name)) +
  geom_line()