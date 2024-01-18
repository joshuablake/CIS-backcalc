source("utils.R")
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
conv <- function(incidence, S) {
  convolve(incidence, rev(S), type = "o")[seq_along(incidence)]
}

base_dir = here::here("data/STATS18596/")
S = readr::read_csv(
    file.path(base_dir, "meanS.csv"),
    show_col_types = FALSE
) |>
    arrange(time) |>
    pull("S")

deconv_exponential <- function(prev, S, expon_prev_inc = TRUE) {
  # Some pre-conditions
  stopifnot(S[1] >= 0.3) # for numeric stability
  stopifnot(all(is.finite(prev) & prev >= 0))
  stopifnot(all(is.finite(S) & S >= 0))
  
  # Remove the prevalence explained by historic incidence
  if (expon_prev_inc) {
    # growth rate of exponential
    r = log(prev[2]) - log(prev[1])
    # incidence up to day 1 based on exponential growth, up to a constant
    historic_log_inc_upto_constant = r * -(1:length(S))
    log_day1_prev_upto_constant = matrixStats::logSumExp(historic_log_inc_upto_constant + log(S))
    # calculate constant
    log_constant = log(prev[1]) - log_day1_prev_upto_constant
    # incidence before day 1
    historic_inc = exp(historic_log_inc_upto_constant[-1] + log_constant)
    # convolve to see who is still positive
    matrix_still_pos = construct_matrix_S(S, length(historic_inc))
    still_pos = pad_to_length(matrix_still_pos %*% historic_inc, length(prev))

    prev_new <- pmax(0, prev - still_pos)
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

prev = load_prev()
results2 = prev |>
  group_by(region, sex, age_group, ethnicityg, .draw) |>
  arrange(daynr, .by_group = TRUE) |>
  mutate(
      incidence = deconv_exponential(prev_predict, S),
  ) |>
  ungroup()

saveRDS(
  results2,
  here::here("outputs/deconv2.rds")
)

age_incidence2 = results2 |>
    mutate(incidence = pmax(0, incidence)) |>
    poststratify(poststrat_table, incidence, age_group)
age_incidence2 |>
    mutate(
        date = daynr + day0,
    ) |>
    ggplot(aes(date, val, colour = age_group, fill = age_group)) +
    stat_lineribbon(alpha = 0.3, .width = 0.95) +
    # facet_wrap(~panel) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
        x = "Date",
        y = "Incidence proportion",
    )
