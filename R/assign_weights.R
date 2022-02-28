#' Assign weights (i.e., likelihoods) to each of the value
#'
#' The \code{assign_weights()} calculate the likelihoods assuming case incidences are Poisson-distributed random numbers
#' @param var A vector of state variables from simulation
#' @param t current time
#' @param data data on the times of the infection transmission process (e.g., infection, symptom onset, or confirmation)
#' @param data_type c("infection", "symptom onset", "confirmation")
#' @export
#' @examples
#'  wt <- assign_weights(var = latent_var, t = t, data = data, data_type = type)
assign_weights <- function (var,
                            t,
                            data,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            error_pdf = c("pois", "negbin", "binom"),
                            negbin_size = 5,
                            binom_prob = 0.8) {

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  if (type == "infection") {
    case_expected <- var[, t, "CE"]
    case_data <- round(unlist(data[t, "daily_infected"]))
  }
  else if (type == "symptom onset") {
    case_expected <- var[, t, "CI"]
    case_data <- round(unlist(data[t, "daily_symptom_onset"]))
  }
  else if (type == "confirmation") {
    case_expected <- var[, t, "CR"]
    case_data <- round(unlist(data[t, "daily_confirmed"]))
  }


  if (!is.na(case_data)) {
    expected_val <- pmax(0, case_expected)# case_expected is a vector of length npart
    if (error_pdf == "pois"){
      log_lik <- dpois(round(case_data), lambda = expected_val, log = T)
    }
    else if (error_pdf == "negbin") {
      log_lik <- dnbinom(round(case_data), size = negbin_size, mu = expected_val, log = T)
    }
    else if (error_pdf == "binom") {
      log_lik <- dbinom(round(case_data), size = round(expected_val), prob = binom_prob, log = T)
    }
  }
  else {
    log_lik <- -Inf
  }
  # cat("sum(is.na(log_lik)) =", sum(is.na(log_lik)), ", sum(log_lik) =", sum(log_lik),"\n")
  return (exp(log_lik)) # convert to normal probability
}
