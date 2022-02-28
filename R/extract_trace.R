#' Extract results from particle filtering and backward sampling (trace)
#' @param theta Parameter values
#' @param y lstate variables
#' @param data data to calculate likelihood against
#' @param data_type infection, symptom onset, or confirmation
#' @export
#' @examples
#' res <- extract_trace()
extract_trace <- function (params = NULL,
                           y = NULL,
                           data = NULL,
                           data_type = c("infection", "symptom onset", "confirmation"),
                           rep = 1,
                           npart = 1000,
                           tend = 200,
                           dt = 0.2,
                           error_pdf = c("pois", "negbin", "binom"),
                           negbin_size = 15,
                           binom_prob = 0.8,
                           backward_sampling = TRUE,
                           systematic_resampling = FALSE,
                           stoch = TRUE) {

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  nstatevar <- length(y)
  if(is.data.frame(y)) {
    nstatevar <- ncol(y)
  }

  res <- pfilter(params = params,
                 y = y,
                 data = data,
                 data_type = type,
                 npart = npart,
                 tend = tend,
                 dt = dt,
                 error_pdf = error_pdf,
                 negbin_size = negbin_size,
                 binom_prob = binom_prob,
                 backward_sampling = backward_sampling,
                 systematic_resampling = systematic_resampling,
                 stoch = stoch)

  # variables to return = latent state variable + beta
  output <- data.frame(matrix(NA, nrow = tend, ncol = nstatevar + 1))
  names(output) <- c(names(y), "Rt")

  for (nm in names(y)) {
    output[, nm] <- res$trace[[nm]]
  }
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["fa"]]
  ba <- params[["ba"]]
  bp <- params[["bp"]]
  R0_dur <- ((1 - fa) + fa * ba) * durI + bp * durP
  output[, "Rt"] <- res$trace[["beta"]] * R0_dur

  return (output)
}
