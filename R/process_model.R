#' Models the process (infection transmission, recovery, etc.).
#' State variables and beta are given as a vector in which the number of elements are the number of particles
#'
#' The \code{process_model()} does SE1E2IR model integration using Euler method
#' @param tstart start time for simulation with default value of 0
#' @param tend end time with default value of 1 (i.e., one day process of SEI1I2R)
#' @param dt A step size
#' @param params Parameters
#' @param y A vector of a vector of state variables
#' @param beta A vector of beta
#' @export
#' @import data.table
#' @examples
#' pop <- process_model(params, y, tbegin = 0, tend = 1,)
# Process model for simulation --------------------------------------------
process_model <- function (params = NULL,
                           y = NULL,
                           tbegin = 0,
                           tend = 1,
                           dt = 0.2,
                           beta = NULL,
                           stoch = TRUE) {

  y[, c("CE", "CI", "CR")] <- 0 # reset to zero to hold values from tbegin to tend

  S <- y[, "S"]
  E <- y[, "E"]
  P <- y[, "P"]
  A <- y[, "A"]
  I <- y[, "I"]
  R <- y[, "R"]
  daily_infected <- y[, "CE"]
  daily_symptom_onset <- y[, "CI"]
  daily_confirmed <- y[, "CR"]

  N <- S + E + P + A + I + R
  eps <- params[["epsilon"]]
  bp <- params[["bp"]]
  ba <- params[["ba"]]
  fa <- params[["fa"]]
  gam <- params[["gamma"]]
  durP <- (1 / params[["delta"]] - 1 / eps) # residence time in P
  zeta <- 1 / durP # rate from P to A or I

  for (i in seq((tbegin + dt), tend, dt)) {
    # FOI is not adjusted by N because of using time-dependent Rt
    FOI <- beta * (bp * P + ba * A + I)
    # message(paste("length of S =", length(S)))
    if (stoch) {
      # note that FOI is multiplied by multiplying 1/S
      # FOI * S/N in ODE  <=> Bin(S, FOI/N) changes to
      # FOI in ODE  <=> Bin(S, FOI/S) to account for time-dependent Rt
      S_to_E <- rbinom(length(S), S, 1 - exp( - FOI / S * dt))
      E_to_P <- rbinom(length(E), E, 1 - exp( - eps * dt))
      P_to_AI <- rbinom(length(P), P, 1 - exp( - zeta * dt))
      P_to_A <- rbinom(length(P_to_AI), P_to_AI, fa)
      P_to_I <- P_to_AI - P_to_A
      I_to_R <- rbinom(length(I), I, 1 - exp( - gam * dt))
      A_to_R <- rbinom(length(A), A, 1 - exp( - gam * dt))
    }
    else {
      S_to_E <- FOI * dt
      E_to_P <- E * eps * dt
      P_to_A <- P * fa * zeta * dt
      P_to_I <- P * (1 - fa) * zeta * dt
      I_to_R <- I * gam * dt
      A_to_R <- A * gam * dt
    }
    # Process model for SEPAIR
    S <- S - S_to_E
    E <- E + S_to_E - E_to_P
    P <- P + E_to_P - P_to_A - P_to_I
    A <- A + P_to_A - A_to_R
    I <- I + P_to_I - I_to_R
    R <- R + I_to_R + A_to_R
    daily_infected <- daily_infected + S_to_E
    daily_symptom_onset <- daily_symptom_onset + P_to_I
    daily_confirmed <- daily_confirmed + A_to_R + I_to_R

  }

  y[, "S"] <- S
  y[, "E"] <- E
  y[, "P"] <- P
  y[, "A"] <- A
  y[, "I"] <- I
  y[, "R"] <- R
  y[, "CE"] <- daily_infected
  y[, "CI"] <- daily_symptom_onset
  y[, "CR"] <- daily_confirmed

  return(y)
}
