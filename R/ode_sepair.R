#' An implementation of a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time \code{$\beta$}, average incubation
#' period \code{1/$\sigma$}, and recovery rate \code{$\gamma$}.
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\delta$} rate of transitioning from exposed to infectious state,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' @return A list of changes in S, E, I, and R
#' @export
ode_sepair <- function(t, y, params) {

  S <- y[["S"]]
  E <- y[["E"]]
  P <- y[["P"]]
  A <- y[["A"]]
  I <- y[["I"]]
  R <- y[["R"]]

  ## set beta based on the predefined Rt
  ## first set the duration of infectiousness correct
  ## account for the relative infectiousness of P and A states
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["fa"]]# fraction of asymptomatic state
  bp <- params[["bp"]] # relative infectiousness of pre-symptomatic state
  ba <- params[["ba"]] # relative infectiousness of asymptomatic state
  R0_dur <- ((1 - fa) + fa * ba) * durI + bp * durP

  ## now extract the parameters
  epsilon <- params[["epsilon"]]
  delta <- params[["delta"]]
  gamma <- params[["gamma"]]

  N <- S + E + P + A + I + R

  beta <- params[["beta"]]
  if (!is.null(params[["R0"]])) {
    beta <- params[["R0"]] / R0_dur
  }
  beta <- beta * S / N # transmission rate is adjusted by the prop of susc

  if (params[["time_dep_Rt"]]) {
    # when Rt is used, NO NEED for adjusting for the prop of susc
     beta <- get_Rt(t) / R0_dur
  }

  dS <- 0
  dE <- 0
  dP <- 0
  dA <- 0
  dI <- 0
  dR <- 0
  # Cumulative variables (CE, infected; CI, symptom onset; CR, confirmed)
  # used to extract numbers over a defined period (e.g., daily, weekly, monthly)
  dCE <- 0
  dCI <- 0
  dCR <- 0

  rate_from_p <- 1 / (1/delta - 1/epsilon) # rate of transitioning from state P
  rate_StoE <- beta * (bp * P + ba * A + I)
  rate_EtoP <- epsilon * E
  rate_PtoA <- rate_from_p * fa * P
  rate_PtoI <- rate_from_p * (1 - fa) * P
  rate_AtoR <- gamma * A
  rate_ItoR <- gamma * I

  dS <- - rate_StoE
  dE <- rate_StoE - rate_EtoP
  dP <- rate_EtoP - rate_PtoA - rate_PtoI
  dA <- rate_PtoA - rate_AtoR
  dI <- rate_PtoI - rate_ItoR
  dR <- rate_AtoR + rate_ItoR
  # Cumulative variables (CE, infected; CI, symptom onset; CR, confirmed)
  # used to extract numbers over a defined period (e.g., daily, weekly, monthly)
  dCE <- rate_StoE
  dCI <- rate_PtoI
  dCR <- rate_AtoR + rate_ItoR

  return (list(c(dS, dE, dP, dA, dI, dR, dCE, dCI, dCR)))
}
