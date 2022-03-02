#' Simulate Rt
#'
#' The function \code{stoch_seir} implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time ($\beta$), average incubation
#' period (1/$\sigma$), and recovery rate ($\gamma$).
#'
#' @param t Time \code{$\beta$} - transmission rate per unit of time,
#' @return A vector of \code{$R_t$}
#' @export
get_Rt <- function(t) {
  Rt_default <- c(rep(1.2, 100), 0.5*sin(0.1*pi*0:32) + 1.2, rep(0.9, 100))
  #time can be real, not integer, during ODE integration
  Rt <- Rt_default[floor(t+1)]
  return (Rt)
}

