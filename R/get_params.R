#' A routine to extract outputs from a stochastic SEIR model
#'
#'
#' @return Sum of squared distance between the model and the data
#' @export
get_params <- function(){
  R0 <- 2
  epsilon <- 1 / 3
  delta <- 1 / 5
  rate_isol <- 1/2.5
  gamma <- 1 / 20
  fA <- 0.3
  bA <- 1
  bP <- 1
  xP <- 0
  xA <- 1

  y0 <- list(S = 1e7, E = 0, P = 0, A = 0, I = 10, X = 0, R = 0, CI = 0)

  params <- list()
  params$init <- y0
  params$R0 <- R0
  params$epsilon <- epsilon
  params$delta <- delta
  params$rate_isol <- rate_isol
  params$gamma <- gamma
  params$fA <- fA
  params$bP <- bP
  params$bA <- bA
  params$xP <- xP
  params$xA <- xA

  # ode stuff
  # params$nsteps <- 1e3
  params$tau <- 0.25
  params$ndays <- 200
  return (params)
}
