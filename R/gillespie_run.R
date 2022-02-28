#' Simulate the model implemented in the Gillespie's direct method
#'
#' @param fun Model
#' @param tend Simulation time
#' @param nrun The number of simulation runs
#' @param y0 Initial states
#' @param params Parameters for the model
#' @return Simulation results
#' @export
gillespie_run <- function(func, tend, nrun, y0, params, report_dt = 1) {
  sim_res <- list()
  for (i in seq_len(nrun)) {
    # cat("i =", i, "\n")
    res <- data.frame(time = 0, t(y0))
    t <- 0
    yt <- y0
    dt <- report_dt # output is recorded at t close to dt,
    while (t < tend) {
      if (params[["time_dep_Rt"]]) {
        params[["R0"]] <- get_Rt(t)
      }
      sim <- func(y = yt, params = params)
      yt <- sim$y
      yt_1 <- c(t, t(yt)) # holds values before report_dt arrives
      t <- t + sim$tau
      if (t > dt){
        res <- rbind(res, yt_1)
        dt <- dt + report_dt # dt
      }
      if (!sim$event_occurred) break
    }
    sim_res[[i]] <- res
    message(paste0("i = ", i))
  }
  return (sim_res)
}
