############
# SIMULATE #
############

#' Generate/Simulate From a Beta-t-GARCH Model
#'
#' @details Simulates observations from a Beta-t-GARCH model.
#'
#' @param n number of observations to simulate
#' @param theta numeric vector of parameter values
#' @param f_0 the initial value of the conditional variance
#' @param burnin the length of the burn-in period for the the simulation
#'
#' @returns
#'   Returns a list of vectors. The elements of the list are (in order):
#'   \enumerate{
#'     \item the simulated observations (`"y_t"`),
#'     \item the conditional variance (`"f_t"`),
#'     \item score innovations (`"s_t"`),
#'     \item the contributions to the log-likelihood (`"llik_t"`) and
#'     \item the series of innovations (`"innovations"`).
#'   }
#'
#' @export
generate <- function(n, theta, f_0 = theta[1] / (1 - theta[3]),
                     burnin = ceiling(n * 0.1)) {
  if (is.null(f_0)) f_0 <- theta[1] / (1 - theta[3])
  # Simulate the innovations
  N <- n + burnin
  e <- rugarch::rdist(distribution = "std", n = N, mu = 0, sigma = 1,
                      shape = theta[5])
  # Simulate the observations
  res <- simulate_lst(e, f_0, theta)
  # Remove burnin period
  if (N != n) {
    res$y_t <- res$y_t[-1*(1:burnin)]
    res$f_t <- res$f_t[-1*(1:burnin)]
    res$s_t <- res$s_t[-1*(1:burnin)]
    res$llik_t <- res$llik_t[-1*(1:burnin)]
    res$innovations <- res$innovations[-1*(1:burnin)]
  }
  return(res)
}
