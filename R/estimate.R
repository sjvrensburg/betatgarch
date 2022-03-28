###################
# Estimation Code #
###################

# Estimate without imposing the empirical version of the Lyapunov condition
# as a constraint.
estimate_ncon <- function(theta_0, y, f_0, lb = c(
  .Machine$double.eps, 0, 0, 0, 2.1), ub = c(Inf, 1, 1, 1, Inf), ...) {
  # Function to pass to nloptr
  eval_f <- function(theta) nll(y = y, f_0 = f_0, theta = theta)
  # You still have to impose the constraint (a + g) >= 0 => (-a - g) <= 0.
  eval_g <- function(theta) -1 * (theta[2] + theta[4])
  eval_jac_g <- function(theta) c(0, -1, 0, -1, 0)

  nloptr::nloptr(
    x0 = theta_0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g,
    eval_jac_g_ineq = eval_jac_g, ...)
}

# Estimate while imposing the empirical version of the Lyapunov condition
# as a constraint.
estimate_con <- function(theta_0, y, f_0, lb = c(
  .Machine$double.eps, 0, 0, 0, 2.1), ub = c(Inf, 1, 1, 1, Inf), ...) {
  # Function to pass to nloptr
  eval_f <- function(theta) nll(y = y, f_0 = f_0, theta = theta)
  # You still have to impose the constraint (a + g) >= 0 => (-a - g) <= 0.
  eval_g <- function(theta) cnstr(y = y, theta = theta)$objective
  eval_jac_g <- function(theta) cnstr(y = y, theta = theta)$jacobian

  nloptr::nloptr(
    x0 = theta_0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g,
    eval_jac_g_ineq = eval_jac_g, ...)
}

#' Fit a Beta-t-GARCH Model Using Maximum Likelihood Estimation
#'
#' @details
#'   This function estimates the Beta-t-GARCH(1,1) model with leverage of
#'   \insertCite{Harvey2008;textual}{betatgarch}. It can, optionally, maximise
#'   the log-likelihood over a restricted parameter space that satisfies an
#'   empirical version of the Lyapunov condition as in
#'   \insertCite{Blasques2018;textual}{betatgarch}.
#'
#'   NOTE: If you receive the error
#'   "\code{student_t_lpdf: Scale parameter is 0, but must be > 0!}"
#'   then try to increase the lower bound on the parameter nu.
#'
#' @param x numeric vector of observations
#' @param n an optional integer that give the number of observation to use as
#'   the pre-sample.
#' @param theta_0 either NULL or a numerical vector that specifies the initial
#'   guess used by `nloptr` to find optimal parameter estimates.
#' @param lb,ub numeric vectors that specify the upper and lower bounds used by
#'   the optimisation routine.
#' @param constrain logical that indicates whether to optimise the model over a
#'   restricted parameter space that satisfies an empirical version of the
#'   Lyapunov condition.
#' @param opts a list that specifies the algorithm and other options used by
#'   `nloptr` to minimise the negative log-likelihood.
#' @param ... additional arguments to pass to `nloptr`.
#'
#' @return
#'   Returns a list that contains both the output from `nloptr` and the result
#'   of filtering the series using the estimates obtained from `nloptr`.
#'
#' @references
#' \insertRef{Harvey2008}{betatgarch}
#'
#' \insertRef{Blasques2018}{betatgarch}
#'
#' @importFrom Rdpack reprompt
#' @export

fit <- function(x, n = 5, theta_0 = NULL, lb = c(
  .Machine$double.eps, 0, 0, 0, 2.1), ub = c(Inf, 1, 1, 1, Inf),
  constrain = TRUE, opts = list("algorithm" = "NLOPT_LD_SLSQP",
    "xtol_rel" = 1.0e-4, maxeval = 100), ...) {
  # Prepare the data
  prep <- prep_data(x, n)
  f_0 <- prep$f_0
  y <- prep$y
  # If not supplied by the use then set reasonable initial guess and
  # optimisation options.
  if (is.null(theta_0)) {
    theta_0 <- c(w = f_0*0.1, a = 0.05, b = 0.9, g = 0, n = 4)
    }
  # Estimate
  if (constrain) {
    est <- estimate_con(theta_0, y, f_0, lb = lb, ub = ub, opts = opts, ...)
  } else {
    est <- estimate_ncon(theta_0, y, f_0, lb = lb, ub = ub, opts = opts, ...)
  }
  if (est$status < 0 | est$status == 5) warning(est$message)
  names(est$solution) <- c("omega", "alpha", "beta", "gamma", "nu")
  # Filter
  flt <- recursion_lst(y, f_0, est$solution)

  return(list(Estimation = est, Filter = flt, Constrain = constrain))
}
