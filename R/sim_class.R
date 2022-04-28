#' @title BetaTGARCHsim: Simulate from a Beta-t-GARCH(1, 1) model
#' @description  The class provides methods to simulate from either an estimated
#' Beta-t-GARCH(1,1) model or a model specified via its coefficients and an
#' initial value for the conditional variance.
#' @returns
#'   Both of the methods \code{generate} and \code{simulate} returns a list
#'   of vectors. The elements of the list are (in order):
#'   \enumerate{
#'     \item the simulated observations (\code{y_t}),
#'     \item the conditional variance (\code{f_t}),
#'     \item score innovations (\code{s_t}),
#'     \item the contributions to the log-likelihood (\code{llik_t}) and
#'     \item the series of innovations (\code{innovations}).
#'   }
BetaTGARCHsim <- R6::R6Class("BetaTGARCHsim",
  lock_objects = FALSE,
  public = list(
    #' @description Initial value of the conditional variance.
    f_0 = function() private$f_0__,

    #' @description Returns a vector of estimated coefficients/parameters.
    coef = function() private$coef__,

    #' @description
    #' Create a simulation object for a Beta-t-GARCH(1, 1) model where the
    #' coefficients and the initial of the conditional variance are either
    #' extracted from an estimated model object or supplied by the user.
    #' @param model \code{NULL} or a \code{BetaTGARCHfit} object
    #' @param f_0 \code{NULL} or the initial value of the conditional variance
    #' @param coef \code{NULL} or the model coefficients
    #' @details
    #' If \code{model = NULL} then the user must supply values for \code{f_0}
    #' and \code{coef}. Moreover, it is assumed that a supplied model was
    #' estimated.
    initialize = function(model = NULL, f_0 = NULL, coef = NULL) {
      if ("BetaTGARCHfit" %in% class(obj)) {
        private$f_0__ <- model$f_0()
        private$coef__ <- model$coef()
      } else {
        private$f_0__ <- f_0
        private$coef__ <- coef
      }
      if (is.null(private$f_0__) || is.null(private$coef__)) {
        stop("Either supply an estimated model or values for f_0 and coef.")
      }
    },

    #' @description
    #' Generates observations from a given vector of innovations.
    #' @param innovations vector of innovations
    #' @param B an integer that gives the number of periods to discard at the
    #'   start of the generated series, i.e., it is the burn-in period.
    generate = function(innovations, B = ceiling(0.1 * length(innovations))) {
      n <- length(innovations)
      N <- n + B
      if (B >= n) {
        stop("The length of the innovations must exceed the burn-in period.")
      }
      # Simulate the observations
      res <- simulate_lst(innovations, private$f_0__, private$coef__)
      # Remove burnin period
      if (N != n) {
        res$y_t <- res$y_t[-1 * (1:B)]
        res$f_t <- res$f_t[-1 * (1:B)]
        res$s_t <- res$s_t[-1 * (1:B)]
        res$llik_t <- res$llik_t[-1 * (1:B)]
        res$innovations <- res$innovations[-1 * (1:B)]
      }
      return(res)
    },

    #' @description
    #' Simulates a specified number of observations from the specified model.
    #' @param n integer that specifies the number of observation to simulate
    #' @param B an integer that gives the number of periods to discard at the
    #'   start of the generated series, i.e., it is the burn-in period.
    simulate = function(n, B = ceiling(0.1 * n)) {
      N <- n + B
      # Simulate the innovations
      shape <- private$coef__[5]
      innovations <- rugarch::rdist(
        distribution = "std", n = N, mu = 0,
        sigma = 1, shape = shape
      )
      # Generate observations
      self$generate(innovations, B = B)
    }
  ),
  private = list(
    f_0__ = NULL,
    coef__ = NULL
  )
)
