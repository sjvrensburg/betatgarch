#' @title BetaTGARCHfit: Estimation And Related Methods
#' @description  The class provides a method to estimate the Beta-t-GARCH(1,1)
#' model and access estimates, standard errors, etc.
BetaTGARCHfit <- R6::R6Class("BetaTGARCHfit",
  lock_objects = FALSE,
  public = list(
    #' @description Original data passed to BetaTGARCHfit
    data = function() private$data__,

    #' @description Initial value of the conditional variance.
    f_0 = function() private$f_0__,

    #' @description Matrix of estimates, robust standard errors, robust t
    #'   statistics and associated p-values
    #' @details
    #' Note that the p-values are expressed \eqn{P(Z \ge |t|)}. Therefore, a
    #' parameter is significantly different from nought if the p-value is less
    #' than \eqn{\alpha.}
    matcoef = function() private$matcoef__,

    #' @description Estimated standard errors
    se = function() private$se__,

    #' @description The raw results returned by \code{NLOPT}
    nlopt_res = function() private$nlopt_res__,

    #' @description Boolean that indicates whether estimation was successful
    convergence = function() private$convergence__,

    #' @description
    #' Create a new \code{BetaTGARCHfit} object.
    #' @param y a numeric vector or univariate zoo/xts time series
    #' @param n positive integer giving the size of the pre-sample
    #' @param f_0 optional, value of conditional variance with which to
    #'   initialise the filter. Default is \code{NULL}.
    #' @details
    #' The first \code{n} observations of \code{y} are used to estimate the
    #' initial value of the conditional variance, \code{f_0}. The first
    #' \code{n} observations also do not enter the likelihood calculation. If
    #' \code{f_0} is given then the value of \code{n} is ignored.
    initialize = function(y, n = 5L, f_0 = NULL) {
      if ((is.null(n) | n <= 0) & is.null(f_0)) {
        stop("n must be positive if f_0 is NULL")
      }
      if (!is.null(f_0)) {
        warning("n not used since f_0 given")
        private$n__ <- 1
        private$f_0__ <- f_0
        private$x__ <- zoo::coredata(y)
      } else {
        prep <- prep_data(zoo::coredata(y), n = n)
        private$n__ <- n
        private$f_0__ <- prep$f_0
        private$x__ <- prep$y
      }
      private$data__ <- y
      private$idx__ <- zoo::index(y)
      invisible(self)
    },

    #' @description
    #' Returns the number of observations used in estimation. This excludes the
    #' \code{n} pre-sample observations used to initialise the filter.
    nobs = function() length(private$x__) - 1,

    #' @description Maximum Likelihood Estimation of The Beta-t-GARCH(1,1) Model
    #' @param x0 initial values passed to \code{nloptr}
    #' @param restrict logical that indicates whether to optimise the model over
    #'   a restricted parameter space that satisfies an empirical version of the
    #'   Lyapunov condition.
    #' @param lb,ub vectors that specify the upper and lower bounds for model
    #'   parameters
    #' @param opts a list that specifies the algorithm and other options used by
    #' @param ... additional arguments passed to \code{nloptr}.
    #'   \code{nloptr} to minimise the negative log-likelihood.
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
    #' @references
    #' \insertRef{Harvey2008}{betatgarch}
    #' \insertRef{Blasques2018}{betatgarch}
    #' @importFrom Rdpack reprompt
    fit = function(x0 = NULL, restrict = TRUE, lb = c(
                     .Machine$double.eps, 0, 0, -1, 2.1
                   ), ub = c(Inf, 1, 1, 1, Inf),
                   opts = list(
                     algorithm = "NLOPT_LD_AUGLAG",
                     xtol_rel = 1.0e-12, maxeval = 1000,
                     local_opts = list(
                       algorithm = "NLOPT_LD_LBFGS",
                       xtol_rel = 1.0e-10, maxeval = 500
                     )
                   ), ...) {
      private$restrict__ <- restrict
      # Define the objective, gradient and constraints.
      eval_f <- function(x) private$eval_f(x)
      eval_g <- function(x) private$eval_g(x)
      eval_jac_g <- function(x) private$eval_jac_g(x)
      # Initial values...
      if (is.null(x0)) {
        x0 <- c(w = private$f_0__ * 0.1, a = 0.05, b = 0.9, g = 0, n = 4)
      }
      # Perform the optimisation
      private$nlopt_res__ <- nloptr::nloptr(
        x0 = x0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g,
        eval_jac_g_ineq = eval_jac_g, opts = opts, ...
      )
      private$convergence__ <- TRUE
      if (private$nlopt_res__$status < 0 | private$nlopt_res__$status == 5) {
        warning(private$nlopt_res__$message)
        private$convergence__ <- FALSE
      }
      self$coef()
      self$vcov()
      self$show()
      self$filter()
      invisible(self)
    },

    #' @description
    #' Filter with the estimated model. This method will return a list with the
    #' following elements:
    #'   \enumerate{
    #'     \item each observation's contribution to the log-likelihood (`"llik_t"`),
    #'     \item the conditional variances (f),
    #'     \item score residuals (s_t) and
    #'     \item standardised residuals (residuals).
    #'   }
    #' Note that this method will update several other fields.
    #' @param parm parameter vector. If \code{NULL} (default) then it uses
    #'   estimated values. Supplying a parameter vector is not recommended since
    #'   \code{filter} will overwrite several fields of the object.
    filter = function(parm = NULL) {
      if (is.null(parm) & is.null(self$coef())) {
        stop("Either supply parameter vector of estimate the model.")
      }
      if (is.null(parm) & !is.null(self$coef())) {
        ans <- recursion_lst(private$x__, private$f_0__, self$coef())
      }
      if (!is.null(parm)) {
        ans <- recursion_lst(private$x__, private$f_0__, parm)
      }
      idx <- private$idx__[private$n__:length(private$data__)]
      ans <- lapply(ans, function(x) {
        y <- zoo::as.zoo(x)
        zoo::index(y) <- idx
        return(y)
      })
      private$f_t__ <- ans$f_t[-1]
      private$score_resid__ <- ans$s_t[-1]
      private$stnd_resid__ <- ans$residuals[-1]
      private$llik_t__ <- ans$llik_t[-1]
      return(ans)
    },

    #' @description Returns the log-likelihood of the estimated model.
    logLik = function() {
      ans <- NULL
      if (!is.null(private$nlopt_res__)) ans <- -1 * private$nlopt_res__$objective
      return(ans)
    },

    #' @description Returns a vector of estimated coefficients/parameters.
    coef = function() {
      ans <- NULL
      if (!is.null(private$nlopt_res__)) {
        ans <- private$nlopt_res__$solution
        names(ans) <- c("omega", "alpha1", "beta1", "gamma1", "shape")
      }
      return(ans)
    },

    #' @description
    #' Calculates the robust variance-covariance matrix of a fitted
    #' Beta-t-GARCH(1,1) model. This method returns a matrix of estimated
    #' covariances between model parameters.
    #'
    #' @param ... additional arguments passed to \code{numDeriv::jacobian}
    vcov = function(...) {
      if (is.null(self$coef())) {
        msg <- paste(
          "Please estimate the model or filter with or filter with",
          "user-supplied model parameters."
        )
        stop(msg)
      }
      H <- numDeriv::jacobian(
        function(x) -1 * nll(private$x__, private$f_0__, x)$gradient,
        self$coef(), ...
      )
      G <- numDeriv::jacobian(
        function(x) recursion_lst(private$x__, private$f_0__, x)$llik_t[-1],
        self$coef(), ...
      )
      rob <- solve(H) %*% t(G) %*% G %*% solve(H)
      rownames(rob) <- names(self$coef())
      colnames(rob) <- names(self$coef())
      private$se__ <- sqrt(diag(rob))
      return(rob)
    },

    #' @description Returns a vector of standardised residuals from the
    #'  estimated model.
    residuals = function() private$stnd_resid__,

    #' @description Returns a vector of the estimated conditional standard
    #'  deviation.
    sigma = function() {
      if (is.null(private$f_t__)) {
        msg <- paste(
          "Please estimate the model or filter with or filter with",
          "user-supplied model parameters."
        )
        stop(msg)
      }
      return(sqrt(private$f_t__))
    },

    #' @description Prints a summary of estimation results to the screen.
    show = function() {
      cat("\n############################", sep = "\n")
      cat("# BETA-t-EGARCH(1,1) Model #", sep = "\n")
      cat("############################", sep = "\n")
      if (private$restrict__) {
        cat(c(
          "\nModel estimated over a restricted parameter",
          "space that satisfies an empirical version of",
          "the Lyapunov condition.\n"
        ),
        sep = "\n"
        )
      }
      cat(paste("Converged:", ifelse(private$convergence__, "TRUE",
        "FALSE"
      ), "\n"), sep = "\n")

      cat(sprintf("\nLog-Likelihood: %.4f\n", self$logLik()),
        sep = "\n"
      )

      matcoef <- cbind(
        " Estimate" = self$coef(),
        " Std. Error" = sqrt(diag(self$vcov()))
      )
      matcoef <- cbind(matcoef, " t value" = matcoef[, 1] / matcoef[, 2])
      private$matcoef__ <- cbind(matcoef, "Pr(>|t|)" = pnorm(abs(matcoef[, 3]),
        lower.tail = FALSE
      ))
      show(private$matcoef__)

      invisible(self)
    }
  ),
  private = list(
    # Private fields
    n__ = NULL,
    idx__ = NULL,
    f_0__ = NULL,
    x__ = NULL,
    f_t__ = NULL,
    llik_t__ = NULL,
    score_resid__ = NULL,
    stnd_resid__ = NULL,
    restrict__ = TRUE,
    data__ = NULL,
    matcoef__ = NULL,
    se__ = NULL,
    nlopt_res__ = NULL,
    convergence__ = FALSE,
    # Private methods used by the method `fit`
    eval_f = function(parm) nll(y = private$x__, f_0 = private$f_0__, theta = parm),
    eval_g = function(parm) {
      if (private$restrict__) {
        ans <- cnstr(y = private$x__, theta = parm)$objective
      } else {
        ans <- -1 * (parm[2] + parm[4])
      }
      return(ans)
    },
    eval_jac_g = function(parm) {
      ans <- c(0, -1, 0, -1, 0)
      if (private$restrict__) {
        ans <- cnstr(y = private$x__, theta = parm)$jacobian
      }
      return(ans)
    }
  )
)
