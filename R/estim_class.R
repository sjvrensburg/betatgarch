#' @title BetaTGARCHfit: Estimation And Related Methods
#' @description  The class provides a method to estimate the Beta-t-GARCH(1,1)
#' model and access estimates, standard errors, etc.
BetaTGARCHfit <- R6::R6Class("BetaTGARCHfit",
  lock_objects = FALSE,
  public = list(
    #' @field data original data passed to BetaTGARCHfit
    data = NULL,
    #' @field matcoef matrix of estimates, robust standard errors, robust t
    #'   statistics and associated p-values
    matcoef = NULL,
    #' @field se estimated standard errors
    se = NULL,
    #' @field nlopt_res the raw results returned by \code{NLOPT}
    nlopt_res = NULL,
    #' @field boolean that indicates whether estimation was successful
    convergence = FALSE,

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
        private$n <- 1
        private$f_0 <- f_0
        private$x <- zoo::coredata(y)
      } else {
        prep <- prep_data(zoo::coredata(y), n = n)
        private$n <- n
        private$f_0 <- prep$f_0
        private$x <- prep$y
      }
      self$data <- y
      private$idx <- zoo::index(y)
      invisible(self)
    },

    #' @description
    #' Returns the number of observations used in estimation. This excludes the
    #' \code{n} pre-sample observations used to initialise the filter.
    nobs = function() length(private$x) - 1,

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
      private$restrict <- restrict
      # Define the objective, gradient and constraints.
      eval_f <- function(x) private$eval_f(x)
      eval_g <- function(x) private$eval_g(x)
      eval_jac_g <- function(x) private$eval_jac_g(x)
      # Initial values...
      if (is.null(x0)) {
        # TODO: Consider changing this to use estimates from a GJR-GARCH.
        x0 <- c(w = private$f_0*0.1, a = 0.05, b = 0.9, g = 0, n = 4)
      }
      # Perform the optimisation
      self$nlopt_res <- nloptr::nloptr(
        x0 = x0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g,
        eval_jac_g_ineq = eval_jac_g, opts = opts, ...
      )
      self$convergence <- TRUE
      if (self$nlopt_res$status < 0 | self$nlopt_res$status == 5) {
        warning(self$nlopt_res$message)
        self$convergence <- FALSE
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
        ans <- recursion_lst(private$x, private$f_0, self$coef())
      }
      if (!is.null(parm)) {
        ans <- recursion_lst(private$x, private$f_0, parm)
      }
      idx <- private$idx[private$n:length(self$data)]
      ans <- lapply(ans, function(x) {
        y <- zoo::as.zoo(x)
        zoo::index(y) <- idx
        return(y)
      })
      private$f_t <- ans$f_t[-1]
      private$score_resid <- ans$s_t[-1]
      private$stnd_resid <- ans$residuals[-1]
      private$llik_t <- ans$llik_t[-1]
      return(ans)
    },

    logLik = function() {
      # TODO: If we change to minimising the average negative log-lik then we
      # must multiply by self.nobs()
      ans <- NULL
      if (!is.null(self$nlopt_res)) ans <- -1 * self$nlopt_res$objective
      return(ans)
    },

    #' @description Returns a vector of estimated coefficients/parameters.
    coef = function() {
      ans <- NULL
      if (!is.null(self$nlopt_res)) {
        ans <- self$nlopt_res$solution
        names(ans) <- c("omega", "alpha1", "beta1", "gamma1", "shape")
      }
      return(ans)
    },

    #' @description Calculates the robust variance-covariance matrix of a fitted
    #'  Beta-t-GARCH(1,1) model. This method returns matrix of estimated
    #'  covariances between model parameters.
    vcov = function(...) {
      if (is.null(self$coef())) {
        msg <- paste("Please estimate the model or filter with or filter with",
                     "user-supplied model parameters.")
        stop(msg)
      }
      H <- numDeriv::jacobian(
        function(x) -1*nll(private$x, private$f_0, x)$gradient,
        self$coef(), ...)
      G <- numDeriv::jacobian(
        function(x) recursion_lst(private$x, private$f_0, x)$llik_t[-1],
        self$coef(), ...)
      rob <- solve(H) %*% t(G) %*% G %*% solve(H)
      rownames(rob) <- names(self$coef())
      colnames(rob) <- names(self$coef())
      self$se <- sqrt(diag(rob))
      return(rob)
    },

    #' @description Returns a vector of standardised residuals from the
    #'  estimated model.
    residuals = function() private$stnd_resid,

    sigma = function() {
      if (is.null(private$f_t)) {
        msg <- paste("Please estimate the model or filter with or filter with",
                     "user-supplied model parameters.")
        stop(msg)
      }
      return(sqrt(private$f_t))
    },

    show = function() {
      cat("\n############################", sep = "\n")
      cat("# BETA-t-EGARCH(1,1) Model #", sep = "\n")
      cat("############################", sep = "\n")
      if (private$restrict) {
        cat(c("\nModel estimated over a restricted parameter",
              "space that satisfies an empirical version of",
              "the Lyapunov condition.\n"),
            sep = "\n")
      }
      cat(paste("Converged:", ifelse(self$convergence, "TRUE",
                                     "FALSE"), "\n"), sep = "\n")

      cat(sprintf("\nLog-Likelihood: %.4f\n", self$logLik()),
          sep = "\n")

      matcoef <- cbind(" Estimate" = self$coef(),
                       " Std. Error" = sqrt(diag(self$vcov())))
      matcoef <- cbind(matcoef, " t value" = matcoef[, 1] / matcoef[, 2])
      self$matcoef <- cbind(matcoef, "Pr(>|t|)" = pnorm(abs(matcoef[, 3]),
                                                        lower.tail = FALSE))
      show(self$matcoef)

      invisible(self)
    },
    #' @description Summary of Estimation Results
    summary = function() invisible(self)
  ),

  private = list(
    # Private fields
    n = NULL,
    idx = NULL,
    f_0 = NULL,
    x = NULL,
    f_t = NULL,
    llik_t = NULL,
    score_resid = NULL,
    stnd_resid = NULL,
    restrict = TRUE,
    # Private methods
    eval_f = function(parm) nll(y = private$x, f_0 = private$f_0, theta = parm),
    eval_g = function(parm) {
      if (private$restrict) {
        ans <- cnstr(y = private$x, theta = parm)$objective
      } else {
        ans <- -1 * (parm[2] + parm[4])
      }
      return(ans)
    },
    eval_jac_g = function(parm) {
      ans <- c(0, -1, 0, -1, 0)
      if (private$restrict) {
        ans <- cnstr(y = private$x, theta = parm)$jacobian
      }
      return(ans)
    }
  )
)
