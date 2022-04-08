# A class that contains most of the pertinent data relating to
# fitted models.
BetaTGARCH <- R6Class(
  "BetaTGARCH", public = list(
    y = NULL, f_0 = NULL, restricted = NULL, convergence = NULL,
    coef = NULL, sigma = NULL, likelihood.contrib = NULL, likelihood = NULL,
    residuals = NULL, score_residuals = NULL,
    initialize = function(x, coef, f_0 = NULL, n = 5, restricted = TRUE,
                          convergence = TRUE) {
      if (is.null(f_0)) {
        prepped <- prep_data(x, n)
        self$y <- prepped$y
        self$f_0 <- prepped$f_0
      } else {
        self$y <- x
        self$f_0 <- f_0
      }
      # TODO: self$coef should be matrix that contains both the estimates and
      #       the standard errors.
      self$coef <- coef
      self$restricted <- restricted
      self$convergence <- convergence
      rec <- recursion_lst(self$y, self$f_0, self$coef)
      self$sigma <- sqrt(rec$f_t)
      self$likelihood.contrib <- rec$llik_t
      self$likelihood <- sum(self$likelihood.contrib)
      self$residuals <- rec$residuals
      self$score_residuals <- rec$s_t
    },
    show = function() {
      cat("\n#############################", sep = "\n")
      cat("# BETA-t-EGARCH(1, 1) Model #", sep = "\n")
      cat("#############################", sep = "\n")
      if (self$restricted) {
        cat(c("\nModel estimated over a restricted parameter",
              "space that satisfies an empirical version of",
              "the Lyapunov condition.\n"),
            sep = "\n")
      }
      cat(paste("Converged:", ifelse(self$convergence, "TRUE",
                                     "FALSE"), "\n"), sep = "\n")
      # NOTE: If we update self$coef to be a matrix then will have to change.
      coef_mat <- matrix(
        self$coef, ncol = 1, dimnames = list(
          Coefficient = c("omega", "alpha", "beta", "gamma", "nu"),
          Estimate = c("")))
      show(coef_mat)
      cat(sprintf("\nLog-Likelihood: %.4f\n", self$likelihood),
          sep = "\n")
      invisible(self)
    }
  )
)
