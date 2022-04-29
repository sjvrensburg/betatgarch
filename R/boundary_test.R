boundary_test <- function(y, theta) {
  upsilon <- elc(y = y, theta = theta)
  upsilon <- sqrt(length(upsilon)) * upsilon
  lmFit <- lm(upsilon~1)
  res <- lmtest::coeftest(lmFit, df = Inf, vcov. = sandwich::NeweyWest)
  return(res)
}
