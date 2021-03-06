#' NASDAQ Composite Index (IXIC) Monthly Returns
#'
#' A dataset that contains the monthly (percentage) log returns to
#' the NASDAQ Composite Index from 1 February 1985 to 1 April 2016.
#' Returns calculated as \deqn{r_t = 100\log{\frac{p_t}{p_{t-1}}},}
#' where \eqn{p_t} is the adjusted closing price at time period
#' \eqn{t.} Data obtained from Yahoo Finance.
#'
#' @format A zoo object with 375 observations.
#'
#' @source \url{https://finance.yahoo.com/}
"nasdaq"
