# Calculate f_0 and create vector y_0, y_1, ..., y_T.
prep_data <- function(x, n = 5) {
  f_0 <- var(x[1:n])
  y <- x[n:length(x)]
  list(f_0 = f_0, y = y)
}
