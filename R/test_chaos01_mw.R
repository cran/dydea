# Function for calculating the 0-1 test for chaos in moving windows. --------------------------------------------------
test_chaos01_mw <- function(data, window_length, skip_window, skip_test01, test01_thresh = 0.05) {
  ## Computes the 0-1 test for chaos in the window that is moving in the data.
  ##
  ## @param data Analyzed data.
  ## @param window_length Length of the window for in which the 0-1 test for chaos will be computed
  ## @param skip_window Length of the skip of the window moving in the data.
  ## @param skip_test01 Length of the skip to take data for calculation the 0-1 test for chaos in the window.
  ## @param test01_thresh The threshold to decide about motion
  ##   ## \itemize{
  ##     \item regular motion - results of the 0-1 test for chaos <= test01_thresh
  ##     \item chaotic motion - results of the 0-1 test for chaos >= 1 - test01_thresh
  ##     }
  ## @return The results of the 0-1 test for chaos in each computed window.

  test01_res <- matrix(0, length(seq(1, length(data) - window_length, skip_window)), 1)
  index <- 1
  for (a in seq(1, length(data) - window_length, skip_window)) {
    test_series <- data[seq(a, a + window_length, skip_test01)]
    test01_res[index] <- Chaos01::testChaos01(test_series, c.gen = "equal", par = "seq")
    index <- index + 1
  }

  # set results of chaotic motion to 1 and regular motion to 0
  test01_res <- abs(test01_res)  # absolute values of the 01 test for chaos
  test01_res[test01_res <= test01_thresh] <- 0  # results less than test01_thresh = 0 (regular behaviour)
  test01_res[test01_res >= 1 - test01_thresh] <- 1  # results greater than 1-test01_thresh = 1 (chaotic behaviour)

	intervals <- cbind(seq(1, length(data) - window_length, skip_window), (seq(1, length(data) - window_length, skip_window)+window_length))
	results <- list("test01_res" = test01_res, "intervals" = intervals)

	return(results)
}
