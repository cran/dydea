# Find chaotic regions in data. ---------------------------------------------------------------------------------
find_chaos <- function(data, window_length, skip_window, skip_test01 = 1, test01_thresh = 0.05, find_thresh = 20) {
  #' Find chaotic motions in the data.
  #'
  #' @param data Analyzed data.
  #' @param window_length Length of the window for in which the 0-1 test for chaos will be computed.
  #' @param skip_window Length of the skip of the window moving in the data.
  #' @param skip_test01 Length of the skip to take data for calculation the 0-1 test for chaos in the window.
  #' @param test01_thresh The threshold to decide about motion.
  #' @param find_thresh Precision of found intervals.
  #' @return The list of optimized chaotic motion borders.
  #' @importFrom graphics lines
  #' @importFrom graphics plot
  #' @importFrom graphics points
  #' @export
  #' @examples
  #' # Calculate the logistic map.
  #' cons <- 0.5
  #' data.len <- 17000
  #' chaos.start <- c(5536, 9768)
  #' vec.x <- matrix(cons, data.len, 1)
  #'
  #' vec.x[1] <- (2^0.5)/2
  #' for (i in 2:data.len){
  #'   # x_n+1 = r*x_n(1-x_n)
  #'   vec.x[i] <- 3.7*vec.x[i-1]*(1-vec.x[i-1])
  #' }
  #' vec.x[1:(chaos.start[1]-1)] <-cons
  #' vec.x[(chaos.start[2]+1):data.len] <-cons
  #' tr1 <- seq(from = cons, to = vec.x[chaos.start[1]], length.out = 2001)
  #' tr2 <- seq(from = vec.x[chaos.start[2]], to = cons, length.out = 2001)
  #' vec.x[(chaos.start[1]-2000):chaos.start[1]] <- tr1
  #' vec.x[chaos.start[2]:(chaos.start[2]+2000)] <- tr2
  #'
  #' # Find chaotic intervals in vec.x and plot results.
  #' chaotic_borders <- find_chaos(vec.x, "skip_window" = 1000,
  #'   "window_length" = 3000, "find_thresh" = 300)

  # The 0-1 test for chaos calculated in mowing window.
  test01_res = test_chaos01_mw(data, window_length, skip_window, skip_test01, test01_thresh)

  # Find borders of chaotic motion.
  chaos_borders <- find_chaotic_borders(test01_res)

  # Optimizing the boundaries of chaotic motions.
  chaos_borders_final <- optimize_chaos(find_thresh, test01_thresh, chaos_borders, test01_res, data, skip_window, window_length,
      skip_test01)

  return(do.call(cbind, chaos_borders_final))
}



# Find the borders of chaotic motion (from the output of the function test_chaos01_mw). -----------------------------
find_chaotic_borders <- function(test01_res) {
  ## Find the borders of chaotic motion from the results of test_chaos01_mw.
  ##
  ## @param test01_res The result of the 0-1 test for chaos calculated in mowing window.
  ## @return The list of optimized chaotic motion borders.

  test01res_ <- test01_res$test01_res
  left_borders_temp <- vector(mode = "numeric", length = 0)
  right_borders_temp <- vector(mode = "numeric", length = 0)

  if (test01res_[1] == 1) { # if the first value = 1, the data are chaotic at the beginning
      left_borders_temp <- 1
  }

  for (a in 1:(length(test01res_) - 1)) { # find changes from chaotic behavior
    if (test01res_[a] == 1 & test01res_[a + 1] < 1) {
        right_borders_temp <- c(right_borders_temp, a)
    } else if (test01res_[a] < 1 & test01res_[a + 1] == 1) {
        left_borders_temp <- c(left_borders_temp, a + 1)
    }
  }

  if (test01res_[length(test01res_) - 1] == 1 & test01res_[length(test01res_)] == 1) { # check if the data ends with chaotic behaviour
      right_borders_temp <- c(right_borders_temp, length(test01res_))
  } else if (test01res_[length(test01res_) - 1] < 1 & test01res_[length(test01res_)] == 1) {
      right_borders_temp <- c(right_borders_temp, length(test01res_))
  }

  return(list(left_borders_temp, right_borders_temp))
}



# Optimization the boundaries found by the function find_chaotic_borders. ------------------------------------------------
optimize_chaos <- function(find_thresh, test01_thresh = 0.05, chaos_borders_temp, test01_res, data, skip_window, window_length,
    skip_test01) {
  ## Returns boundaries of found chaotic motion.
  ##
  ## @param find_thresh Precision of found intervals.
  ## @param test01_thresh The threshold to decide about motion.
  ## @param chaos_borders_temp Borders of chaotic motion found by function find_chaotic_borders.
  ## @param test01_res The results of the 0-1 test for chaos in each computed window.
  ## @param data Analyzed data.
  ## @param skip_window Length of the skip of the window moving in the data.
  ## @param window_length Length of the window for in which the 0-1 test for chaos will be computed
  ## @param skip_test01 Length of the skip to take data for calculation the 0-1 test for chaos in the window.
  ## @return The list of optimized chaotic motion borders.

  if ((length(chaos_borders_temp[[1]]) == 0) && (length(chaos_borders_temp[[2]]) == 0)) { # if there is no chaotic interval found
      chaos_borders_final <- chaos_borders_temp  #vector(mode='numeric', length=0)
  } else {
    if ((chaos_borders_temp[[1]][1] == 1) && (chaos_borders_temp[[2]][1] == length(test01_res$test01_res))) { # if there is only chaotic interval found
      chaos_borders_final <- c(1, length(data))
    } else {
      chaos_borders_final <- optimize_chaos_run(find_thresh, test01_thresh, chaos_borders_temp, data,
          skip_window, window_length, skip_test01)
    }
  }
  return(chaos_borders_final)
}


# Optimization the boundaries found by the function find_chaotic_borders. ------------------------------------------------------
optimize_chaos_run <- function(find_thresh, test01_thresh = 0.05, chaos_borders_temp, data, skip_window, window_length, skip_test01) {
  ## Optimization of chaotic motion borders based on bisection method.
  ##
  ## @param find_thresh Precision of found intervals.
  ## @param test01_thresh The threshold to decide about motion.
  ## @param chaos_borders_temp Borders of chaotic motion found by function find_chaotic_borders
  ## @param data Analyzed data.
  ## @param skip_window Length of the skip of the window moving in the data.
  ## @param window_length Length of the window for in which the 0-1 test for chaos will be computed.
  ## @param skip_test01 Length of the skip to take data for calculation the 0-1 test for chaos in the window.
  ## @return The list of optimized chaotic motion borders.

  right_borders_final <- vector(mode = "numeric", length = 0)
  left_borders_final <- vector(mode = "numeric", length = 0)

  length_of_data = length(data)
  intervals_middles <- seq(1, length(data) - window_length, skip_window) + round(window_length/2)

  if (chaos_borders_temp[[1]][1] == 1) { # data starts with chaotic behaviour, then find the end of that behaviour
    left_borders_final <- 1

    reg_int2check <- c(intervals_middles[chaos_borders_temp[[2]][1]], intervals_middles[chaos_borders_temp[[2]][1] +
        1]) # the interval in which change of behaviour will be looking for

    middle <- round((reg_int2check[2] + reg_int2check[1])/2)
    interval <- c(round(middle - window_length/2), round(middle + window_length/2))
    size_interval <- reg_int2check[2] - reg_int2check[1]

    while (size_interval > find_thresh) { # find the change of behaviour (based on bisection method)
      test_series <- data[seq(interval[1], interval[2], skip_test01)]
      res <- Chaos01::testChaos01(test_series, c.gen = "equal", par = "seq")
      if (res > (1 - test01_thresh)) {
        reg_int2check <- c(middle, reg_int2check[2])
      } else {
        reg_int2check <- c(reg_int2check[1], middle)
      }
      middle <- round((reg_int2check[2] + reg_int2check[1])/2)
      interval <- c(round(middle - window_length/2), round(middle + window_length/2))
      size_interval <- reg_int2check[2] - reg_int2check[1]
    }
    right_borders_final <- reg_int2check[2]
  }

  if (chaos_borders_temp[[1]][1] == 1) { # if data starts with chaotic behaviour, then start optimization from the next found chaotic region
      aa <- 2
  } else {
      aa <- 1
  }

  # if data ends with chaotic behaviour, don't optimize the last region in the next while loop
  if (chaos_borders_temp[[2]][length(chaos_borders_temp[[2]])] == length(intervals_middles)) {
      bb <- length(chaos_borders_temp[[2]]) - 1
  } else {
      bb <- length(chaos_borders_temp[[2]])
  }

  # optimization of chaotic intervals (based on bisection method)
  while (aa <= bb) {

    # the right side of the interval
    reg_int2check <- c(intervals_middles[chaos_borders_temp[[2]][aa]], intervals_middles[chaos_borders_temp[[2]][aa] +
        1])
    middle <- round((reg_int2check[2] + reg_int2check[1])/2)
    interval <- c(round(middle - window_length/2), round(middle + window_length/2))
    size_interval <- reg_int2check[2] - reg_int2check[1]
    while (size_interval > find_thresh) {
      test_series <- data[seq(interval[1], interval[2], skip_test01)]
      res <- Chaos01::testChaos01(test_series, c.gen = "equal", par = "seq")
      if (res > (1 - test01_thresh)) {
          reg_int2check <- c(middle, reg_int2check[2])
      } else {
          reg_int2check <- c(reg_int2check[1], middle)
      }
      middle <- round((reg_int2check[2] + reg_int2check[1])/2)
      interval <- c(round(middle - window_length/2), round(middle + window_length/2))
      size_interval <- reg_int2check[2] - reg_int2check[1]
    }
    right_borders_final <- c(right_borders_final, reg_int2check[2])

    # the left side of the interval
    reg_int2check <- c(intervals_middles[chaos_borders_temp[[1]][aa] - 1], intervals_middles[chaos_borders_temp[[1]][aa]])
    middle <- round((reg_int2check[2] + reg_int2check[1])/2)
    interval <- c(round(middle - window_length/2), round(middle + window_length/2))
    size_interval <- reg_int2check[2] - reg_int2check[1]
    while (size_interval > find_thresh) {
      test_series <- data[seq(interval[1], interval[2], skip_test01)]
      res <- Chaos01::testChaos01(test_series, c.gen = "equal", par = "seq")
      if (res < (1 - test01_thresh)) {
        reg_int2check <- c(middle, reg_int2check[2])
      } else {
        reg_int2check <- c(reg_int2check[1], middle)
      }
      middle <- round((reg_int2check[2] + reg_int2check[1])/2)
      interval <- c(round(middle - window_length/2), round(middle + window_length/2))
      size_interval <- reg_int2check[2] - reg_int2check[1]
    }
    left_borders_final <- c(left_borders_final, reg_int2check[1])

    aa = aa + 1
  }

  # if data ends as chaotic, optimize only left border of last found interval
  if (chaos_borders_temp[[2]][length(chaos_borders_temp[[2]])] == length(intervals_middles)) {

    right_borders_final <- c(right_borders_final, length_of_data)

    reg_int2check <- c(intervals_middles[chaos_borders_temp[[1]][length(chaos_borders_temp[[1]])] - 1], intervals_middles[chaos_borders_temp[[1]][length(chaos_borders_temp[[1]])]])
    middle <- round((reg_int2check[2] + reg_int2check[1])/2)
    interval <- c(round(middle - window_length/2), round(middle + window_length/2))
    size_interval <- reg_int2check[2] - reg_int2check[1]
    while (size_interval > find_thresh) {
      test_series <- data[seq(interval[1], interval[2], skip_test01)]
      res <- Chaos01::testChaos01(test_series, c.gen = "equal", par = "seq")
      if (res < (1 - test01_thresh)) {
        reg_int2check <- c(middle, reg_int2check[2])
      } else {
        reg_int2check <- c(reg_int2check[1], middle)
      }
      middle <- round((reg_int2check[2] + reg_int2check[1])/2)
      interval <- c(round(middle - window_length/2), round(middle + window_length/2))
      size_interval <- reg_int2check[2] - reg_int2check[1]
    }
    left_borders_final <- c(left_borders_final, reg_int2check[1])

  }

  return(list(left_borders_final, right_borders_final))
}
