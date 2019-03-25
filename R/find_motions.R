# Find chaotic and regular regions in data. ------------------------------------------------------------------------
find_motions <- function(data, window_length, skip_window, skip_test01 = 1, test01_thresh = 0.05, find_thresh = 20) {
  #' Find regular and chaotic motions in the data and plots the results.
  #'
  #' @param data Analyzed data.
  #' @param window_length Length of the window for in which the 0-1 test for chaos will be computed
  #' @param skip_window Length of the skip of the window moving in the data.
  #' @param skip_test01 Length of the skip to take data for calculation the 0-1 test for chaos in the window.
  #' @param test01_thresh The threshold to decide about motion.
  #' @param find_thresh Precision of found intervals.
  #' @return The list of optimized regular and chaotic motion borders.
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
  #' # Find chaotic and regular intervals in vec.x and plot results.
  #' find_motions(vec.x, "skip_window" = 1000, "window_length" = 3000, "find_thresh" = 300)

  # The 0-1 test for chaos calculated in mowing window.
  test01_res = test_chaos01_mw(data, window_length, skip_window, skip_test01, test01_thresh)

  # Find borders of regular motion.
  reg_borders <- find_reg_borders(test01_res)

  # Find borders of chaotic motion.
  chaos_borders <- find_chaotic_borders(test01_res)

  # Optimizing the boundaries of regular motions.
  reg_borders_final <- optimize_reg(find_thresh, test01_thresh, reg_borders, test01_res, data, skip_window, window_length,
      skip_test01)

  # Optimizing the boundaries of chaotic motions.
  chaos_borders_final <- optimize_chaos(find_thresh, test01_thresh, chaos_borders, test01_res, data, skip_window, window_length,
      skip_test01)

  # plot results
  plot_borders(data, window_length, skip_window, test01_res, chaos_borders_final, reg_borders_final)
  
  return(list("regular" = do.call(cbind, reg_borders_final), "chaotic" = do.call(cbind, chaos_borders_final)))
}
