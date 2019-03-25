# plot results ------------------------------------------------------------------------------------------------------
plot_borders <- function(data, window_length, skip_window, test01_res, chaos_borders_final, reg_borders_final) {
  ## Plot the chaotic and regular motions along with the data and results of the 0-1 test for chaos computed
  ## using moving window.
  ##
  ## @param data Analyzed data.
  ## @param window_length Length of the window for in which the 0-1 test for chaos will be computed.
  ## @param skip_window Length of the skip of the window moving in the data.
  ## @param test01_res The results of the 0-1 test for chaos in each computed window.
  ## @param chaos_borders_final The found borders of chaotic motion.
  ## @param reg_borders_final The found borders of regular motion.
  ## @importFrom graphics lines
  ## @importFrom graphics plot
  ## @importFrom graphics points

  test01_res <- test01_res$test01_res
  plot(seq(1, length(data) - window_length, skip_window) + window_length/2, test01_res, col = "green", xlab = "",
      ylim = c(-0.1, 1))
  lines(data, col = "red")
  lines(seq(1, length(data) - window_length, skip_window) + window_length/2, test01_res, col = "green")
  if (length(chaos_borders_final[[1]]) > 0) {
    points(chaos_borders_final[[1]], matrix(1, length(chaos_borders_final[[1]]), 1), col = "blue", pch = 16)
    points(chaos_borders_final[[2]], matrix(1, length(chaos_borders_final[[2]]), 1), col = "blue", pch = 16)
  }
  if (length(reg_borders_final[[1]]) > 0) {
    points(reg_borders_final[[1]], matrix(0, length(reg_borders_final[[1]]), 1), col = "blue", pch = 16)
    points(reg_borders_final[[2]], matrix(0, length(reg_borders_final[[2]]), 1), col = "blue", pch = 16)
  }
}
