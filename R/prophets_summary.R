#' Generate a summary of all survival statistics
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @return a table with PFSr analysis results using four different methods
#' @export
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' prophets_summary(data = input, delta = 1.3)
prophets_summary <- function(data, 
                             #points = NULL,
                             # conf.int = FALSE,
                             n.boot = 2000,
                             delta = 1) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status", "ratio"))
  
  res_countPFSr <- count_PFSr(data = data, delta = delta)
  res_kaplanMeierPFSr <- kaplanMeier_PFSr(data = data, delta = delta, plot = FALSE)[["PFSr_estimator"]]
  res_parametricPFSr <- parametric_PFSr(data = data, delta = delta)
  res_midrankPFSr <- midrank_PFSr(data = data, delta = delta)
  res_kernelPFSr <- kernelKM_PFSr(data = data, delta = delta, conf.int = T, n.boot = n.boot)
  # Format the results from the kernelKM method as data.frame
 # if(all(c("low.points", "upp.points") %in% names(res_kernelPFSr))) {
    conf.low = res_kernelPFSr$low.points
    conf.high = res_kernelPFSr$upp.points
  # } else {
  #   conf.low = NA
  #   conf.high = NA
  #   }
  df_res_kernelPFSr <- data.frame(
    method = "kernelKM",
    delta = delta,
    estimate = res_kernelPFSr$surv.points,
    conf.low = conf.low,
    conf.high = conf.high
  )
  res <- rbind(res_countPFSr, res_kaplanMeierPFSr, res_parametricPFSr, res_midrankPFSr, df_res_kernelPFSr) 
  return(res)
}