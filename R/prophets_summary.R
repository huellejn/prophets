#' Title
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @return a table with PFSr analysis results using four different methods
#' @export
#'
#' @examples
#' data(input)
#' prophets_summary(data = data, delta = 1.3)
prophets_summary <- function(data, 
                             delta = 1) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status", "ratio"))
  
  res_countPFSr <- count_PFSr(data = data, delta = delta)
  res_kaplanMeierPFSr <- kaplanMeier_PFSr(data = data, delta = delta, plot = FALSE)[["PFSr_estimator"]]
  res_parametricPFSr <- parametric_PFSr(data = data, delta = delta)
  res_midrankPFSr <- midrank_PFSr(data = data, delta = delta)
  res <- rbind(res_countPFSr, res_kaplanMeierPFSr, res_parametricPFSr, res_midrankPFSr) 
  return(res)
}