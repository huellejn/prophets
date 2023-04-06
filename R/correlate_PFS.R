#' Calculate the correlation between PFS1 and PFS2
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @return Statistics of the correlation between PFS1 and PFS2
#' @export
#' @importFrom parfm parfm
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' data(input)
#' correlate_PFS(data = data, delta = 1.3)
correlate_PFS <- function(data, 
                          delta = 1) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status"))
  
  data_long <- melt_data(data)
  
  gfm <- parfm::parfm(survival::Surv(PFS, status) ~ ct.line, cluster = "idx", data = data_long,
               dist = "weibull", frailty = "gamma")
  
  return(gfm)
  
}