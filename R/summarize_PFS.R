#' Summarise statistics of the PFS values
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) ) and 'ratio' (numeric, ratio PFS2/PFS1)
#'
#' @return A table with summary statistics on the PFS values.
#' @export
#' @importFrom survival Surv survfit
#' @importFrom dplyr select
#'
#' @examples
#' data(input)
#' summarize_PFS(data = data)
summarize_PFS <- function(data) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status", "ratio"))
  
  survPFS1 <- Surv(data$PFS1, rep(1, times = nrow(data)) )
  survPFS2 <- Surv(data$PFS1, data$status)
  survPFSr <- Surv(data$ratio, data$status)
  
  summary_PFS1 <- glance(survfit(survPFS1 ~ 1)) 
  summary_PFS2 <- glance(survfit(survPFS2 ~ 1)) 
  summary_PFSr <- glance(survfit(survPFSr ~ 1))
  
  summary <- rbind(
    data.frame(variable = "PFS1", summary_PFS1),
    data.frame(variable = "PFS2", summary_PFS2),
    data.frame(variable = "PFSr", summary_PFSr) 
  ) 
  summary <- summary[ , c("variable", "records", "events", "median", "conf.low", "conf.high")]
  
  return(summary)
  
}