#' Test if the PFS values follow a Weibull distribution
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param tidy a logical value indicating if the data should be returned in the tidy data format
#'
#' @return a table with Weibull survival statistics
#' @export
#' @importFrom survival survreg Surv
#' @importFrom broom tidy
#'
#' @examples
#' data(input)
#' test_weibull_distribution(data = input, tidy = TRUE)
test_weibull_distribution <- function(data, tidy = TRUE) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status"))
  
  wei1 <- survival::survreg(survival::Surv(data$PFS1, rep(1, times = nrow(data)) ) ~ 1, dist = 'weibull', data = data)
  wei2 <- survival::survreg(survival::Surv(data$PFS2, data$status) ~ 1, dist = 'weibull', data = data)
  
  if(tidy == TRUE) {
    wei1 <- broom::tidy(wei1)
    wei2 <- broom::tidy(wei2)
  }
  
  l_res <- list(
    PFS1 = wei1,
    PFS2 = wei2
  )
  
  return(
    l_res
  )
}