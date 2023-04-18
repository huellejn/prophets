#' Analysis of progression free survival ratio (PFSr) using the count-based method
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#' @param prob a numeric value. The probability
#'
#' @return Count-based statistics of PFSr analysis 
#' @export
#' @importFrom dplyr filter mutate select
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom broom tidy
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' count_PFSr(data = input)
count_PFSr <- function(data, 
                       delta = 1, 
                       prob = 0.5 ) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("status", "ratio"))

  data <- dplyr::filter(data, .data$status == 1 | (.data$status == 0 & .data$ratio >= delta))
  pos <- sum(data$ratio >= delta)
  neg <- sum(data$ratio < delta)
  res <- stats::binom.test(pos, (neg+pos), p = prob) %>%
    broom::tidy() %>%
    dplyr::mutate(
      method = "count-based",
      delta = delta
    ) %>%
    dplyr::select(method, delta, estimate, conf.low, conf.high)
  
  return(res)
}
