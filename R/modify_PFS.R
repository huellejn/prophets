#' Modified progression free survival ratio (mPFSr)
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.  
#' @param min_pfs2 a numeric value. The value 'min_pfs2' is the minimum PFS2 in month that can be considered satisfying irrespective of PFS1.
#'
#' @return dataframe with modified PFS and PFSr.
#' @export
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
#' @examples
#' data(input)
#' modify_PFS(data = input, delta = 1.3, min_pfs2 = 6)
modify_PFS <- function(data, delta = 1, min_pfs2 = 6) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2"))
  
  data <- dplyr::mutate(data, 
      PFS1  = pmax(.data$PFS1, 2.0),
      ratio = .data$PFS2/.data$PFS1,
      PFS2  = ifelse(.data$PFS2 > min_pfs2 & .data$ratio < delta, .data$PFS1 * delta + 0.25, .data$PFS2),
      ratio = .data$PFS2/.data$PFS1
    ) 
  return(data)
}