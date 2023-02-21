#' Title
#'
#' @param data
#' @param del
#' @param modified
#' @param min_pfs2
#'
#' @return
#' @export
#'
#' @examples
prophets_summary <- function(data, del, modified=F, min_pfs2=6){
  v <- count_PFSr(data = data, delta = del, modified = modified, min_pfs2=min_pfs2, tidy = T)
  w <- kaplanMeier_PFSr(data = data, plot = F, tidy=T, delta = del, modified = modified, min_pfs2=min_pfs2)
  w <- as_tibble(w$`PFSratio estimator`)
  x <- parametric_PFSr(data = data, modified = modified, tidy=T, min_pfs2=min_pfs2, delta = del)
  y <- midrank_PFSr(data = data, modified = modified, tidy=T, min_pfs2=min_pfs2, delta = del)
  z <- bind_rows(v,w,x,y) %>%
    as.data.frame()
  tabulated_data <- z %>%
    gt(caption = "Summary of PFSratio methods results") %>%
    fmt_number(columns = 2:5, decimals = 2)
  return(list(z,tabulated_data))
}
