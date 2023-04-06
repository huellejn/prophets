#' Plot cumulative hazard ratio
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param selected_PFS a vector of selected PFS values, by default c("PFS1", "PFS2")
#'
#' @return A plot of the cumulative hazard ratio.
#' @export
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival Surv
#'
#' @examples
#' data(input)
#' plot_cumHaz(data = data, selected_PFS = "PFS2")
plot_cumHaz <- function(data,  
                        selected_PFS = c("PFS1", "PFS2") 
                        ) {
  # Checks
  if(! all(selected_PFS %in% c("PFS1", "PFS2")) ) 
    stop("The selected PFS has to be PFS1, PFS2, or a vector of both values.")
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status"))

  # Melt data
  data_long <- melt_data(data)
  
  if(length(selected_PFS) == 1) {
    if(selected_PFS == "PFS1") {
      data_long <- data_long[data_long$ct.line == 1, ]
    }
    if(selected_PFS == "PFS2") {
      data_long <- data_long[data_long$ct.line == 2, ]
    }
  }
  
  km <- survminer::surv_fit(survival::Surv(data_long$PFS, data_long$status) ~ ct.line, data = data_long)
  
  p <- survminer::ggsurvplot(
    km, 
    data = data_long, 
    fun = "cumhaz",
    title = paste("Progression Free Survival", paste(unique(data_long$ct.line), collapse="+")),
    xlab = "Time (months)", ylab = "Cumulative hazard",
    legend = "none",
    conf.int = FALSE,
    ggtheme = theme_survminer(font.main = 18),
    palette = ifelse(selected_PFS == "PFS1", "#224DA9", ifelse(selected_PFS == "PFS2", "#E5703D", c("#224DA9", "#E5703D"))),
    break.time.by = 3,
    xlim = c(0, 24), 
    ylim = c(0, 4.5)
  ) 
  return(p)
  
}