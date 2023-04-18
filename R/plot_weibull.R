#' Plot Weibull distribution for PFS values
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#'
#' @return A plot of the Weibull distribution for PFS values
#' @export
#' @importFrom survival Surv
#' @importFrom survminer surv_fit ggsurvplot theme_survminer
#' @importFrom tidyr gather
#' @importFrom ggplot2 geom_line
#'
#' @examples
#' data(input)
#' plot_weibull(data = input)
plot_weibull <- function(data) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status"))
  
  # Melt data
  data_long <- melt_data(data)
  
  # Generate survfit object
  fKM <- survminer::surv_fit(survival::Surv(data_long$PFS, data_long$status) ~ as.factor(ct.line), data = data_long)
  
  p <- survminer::ggsurvplot(fKM, data = data_long,
                 title = "Weibull survival plot",
                 xlab = "Time in months", 
                 risk.table = FALSE, size = 0.7,
                 surv.median.line = "hv",
                 legend.labs = c("PFS1" ,"PFS2"), 
                 legend.title = "",
                 palette = c("#224DA9", "#E5703D"),
                 ggtheme = survminer::theme_survminer(font.main = 18),
                 break.time.by = 3)
  
  # Prediction lines
  
  ## Get prediction lines
  wei <- test_weibull_distribution(data, tidy = FALSE)
  pred.line1 = predict(wei[["PFS1"]], type = "quantile", p = seq(.01,.99, by = .01))[1,]
  pred.line2 = predict(wei[["PFS2"]], type = "quantile", p = seq(.01,.99, by = .01))[1,]
  
  ## Format data
  data_predLines = data.frame(y = seq(.99,.01, by = -.01), line1 = pred.line1, line2 = pred.line2) %>%
    tidyr::gather(key = "line", value = "time", -y)
  
  ## Add prediction lines to the plot
  p$plot = p$plot + 
    ggplot2::geom_line(data = data_predLines, aes(x = time, y = y, group = line))
  
  return(p)
  
}