#' Analysis of progression free survival ratio (PFSr) using the Kaplan-Meier method
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#' @param plot a logical value. Indication if a plot should be generated.
#'
#' @return Kaplan-Meier based statistics of PFSr and optionally a Kaplan-Meier plot
#' @export
#' @importFrom dplyr mutate select
#' @importFrom survival Surv survfit 
#' @importFrom gtsummary tbl_survfit
#' @importFrom survminer ggsurvplot
#' @importFrom broom glance
#' @importFrom ggplot2 geom_label geom_segment
#' @importFrom grid unit
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' kaplanMeier_PFSr(data = input, delta = 1.3)
kaplanMeier_PFSr <- function(data, 
                             delta = 1,
                             plot = FALSE) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status", "ratio"))
  
  obj_survfit <- survival::survfit(survival::Surv(ratio, data$status) ~ 1, data = data)
  
  res_summary <- broom::glance(obj_survfit) %>%
    dplyr::select(records, events, median, conf.low, conf.high)
  
  res_stats <- tbl_survfit(
    obj_survfit,
    times = delta)$meta_data$df_stats[[1]] %>% 
    dplyr::mutate(method = "Kaplan-Meier") %>% 
    dplyr::select(method, delta = time, estimate, conf.low, conf.high)
  
  res_list <- list(PFSr_summary = res_summary, PFSr_estimator = res_stats)
  
   if(plot == TRUE) {
     
      plot_km <- survminer::ggsurvplot(obj_survfit, data = data,
                            legend = "none",
                            title = "Survival function estimate of the PFSratio\nvia the Kaplan-Meier method",
                            xlab = "PFSratio",
                            ylab = "Probability",
                            break.time.by = 0.2, 
                            size = 0.7,
                            palette = "#224DA9",
                            xlim = c(0, 3),
                            ylim = c(0, 1),
                            ggtheme = theme_survminer(font.main = 18),
                            risk.table = FALSE)
      
      plot_km <- plot_km$plot + 
        ggplot2::geom_segment(aes(x = 0, y = res_stats$estimate, xend = delta, yend = res_stats$estimate), 
                     colour = "#E5703D", 
                     linetype = "dashed",
                     size = 0.3, 
                     alpha = 0.8) +
        ggplot2::geom_segment(aes(x = delta, y = 0, xend = delta, yend = res_stats$estimate), 
                     colour = "#E5703D", 
                     linetype = "dashed",
                     size = 0.3, 
                     alpha = 0.8) + 
        ggplot2::geom_label(label = paste("At ∂ = ", delta, ", S(", delta, ") :", round(res_stats$estimate,2), sep = ""), 
                   x = delta+1, 
                   y = res_stats$estimate + res_stats$estimate/2,
                   label.padding = grid::unit(0.25, "lines"), 
                   label.size = 0.35,
                   color = "black",
                   fill="white"
        )
      
      res_list[["kaplanMeier_plot"]] <- plot_km
      
      } 
  
  return(res_list)
}
