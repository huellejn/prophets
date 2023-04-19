#' Plot the correlation of PFS values
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#' @param log_scale a logical value indicating if the scale should be log transformed.
#'
#' @return A plot of the correlation between PFS1 and PFS2.
#' @export
#' @importFrom parfm tau
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw theme labs annotate scale_x_continuous scale_y_continuous
#' @importFrom dplyr mutate 
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' plot_correlation_PFS(data = input, delta = 1.3, log_scale = TRUE)
plot_correlation_PFS <- function(data, 
                                 delta = 1,
                                 log_scale = FALSE) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("status", "ratio"))
  
  glm <- correlate_PFS(data, delta = delta)
  
  tau <- parfm::tau(glm)
  
  p <- data %>% 
    dplyr::mutate(pfsr_delta = ifelse(ratio > {{delta}}, "above", "below")) %>% 
    ggplot2::ggplot(aes(x = PFS1, y = PFS2, color = pfsr_delta, fill = pfsr_delta, size = ratio)) + 
    geom_point(alpha = 0.7)+
    scale_color_manual(values = c("#E5703D", "grey60")) +
    theme_bw()+ 
    theme(legend.position = "none") +
    labs(
      title = "Correlation between PFS1 and PFS2", 
      subtitle = paste0("Colored dots represent cases with PFSratio >", delta),
      x = "PFS1 (Months)", 
      y = "PFS2 (Months)") +
    annotate(geom = "text", x = mean(data$PFS1) + 1, y = max(data$PFS2), label = paste("Kendall's Tau: ", round(tau,3) )) 
  
  if(log_scale == TRUE){
    p <- p + 
      scale_x_continuous(trans = 'log10') +  
      scale_y_continuous(trans = 'log10')
  }
  
  return(p)
}
