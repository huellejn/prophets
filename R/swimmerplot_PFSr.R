#' Swimmer plot for individual PFS values
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @importFrom dplyr arrange mutate
#' @importFrom ggplot2 ggplot geom_segment geom_point theme theme_bw scale_color_manual scale_size_manual scale_shape_manual geom_vline labs element_text coord_flip guide_legend element_blank
#' @return A plot of individual PFS1 and PFS2 values, sorted by the PFSr.
#' @export
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' swimmerplot_PFSr(data = input, delta = 1.3)
swimmerplot_PFSr <- function(data, 
                      delta = 1) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2", "status", "ratio"))
  
  # sort by ratio and add a rank to determine the order in the plot
  data <- data %>%
    dplyr::arrange(ratio, desc = FALSE) %>%
    dplyr::mutate(
      rank = 1:nrow(.),
      status = as.logical(status)
    )
  
  cut_off_delta <- min(data$rank[data$ratio >= delta])
  
  p <- ggplot2::ggplot(data, aes(x = rank, y = (PFS1 + PFS2) )) +
    ggplot2::geom_segment(aes(x = rank, xend = rank, y = 0, yend = PFS1, colour = "PFS1"),
                 linewidth = 3, alpha = 0.6) +
    ggplot2::geom_segment(aes(x = rank, xend = rank, y = PFS1, yend = (PFS1 + PFS2), colour = "PFS2"),
                 linewidth = 3, alpha = 0.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("PFS1" = "#224DA9", "PFS2" = "#E5703D"), name = "PFS") + 
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = NA))) + 
    ggplot2::geom_point(aes(shape = status, size = status, fill = status), stroke = 0.2, alpha = 1, show.legend = T) + 
    ggplot2::scale_size_manual(values = c(2.4,2.2)) + 
    ggplot2::scale_shape_manual(values = c(23,21)) +
    ggplot2::scale_fill_manual(values = c('#FF033E','grey80')) +
    ggplot2::geom_vline(aes(xintercept = cut_off_delta), linetype = 2, color = "#FF033E") +
    ggplot2::labs(title = "Individual Progression Free Survival 1 & 2",
         x = "Patients ranked by descending PFSratio",
         y = "Time in months",
         caption = paste0("Red dots denote censored data.\nPatients above the dotted line have a PFSr >=", delta, ".")) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(),
      plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0, size = 12)
    ) +
    ggplot2::coord_flip() 
  
  return(p)
  
}