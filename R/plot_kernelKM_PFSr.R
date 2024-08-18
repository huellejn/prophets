#' Plot results from the kernel conditional Kaplan-Meier statistics
#'
#' @param res_kernelKM_PFSr results from the kernelKM_PFSr function, min requirement is a list with elements ratio and surv
#'
#' @return a ggplot object displaying a Kaplan-Meier curve for survival data calculated using the kernel conditional KM method
#' @export
#' @importFrom ggplot2 geom_step ggtitle xlab ylab scale_x_continuous geom_ribbon annotate geom_label
#' @importFrom survminer theme_survminer
#'
#' @examples
#' res <- kernelKM_PFSr(data = input, ratio.name = "ratio", status.name = "status", time0.name = "PFS1", conf.int = TRUE, points = 2)
#' plot_kernelKM_PFSr(res_kernelKM_PFSr = res)
plot_kernelKM_PFSr <- function(res_kernelKM_PFSr) {
  
  dat = data.frame(ratio = res_kernelKM_PFSr$ratio, surv = res_kernelKM_PFSr$surv)
  
  ## Generate the base plot
  p = ggplot(dat, aes(x = ratio, y = surv)) +
    geom_step(color = "#224DA9", linewidth = 0.7) +
    ggtitle("Survival function estimate of the PFSratio\nvia the kernel-based Kaplan-Meier method") +
    xlab("PFSratio") + ylab("Probability") +
    theme_survminer() +
    scale_x_continuous(breaks = seq(0, max(dat$ratio)+3, by = 1)) 
  
  ## Add confidence interval if conf.int is TRUE
  if( all(c("low", "upp") %in% names(res_kernelKM_PFSr)) ) {
    
    dat = data.frame(dat, low = res_kernelKM_PFSr$low, upp = res_kernelKM_PFSr$upp)
    
    p = p +
      geom_ribbon(data = dat, aes(ymin = low, ymax = upp), fill = "#224DA9", linetype = 2, alpha = 0.1)
  }
  
  ## Add points if points was not null
  if( all(c("points", "surv.points") %in% names(res_kernelKM_PFSr)) ) {
    
    val_points = res_kernelKM_PFSr$points
    val_surv.points = res_kernelKM_PFSr$surv.points
    
    p = p + 
      annotate("segment", x = 0, xend = val_points, y = val_surv.points, yend = val_surv.points,
               col =  "#E5703D", 
               alpha = .8,
               linetype = "dashed",
               linewidth = .3) +
      annotate("segment", x = val_points, xend = val_points, y = 0, yend = val_surv.points,
               col =  "#E5703D", 
               alpha = .8,
               linetype = "dashed",
               linewidth = .3) +
      geom_label(
        label = paste("At ∂ = ", val_points, ", S(", val_points, "): ", round(val_surv.points, 2), sep = ""),
        x = val_points + 0.5,
        y = val_surv.points + val_surv.points / 4,
        label.padding = unit(0.25, "lines"),
        label.size = 0.35,
        color = "black",
        fill = "white"
      )
    
  }
  return(p)
}
