#' Title
#'
#' @param data
#' @param delta
#' @param plot
#' @param modified
#' @param min_pfs2
#' @param log_scale
#'
#' @return
#' @export
#'
#' @examples
test_pfs_correlation <- function(data, delta=1,  plot=F,
                                 modified=F, min_pfs2=6, log_scale=F){

  x <- c("PFS1", "PFS2", "status")
  if(all(x %in% colnames(data))==F)
    stop("Error: some required column is missing (PFS1, PFS2 or status)")

  if(modified==T){
    data = data %>%
      mutate(
        PFS1  = pmax(PFS1,2.0),
        ratio = round(PFS2/PFS1,3),
        PFS2  = ifelse(PFS2>min_pfs2 & ratio < delta, PFS1*delta+0.25, PFS2),
        ratio = round(PFS2/PFS1,3)
      )
  } else{
    data = data
  }
  ana_long <- data %>%
    mutate(status1 = 1,
           status2=status,
           idx=1:nrow(.)) %>%
    select(idx, PFS1, status1, PFS2, status2) %>%
    filter(!is.na(status2)) %>%
    pivot_longer(cols = c("PFS1","PFS2","status1","status2"),
                 names_to = c(".value", "ct.line"),
                 names_pattern = "^([A-Za-z]+)(\\d+)")
  gfm <- parfm(Surv(PFS, status) ~ ct.line, cluster = "idx", data = ana_long,
               dist = "weibull", frailty = "gamma")
  tau <- tau(gfm)
  if(plot==T){
    correl_plot <- data %>%
      mutate(pfsr_delta=ifelse(ratio > delta, "above", "below")) %>%
      ggplot(aes(x=PFS1, y =PFS2, color=pfsr_delta, fill=pfsr_delta, size=ratio)) +
      geom_point(alpha=0.7)+
      scale_color_manual(values=c( "#E5703D", "grey60")) +
      theme_bw()+
      theme(legend.position = "none") +
      xlab("PFS1 (Months)") + ylab("PFS2 (Months)") +
      ggtitle("Correlation between PFS1 and PFS2", subtitle = paste("Colored dots represent cases with PFSratio > ", delta)) +
      annotate(geom="text", x=mean(data$PFS1) + 1, y=max(data$PFS2), label= paste("Kendall's Tau: ", round(tau,3)))
    if(log_scale==T){
      correl_plot <- correl_plot +
        scale_x_continuous(trans='log10') +  scale_y_continuous(trans='log10')
      return(list(gfm, correl_plot))
    }else{
      return(list(gfm, correl_plot))
    }
  }else{
    gfm
  }
}
