#' Title
#'
#' @param data
#' @param id
#' @param plot
#' @param modified
#' @param min_pfs2
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
plot_PFSr <- function(data, id = cc, plot=T,
                      modified=F, min_pfs2=6,
                      delta=1){

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
  data <- data %>% #dplyr::select(PFS1, PFS2, status, ratio, cc) %>%
    mutate(cc = as.factor({{id}}),
           cc = fct_reorder(cc, ratio),
           status=as.factor(status)
    )
  myLoc <-  min(which(levels(data$cc)  %in%  data$cc[data$ratio>=delta]))
  #data <- data %>% pivot_longer(cols = c(PFS1, PFS2), names_to = "pfs", values_to = "PFS1/2")
  pfs_plot =  ggplot(data, aes(x=cc, y=PFS1 + PFS2)) +
    geom_segment(aes(x=cc, xend=cc, y=0, yend=PFS1, colour="PFS1"),
                 size=3, alpha=0.6) +
    geom_segment(aes(x=cc, xend=cc, y=PFS1, yend=PFS1+PFS2, colour="PFS2"),
                 size=3, alpha=0.5,show.legend = F) +
    scale_color_manual(values=c("PFS1"="#224DA9", "PFS2"="#E5703D"), name="PFS") +
    guides(colour = guide_legend(override.aes = list(shape = NA))) +
    geom_point(aes(shape=status, size=status,
                   fill=status), stroke=0.2, alpha=1, show.legend = T) +
    scale_size_manual(values=c(2.4,2.2)) +
    scale_shape_manual(values=c(23,21))+
    scale_fill_manual(values=c('#FF033E','grey80')) +
    geom_vline(aes(xintercept= myLoc), linetype=2, color="#FF033E") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          # panel.grid.major = element_line(color = "grey80",
          #                                 size = 0.2,
          #                                 linetype = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    xlab("Patients ranked by descending PFSratio") +
    ylab("PFS1 + PFS2 (Months)") +
    coord_flip() +
    ggtitle("Individual Progression Free Survival 1 & 2") +
    labs(caption=paste("Red dots denote censored data. Patients above the dotted line have PFS2/PFS1 >=", delta, sep = "")) +
    theme(plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0, size = 12))
  pfs1 <- glance(survfit(Surv(PFS1, rep(1, nrow(data)))~1, data=data)) %>%
    dplyr::select(records, events, median, conf.low, conf.high)
  tab_pfs1 <- pfs1 %>%
    gt(caption = "Summary of PFS1") %>%
    fmt_number(columns = 2:5,decimals = 2)
  pfs2 <- glance(survfit(Surv(PFS2, as.numeric(as.character(data$status)))~1, data=data))%>%
    dplyr::select(records, events, median, conf.low, conf.high)
  tab_pfs2 <- pfs2 %>%
    gt(caption = "Summary of PFS2") %>%
    fmt_number(columns = 2:5,decimals = 2)
  pfsr <- glance(survfit(Surv(ratio, as.numeric(as.character(data$status)))~1, data=data)) %>%
    dplyr::select(records, events, median, conf.low, conf.high)
  tab_pfsr <- pfsr %>%
    gt(caption = "Summary of PFSratio") %>%
    fmt_number(columns = 2:5,decimals = 2)
  if(plot==T){
    return(list(pfs_plot, `Summary of PFSratio` = pfsr, `Summary of PFS1`=pfs1, `Summary of PFS2`= pfs2))
  }else{
    return(list(tab_pfs1, tab_pfs2, tab_pfsr))
  }
}
