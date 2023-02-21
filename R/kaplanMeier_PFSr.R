#' Analysis of PFSr using the Kaplan-Meier method
#'
#' @param data
#' @param tidy
#' @param modified
#' @param min_pfs2
#' @param plot
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
kaplanMeier_PFSr <- function(data, tidy=T,
                             modified=FALSE, min_pfs2=6,
                             plot=FALSE,
                             delta=1 #, tidy=FALSE
){

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
      )} else{
        data = data
      }
  pfsr <- glance(survfit(Surv(ratio, as.numeric(as.character(data$status)))~1, data=data)) %>%
    dplyr::select(records, events, median, conf.low, conf.high)
  km_obJ <- survfit(Surv(ratio, status) ~ 1, data=data)
  aa <- tbl_survfit(
    km_obJ,
    times = delta)$meta_data$df_stats[[1]] %>%
    dplyr::select(time, estimate, conf.low, conf.high) %>%
    rename(delta="time") %>%
    mutate(Method = "Kaplan Meier") %>%
    relocate(Method, .before = delta)

  aa_tab <- aa %>%
    gt(caption = "PFSr estimator") %>%
    fmt_number(columns = 2:5,decimals = 2)
  pfsr_tab <- pfsr %>%
    gt(caption = "Summary of PFSratio") %>%
    fmt_number(columns = 2:5,decimals = 2)

  if(tidy==T){
    if(plot==T){
      plot_pfs <- ggsurvplot(km_obJ, data=data,
                             legend="none",
                             title= paste("Survival function estimate of the ",
                                          if(modified==T){"modified "}else{""}, "PFSratio\nvia the Kaplan-Meier method", sep = ""),
                             xlab = "PFSratio",
                             ylab = "Probability",
                             break.time.by = 0.2, size=0.7,
                             palette = c("#224DA9"),
                             xlim=c(0, 3),
                             ylim=c(0.0, 1),
                             #surv.median.line = "hv",
                             ggtheme = theme_survminer(font.main = 18),
                             risk.table=FALSE)
      plot_pfs <- plot_pfs$plot +
        geom_segment(aes(x = 0, y = aa$estimate, xend = delta, yend = aa$estimate),
                     colour= "#E5703D", linetype="dashed",
                     size=0.3,
                     alpha=0.8) +
        geom_segment(aes(x = delta, y = 0, xend = delta, yend = aa$estimate),
                     colour= "#E5703D", linetype="dashed",
                     size=0.3,
                     alpha=0.8) +
        geom_label(label=paste("At ∂ = ", delta, ", S(", delta, ") :", round(aa$estimate,2), sep = ""),
                   x=delta+1, #nudge_y = 2.25, nudge_x = 2.25,
                   y=aa$estimate+aa$estimate/2,
                   label.padding = unit(0.25, "lines"),
                   label.size = 0.35,
                   color = "black",
                   fill="white"
        )
      # geom_vline(xintercept = delta, linetype="dashed", color = "black") +
      # geom_hline(yintercept=aa$estimate, linetype="dashed", color = "black")

      return(list(plot_pfs, `Summary of PFSratio` = pfsr, `PFSratio estimator`= aa))
    }else{
      return(list(`Summary of PFSratio` = pfsr, `PFSratio estimator`= aa))
    }
  }
  else{
    return(list(pfsr_tab, aa_tab))
  }
}

