#' Test if PFS1 and PFS2 follow a Weibull distribution
#'
#' @param data
#' @param plot
#' @param modified
#' @param min_pfs2
#'
#' @return
#' @export
#'
#' @examples
test_individual_pfs <- function(data, plot=F,
                                modified=F, min_pfs2=6){

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
  Wei_PFS1 <- tidy(survreg(Surv(PFS1, rep(1, nrow(data))) ~ 1, dist='weibull', data=data))
  tabul_Wei_PFS1 <- Wei_PFS1 %>%
    gt(caption = "PFS1 Weibull model") %>%
    fmt_number(columns = 2:5,decimals = 2)

  Wei_PFS2 <- tidy(survreg(Surv(PFS2, status) ~ 1, dist='weibull', data=data))
  tabul_Wei_PFS2 <- Wei_PFS2 %>%
    gt(caption = "PFS2 Weibull model") %>%
    fmt_number(columns = 2:5,decimals = 2)

  if(plot==T){
    test <- data %>% pivot_longer(cols = c("PFS1", "PFS2"), names_to = "PFS", values_to = "val")
    try <- survfit(Surv(val,status)~PFS,data=test)
    weib_plot <- ggsurvplot(try,
                            data=test,
                            fun="cloglog", size=0.7,
                            title="Weibull regression diagnostic plot",
                            risk.table = FALSE,
                            ggtheme = theme_survminer(font.main = 18),
                            legend.labs=c("PFS1" ,"PFS2"),
                            legend.title="",
                            #xlim = c(0,24),
                            palette = c("#224DA9", "#E5703D"), #surv.median.line = "hv",
                            xlab = "PFS (Log Scale)", ylab="Log Cumulative Hazard",
                            break.time.by = 1)
    return(list(`PFS1 Weibull model`=Wei_PFS1, `PFS2 Weibull model`=Wei_PFS2, weib_plot))
  }else{
    return(list(tabul_Wei_PFS1, Wei_PFS1, tabul_Wei_PFS2, Wei_PFS2))
  }
}
