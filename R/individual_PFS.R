#' Title
#'
#' @param data
#' @param selected_PFS
#'
#' @return
#' @export
#'
#' @examples
individual_PFS <- function(data,
                           #ratio=ratio, PFS1=PFS1, PFS2=PFS2, status=status,
                           selected_PFS){
  if(selected_PFS==1){
    km_PFS1 <- survfit(Surv(PFS1, rep(1, nrow(data))) ~ 1, data = data)
    x <- ggsurvplot(
      km_PFS1, data=data,
      legend="none",
      title= "Progression Free Survival 1",
      xlab = "PFS1",
      ylab = "Cumulative hazard",
      break.time.by = 3, palette = c("#224DA9"),
      xlim = c(0,24),
      ylim = c(0,4.5), ggtheme = theme_survminer(font.main = 18),
      conf.int = FALSE,
      fun = "cumhaz")

    Wei_PFS1 <- survreg(Surv(PFS1, rep(1, nrow(data))) ~ 1, dist='weibull', data=data)

    return(list(x, summary(Wei_PFS1)))
  }else{
    km_PFS2 <- survfit(Surv(PFS2, status) ~ 1, data = data)
    x <- ggsurvplot(
      km_PFS2, data=data,
      legend="none",
      title= "Progression Free Survival 2",
      xlab = "PFS2",
      ylab = "Cumulative hazard",
      break.time.by = 3, palette = c("#E5703D"),
      xlim = c(0,24),
      ylim = c(0,4.5), ggtheme = theme_survminer(font.main = 18),
      conf.int = FALSE,
      fun = "cumhaz")
    Wei_PFS2 <- survreg(Surv(PFS2, status) ~ 1, dist='weibull', data=data)
    return(list(x, summary(Wei_PFS2)))
  }

}
