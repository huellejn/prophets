#' Title
#'
#' @param data
#' @param modified
#' @param min_pfs2
#' @param id
#'
#' @return
#' @export
#'
#' @examples
weibull_survplot <- function(data,
                             modified=F, min_pfs2=6, id=id){

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
    mutate(status1 = 1, status2=status,
           cc = as.factor({{id}})
    ) %>%
    select(cc, PFS1, status1, PFS2, status2) %>%
    pivot_longer(cols = c("PFS1","PFS2","status1","status2"),
                 names_to = c(".value", "ct.line"),
                 names_pattern = "^([A-Za-z]+)(\\d+)")

  fKM <- survfit(Surv(PFS,status) ~ as.factor(ct.line), data=ana_long)
  sWei1 <- survreg(Surv(PFS1, rep(1, nrow(data))) ~ 1,dist='weibull',data=data)
  pred.line1 = predict(sWei1, type="quantile",p=seq(.01,.99,by=.01))[1,]
  sWei2 <- survreg(Surv(PFS2, status) ~ 1,dist='weibull',data=data)
  pred.line2 = predict(sWei2, type="quantile",p=seq(.01,.99,by=.01))[1,]
  df = data.frame(y=seq(.99,.01,by=-.01), line1=pred.line1, line2=pred.line2)
  df_long = gather(df, key= "line", value="time", -y)
  p = ggsurvplot(fKM, data = ana_long,
                 title="Weibull survival plot",
                 risk.table = FALSE, size=0.7,
                 legend.labs=c("PFS1" ,"PFS2"),
                 legend.title="",
                 #xlim = c(0,24),
                 palette = c("#224DA9", "#E5703D"),
                 surv.median.line = "hv",
                 xlab = "PFS (Months)", ggtheme = theme_survminer(font.main = 18),
                 break.time.by = 3)
  p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=line))
  p
}
