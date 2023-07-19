<!-- badges: start -->

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->

<img src="https://github.com/fedenichetti/myrepo/blob/main/imgfile1.png?raw=true" width="200" height="240" align="right"/>

# *PROPHETS*: PROgression-free survival ratio as primary endpoint in PHasE 2 TrialS in oncology

Authors: Federico Nichetti, Jennifer Hüllein

## Description

Progression Free Survival Ratio (PFSratio), as defined as the ratio between PFS on investigational treatment (PFS2) and PFS on the last prior therapy (PFS1), is a popular endpoint in precision oncology (PO) studies. This endpoint has some advantages compared to those usually adopted in phase II trials, like overall response rate or PFS, since it allows:

a)  reduction in heterogeneity, as each patient serves as his/her own control,

b)  design of single-arm trials with limited sample size, and

c)  provision of a clinically relevant estimate of benefit of a new treatment.

In recent years, a detailed exploration of PFSratio statistical properties has been performed, to correctly adopt it for the design and subsequent analysis of a trial, but the methods remained scattered and difficult to interpret for clinicians. Moreover, several PO trials explored PFSratio as a measure of benefit, but different methods have been adopted, making the results poorly comparable.

`prophets` is an R-Package which is developed to collect and implement different methods for the design and analysis of trials that use PFSratio as their primary endpoint. It provides a convenient wrapper around existing methods to calculate and plot PFSratio-based results.

The prophets pipeline is here summarized, from exploration of paired failure times to PFSratio analysis and PFSratio-based power/sample size calculation for trial design.

Detailed descriptions of each method can be found in the literature cited in the documentation.

## Installation

The developmental version can be installed using the `devtools` R-Package:

``` r
library(devtools)

devtools::install_github("https://github.com/huellejn/prophets")
```

## Bug Reports and Feature Requests

If you encounter any bugs or have any specific feature requests, please file an [Issue](https://github.com/huellejn/prophets/issues).

## *PROPHETS* pipeline

The following lines describe a standardised pipeline for PFSratio-based analysis, using the prophets functions for data analysis and plotting.

``` r
# load the package
library(prophets)
```

Within the package, an example dataset of 65 cancer patients treated with 1st and 2nd line chemotherapy is available. Here, we load these data, which include only 3 columns:

-   PFS1 and PFS2, namely progression free survival (in months) to treatment-1 and -2, as calculated as (`date of disease progression` - `date of treatment start`) / `365.25/12`.
-   status, i.e. censoring status to PFS2, reported as `1` if disease progression (PD) or death was observed on treatment 2, or `0` (censored) otherwise. Note that disease progression to treatment-1 is assumed as always occurred. It is up to the investigator whether to include cases that suspended treatment-1 for reasons other than PD.

``` r
data(input)
head(input)
```

PFSratio is calculated as PFS2/PFS1:

``` r
input$ratio <- input$PFS2 / input$PFS1
head(input)
```

Optionally, the investigator may be interested in using a modified PFSratio (mPFSratio), as described previously. Within `prophets`, the `modify_PFS` function is available to obtain the mPFSratio, where PFS1 \< 2 months is transformed to 2 months and PFS2 is transformed as *PFS1 \* 𝛿 + 0.25* when > an investigator set threshold (e.g. 6 months). 

``` r
data_modifiedPFS <- modify_PFS(input, delta = delta, min_pfs2 = 6)
head(input)
```

Delta (𝛿) is defined as the PFSratio cutoff that defines the investigational treatment as “effective”. For the following sections, we will set 𝛿 = 1.3. 

``` r
delta = 1.3
```

For the following sections, the standard PFSratio is adopted.

### Exploration of paired failure times 

Before PFSratio analysis, individual assessment of PFS1 and PFS2 is essential. In detail, it is of primary importance to assess the pattern of each PFS, the censoring rate of PFS2, and possibly the correlation between the two-paired PFS outcomes.

First, a `ggplot`-based function to generate a swimmer plot is implemented to visualize PFS1 and PFS2, highlighting censored cases with red dots. Cases are arranged according to PFSratio, with those with PFSratio > delta highlighted above the red dotted line.

``` r
swimmerplot_PFSr(
    input, 
    delta = delta)
```
<img src="https://github.com/fedenichetti/myrepo/blob/main/swimmer.png?raw=true">

### Correlation between PFS1 and PFS2

The calculation of the correlation index between PFS1 and PFS2 is crucial, as it affects the performance of PFSratio in the analysis phase and sample size calculation in study design. 
The `plot_correlation_PFS` function generates ad scatter plot with PFS1 and PFS2 on the x and y axis, respectively, and with the correlation index (Kendall's Tau) in the plot's right upper quadrant. 

``` r
plot_correlation_PFS(
  input, 
  delta = delta,
  log_scale = FALSE
)
```
<img src="https://github.com/fedenichetti/myrepo/blob/main/scatter.png?raw=true">

### Cumulative hazard ratio

Moreover, in order to select the most appropriate method for PFSratio-based analysis, it is necessary to verify if PFS1 and PFS2 follow a Weibull distribution. A way to do this is to plot the logarithm of the cumulative hazard, i.e. log[- log(S(t)] where log stands for the natural logarithm and S(t) is the Kaplan Meier survival estimate for PFS, against log(survival time). The Weibull assumption is correct if this plot gives two parallel and approximately straight lines for PFS1 and PFS2. The `plot_cumHaz` does the job:

``` r
plot_cumHaz(
  input,  
  selected_PFS = c("PFS1", "PFS2") 
)
```
<img src="https://github.com/fedenichetti/myrepo/blob/main/cumhaz.png?raw=true">

## Methods for PFSratio-based analysis

Whatever the choice of 𝛿, a study that is based on PFSratio as an efficacy endpoint is aimed at estimating the probability that this ratio is equal to or greater than 𝛿, i.e. P(PFSratio≥𝛿), or S~PFSratio~(𝛿) more compactly. S~PFSratio~(𝛿) (and its confidence interval, CI) thus represents the probability of having a ratio > 𝛿, which can be interpreted as the fraction of patients from a given cohort with a clinically relevant improvement of PFS2 relative to PFS1.
As described extensively in our work, there are multiple methods that can be used to estimate S~PFSratio~(𝛿). To date, 4 methods are implemented in `prophets`, namely the count-based (`res_countPFSr`), kaplan-meier (`res_kaplanMeierPFSr`), parametric (`res_parametricPFSr`) and midrank (`res_midrankPFSr`). Here, the kaplan-meier is used for the example, together with its respective plot. 

### Kaplan-Meier-based method

``` r
res_kaplanMeierPFSr <- kaplanMeier_PFSr(
  data = input,
  delta = delta,
  plot = TRUE
)

res_kaplanMeierPFSr
```
<img src="https://github.com/fedenichetti/myrepo/blob/main/km.png?raw=true">


Finally, a convenient function is available to perform the analysis with all available methods, summarizing results in a `tibble` format. 

### Summary of all methods

``` r
res_summary <- prophets_summary(
  data = input,
  delta = delta
)

res_summary
```
# Citation

The main paper associated with this R-Package is:
Nichetti F.,  Hüllein J., et al. Benchmarking Progression-Free Survival Ratio as primary endpoint in precision oncology clinical trials. ...


