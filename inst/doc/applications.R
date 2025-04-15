## -----------------------------------------------------------------------------
#| warning: false
#| eval: false
# devtools::install_github('tingtingzhan/hyper.gam')


## -----------------------------------------------------------------------------
#| warning: false
#| eval: false
# utils::install.packages('hyper.gam')


## -----------------------------------------------------------------------------
#| message: false
library(hyper.gam)
library(survival)


## -----------------------------------------------------------------------------
#| echo: false
library(knitr) # for tables in this vignette
op = par(no.readonly = TRUE)
#options(mc.cores = 1L) # for CRAN submission


## ----echo = FALSE, results = 'asis'-------------------------------------------
c(
  '', 'Forward pipe operator', '`?base::pipeOp` introduced in `R` 4.1.0', 
  '`attr`', 'Attributes', '`base::attr`; `base::attributes`',
  '`contour`', 'Contours', '`graphics::contour`; `hyper.gam::contour.hyper_gam`',
  '`coxph`', 'Cox model', '`survival::coxph`',
  '`gam`', 'Generalized additive models', '`mgcv::gam`',
  '`groupedHyperframe`, `hypercolumn`', '(Hyper column of) (grouped) hyper data frame', ' `groupedHyperframe::as.groupedHyperframe`; `spatstat.geom::hyperframe`', 
  '`htmlwidget`', 'HTML Widgets', '`` ?htmlwidgets::`htmlwidgets-package` ``; `plotly::plotly`',
  '`persp`', 'Perspective plot', '`graphics::persp`; `hyper.gam::persp.hyper_gam`',
  '`PFS`', 'Progression/recurrence free survival', '<https://en.wikipedia.org/wiki/Progression-free_survival>',
  '`quantile`', 'Quantile', '`stats::quantile`',
  '`S3`, `generic`, `methods`', '`S3` object oriented system',  '`base::UseMethod`; `utils::methods`; `utils::getS3method`; <https://adv-r.hadley.nz/s3.html>',
  '`Surv`', 'Survival object', '`survival::Surv`'
) |>
  matrix(nrow = 3L, dimnames = list(c('Term / Abbreviation', 'Description', 'Reference'), NULL)) |>
  t.default() |>
  as.data.frame.matrix() |> 
  kable()


## -----------------------------------------------------------------------------
data(Ki67, package = 'groupedHyperframe')
Ki67


## -----------------------------------------------------------------------------
#| message: false
s = Ki67 |>
  aggregate_quantile(by = ~ patientID, probs = seq.int(from = .01, to = .99, by = .01))
s |> head()


## -----------------------------------------------------------------------------
m0 = hyper_gam(PFS ~ logKi67.quantile, data = s)


## -----------------------------------------------------------------------------
#| eval: false
#| fig-width: 5
#| fig-height: 5
# integrandSurface(m0)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| warning: false
par(mar = c(2, 2, 0, 0))
persp(m0)
par(op)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| warning: false
par(mar = c(4, 5, 1, 0))
contour(m0)
par(op)


## -----------------------------------------------------------------------------
set.seed(145); QI = m0 |> 
  kfoldPredict.hyper_gam(k = 10L, mc.cores = 1L)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3.5
#| warning: false
par(mar = c(4, 5, 1, 0))
boxplot(QI ~ attr(QI, 'fold'), xlab = 'Fold')
par(op)


## -----------------------------------------------------------------------------
suppressWarnings(sQI <- cbind(s, spatstat.geom::hyperframe(QI = QI)) |> 
                   as.data.frame())
coxph(PFS ~ QI, data = sQI) |> summary()


## -----------------------------------------------------------------------------
m1 = hyper_gam(PFS ~ logKi67.quantile, data = s, nonlinear = TRUE)


## ----eval=FALSE,fig.width=5, fig.height=5-------------------------------------
# integrandSurface(m1)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| warning: false
par(mar = c(2, 2, 0, 0))
persp(m1)
par(op)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| warning: false
par(mar = c(4, 5, 1, 0))
contour(m1)
par(op)


## -----------------------------------------------------------------------------
set.seed(145); nlQI = m1 |> kfoldPredict.hyper_gam(k = 10L, mc.cores = 1L)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3.5
#| warning: false
par(mar = c(4, 5, 1, 0))
boxplot(nlQI ~ attr(nlQI, 'fold'), xlab = 'Fold')
par(op)


## -----------------------------------------------------------------------------
suppressWarnings(s_nlQI <- cbind(s, spatstat.geom::hyperframe(nlQI = nlQI)) |> 
                   as.data.frame())
coxph(PFS ~ nlQI, data = s_nlQI) |> summary()

