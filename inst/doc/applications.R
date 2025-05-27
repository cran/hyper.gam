## -----------------------------------------------------------------------------
#| warning: false
#| eval: false
# remotes::install_github('tingtingzhan/groupedHyperframe')
# remotes::install_github('tingtingzhan/hyper.gam')


## -----------------------------------------------------------------------------
#| warning: false
#| eval: false
# utils::install.packages('groupedHyperframe')
# utils::install.packages('hyper.gam')


## -----------------------------------------------------------------------------
#| message: false
library(groupedHyperframe)
library(hyper.gam)
library(survival)


## -----------------------------------------------------------------------------
#| echo: false
op = par(no.readonly = TRUE)
#options(mc.cores = 1L) # for CRAN submission


## -----------------------------------------------------------------------------
data(Ki67, package = 'groupedHyperframe')
Ki67


## -----------------------------------------------------------------------------
#| message: false
Ki67q = Ki67 |>
  aggregate_quantile(by = ~ patientID, probs = seq.int(from = .01, to = .99, by = .01))


## -----------------------------------------------------------------------------
Ki67q |> head()


## -----------------------------------------------------------------------------
m0 = hyper_gam(PFS ~ logKi67.quantile, data = Ki67q)


## -----------------------------------------------------------------------------
m1 = hyper_gam(PFS ~ logKi67.quantile, data = Ki67q, nonlinear = TRUE)


## -----------------------------------------------------------------------------
#| eval: false
#| fig-width: 5
#| fig-height: 5
# m0 |> integrandSurface() # please interact with it on your local computer


## -----------------------------------------------------------------------------
#| eval: true
#| label: fig-integrandSurface
#| fig-width: 5
#| fig-height: 5
#| fig-align: left
#| fig-cap: 'Nonlinear integrand surface, integrand paths and their projections to the "floor" and "backwall"'
m1 |> integrandSurface(n = 101L)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-show: hide
m0 |> persp() # a static figure


## -----------------------------------------------------------------------------
#| warning: false
#| fig-show: hide
m0 |> contour() # a static figure


## -----------------------------------------------------------------------------
#| warning: false
#| fig-show: hide
m1 |> persp() # a static figure


## -----------------------------------------------------------------------------
#| warning: false
#| fig-show: hide
m1 |> contour() # a static figure


## -----------------------------------------------------------------------------
set.seed(16); id = Ki67q |> nrow() |> seq_len() |> caret::createDataPartition(p = .8)
Ki67q_0 = Ki67q[id[[1L]],] # training set
Ki67q_1 = Ki67q[-id[[1L]],] # test set


## -----------------------------------------------------------------------------
m1a = hyper_gam(PFS ~ logKi67.quantile, nonlinear = TRUE, data = Ki67q_0)


## -----------------------------------------------------------------------------
m1a$linear.predictors |> 
  head()


## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: "Optimistically biased!!"
Ki67q_0[,c('PFS', 'age', 'race')] |> 
  as.data.frame() |> # invokes spatstat.geom::as.data.frame.hyperframe()
  data.frame(nlQI = m1a$linear.predictors) |>
  coxph(formula = PFS ~ age + nlQI, data = _)


## -----------------------------------------------------------------------------
m1a |>
  predict(newdata = Ki67q_1) |> 
  head()


## -----------------------------------------------------------------------------
Ki67q_1[,c('PFS', 'age', 'race')] |> 
  as.data.frame() |> # invokes spatstat.geom::as.data.frame.hyperframe()
  data.frame(nlQI = predict(m1a, newdata = Ki67q_1)) |>
  coxph(formula = PFS ~ age + nlQI, data = _)


## -----------------------------------------------------------------------------
set.seed(145); QI = m0 |> 
  kfoldPredict.hyper_gam(k = 10L, mc.cores = 1L)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3.5
#| warning: false
#| eval: false
# boxplot(QI ~ attr(QI, 'fold'), xlab = 'Fold')


## -----------------------------------------------------------------------------
set.seed(145); nlQI = m1 |> kfoldPredict.hyper_gam(k = 10L, mc.cores = 1L)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3.5
#| warning: false
#| eval: false
# boxplot(nlQI ~ attr(nlQI, 'fold'), xlab = 'Fold')


## -----------------------------------------------------------------------------
#| warning: false
Ki67q |>
  cbind(spatstat.geom::hyperframe(QI = QI)) |> 
  as.data.frame() |>
  coxph(formula = PFS ~ QI) |> summary()


## -----------------------------------------------------------------------------
#| warning: false
Ki67q |>
  cbind(spatstat.geom::hyperframe(nlQI = nlQI)) |> 
  as.data.frame() |>
  coxph(formula = PFS ~ nlQI) |> 
  summary()

