## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: "Environment on author's computer"
#| label: author-env
Sys.info()[c('sysname', 'release', 'machine')]
R.version


## -----------------------------------------------------------------------------
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

