


#' @importFrom nlme getData
#' @export
getData.gam <- function(object) {
  
  if (!is.data.frame(data <- object$data)) stop('re-run input gam with `control = list(keepData = TRUE)`')
  
  # see code name in ?nlme:::getData.lme
  if (length(naAct <- object[['na.action']])) {
    if (inherits(naAct, what = 'omit')) {
      data <- data[-naAct, , drop = FALSE]
    } else if (inherits(naAct, what = 'exclude')) {
      # in ?nlme:::getData.lme
      # tzh does not understand this yet
      # do nothing!!!
    } else if (length(naAct_f <- eval(object$call$na.action))) {
      data <- naAct_f(data)
    } # else do nothing
  }
  
  return(data)
  
}

