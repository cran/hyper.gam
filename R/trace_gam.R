

#' @title Trace on \link[mgcv]{gam} Model
#' 
#' @description
#' Trace on the surface of a \link[mgcv]{gam} model.
#' 
#' @param object a \link[mgcv]{gam} model
#' 
#' @param formula one-sided \link[stats]{formula}
#' 
#' @param newdata \link[base]{data.frame}, the test set
#' 
#' @param newid \link[base]{integer} scalar or \link[base]{vector},
#' row indices of `newdata` to be visualized. 
#' Default `1:3`, i.e., the first three test subjects.
#' 
#' @param ... additional parameters, to be passed into 
#' the function \link[plotly]{plot_ly}
#' 
#' @returns 
#' The function [trace_gam()] returns a pretty \CRANpkg{htmlwidgets} created by **R** package \CRANpkg{plotly}.
#' 
#' @note
#' The maintainer is not aware of any functionality of projection of arbitrary **curves** in package \CRANpkg{plotly}.
#' 
#' @examples
#' library(mgcv)
#' # ?s
#' # ?ti
#' 
#' @keywords internal
#' @importFrom mgcv predict.gam
#' @importFrom plotly plot_ly hide_colorbar
#' @importFrom stats setNames
#' @export
trace_gam <- function(
    object,
    formula, 
    newdata = data,
    newid = min(3L, nrow(newdata)) |> seq_len(), 
    ...
) {
  
  if (!inherits(object, what = 'gam')) stop('input needs to be `gam`')
  
  if (!is.call(formula) || formula[[1L]] != '~' || length(formula) != 2L) stop('`formula` must be one-sided formula')
  xnm <- formula[[2L]] # right-hand-side
  if (!is.symbol(xnm)) stop('Right-hand-side ', xnm |> deparse1() |> col_magenta(), ' must be a symbol')
  
  data <- object$data
  
  X <- data[[paste(xnm, 'y', sep = '.')]]
  x. <- as.double(colnames(X))
  nx <- length(x.)
  
  newX <- newdata[[paste(xnm, 'y', sep = '.')]]
  if (!is.matrix(newX)) stop('`newdata` does not contain a matrix column of functional predictor values')
  newx. <- newX |> colnames() |> as.double()
  if (!all.equal.numeric(newx., x.)) stop('grid of training and test data must be exactly the same')
  
  #l <- unique.default(data$L)
  l <- unique.default(data[[paste(xnm, 'L', sep = '.')]])
  if (length(l) != 1L) stop('wont happen')
  
  if (!length(newid)) stop('must have `newid`')
  
  if (!is.integer(newid) || anyNA(newid) || any(newid > nrow(newX))) stop('illegal `newid`')
  
  d <- data.frame(
    x = x.,
    y = newX[newid, , drop = FALSE] |> t.default() |> c(),
    id = rep(newid, each = nx),
    L = l
  )
  
  d_ <- d |>
    setNames(nm = c(
      paste(xnm, 'x', sep = '.'),
      paste(xnm, 'y', sep = '.'),
      'id',
      paste(xnm, 'L', sep = '.')
    ))
  d$z <- predict.gam(object, newdata = d_, se.fit = FALSE, type = 'link')
  
  plot_ly(
    data = d, 
    x = ~x, y = ~y, z = ~z, name = ~id, color = ~id,
    showlegend = FALSE, # need!
    type = 'scatter3d', 
    mode = 'markers',
    # mode = 'lines+markers', # should **not** plot a connected line!!
    ...
  ) |>
    hide_colorbar() # need!
  
}



