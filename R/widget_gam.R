
# read ?mgcv::vis.gam carefully
# why ?mgcv::vis.gam does not require (an equivalent of) the parameter `formula`


#' @title Alternative of \link[mgcv]{vis.gam}
#' 
#' @description
#' An interactive \CRANpkg{htmlwidgets} of the 
#' \link[graphics]{persp}ective plot for 
#' \link[mgcv]{gam} model(s)
#' using the package \CRANpkg{plotly}.
#' 
#' @param ... one or more \link[mgcv]{gam} models
#' based on *a same data set*.
#' 
#' @param formula one-sided \link[stats]{formula}
#' 
#' @param newdata (for future expansion)
#' 
#' @param proj_xy \link[base]{logical} scalar, whether to show 
#' the projection to the \eqn{(x,y)}-plain, default value is `TRUE`
#' 
#' @param proj_xz \link[base]{logical} scalar, whether to show
#' the projection to the \eqn{(x,z)}-plain, default value is `TRUE`
#' 
#' @param n \link[base]{integer} scalar, fineness of visualization,
#' default value is `501L`. See parameter `n.grid` of the function \link[mgcv]{vis.gam}.
#' 
#' @param newid \link[base]{integer} scalar or \link[base]{vector},
#' row indices of `newdata` to be visualized. 
#' Default `1:2`, i.e., the first two test subjects.
#' Use `newid = NULL` to disable the visualization of `newdata`.
#' 
#' @param ylim \link[base]{length}-2 \link[base]{double} \link[base]{vector},
#' range on \eqn{q}-axis. Default is the range of \eqn{X} and \eqn{X^{\text{new}}} combined.
#' 
#' @param surface_col \link[base]{length}-2 \link[base]{character} \link[base]{vector},
#' color of the integrand surface(s), for lowest and highest surface values
#' 
#' @returns 
#' The function [widget_gam()] returns a pretty \CRANpkg{htmlwidgets} created by **R** package \CRANpkg{plotly}.
#' 
#' @note
#' The maintainer is not aware of any functionality of projection of arbitrary curves in package \CRANpkg{plotly}.
#' Currently, the projections are hard coded.
#' 
#' @examples
#' library(mgcv)
#' # ?gam
#' set.seed(2)
#' dat = gamSim(1, n = 400, dist = 'normal', scale = 2)
#' gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#' 
#' @keywords internal
#' @importFrom mgcv predict.gam
#' @importFrom plotly plot_ly add_paths add_surface
#' @importFrom stats setNames
#' @export
widget_gam <- function(
    ...,
    formula, 
    newdata = data,
    proj_xy = TRUE, 
    proj_xz = TRUE,
    n = 501L,
    newid = min(3L, nrow(newdata)) |> seq_len(), 
    ylim = range(X[is.finite(X)], newX[is.finite(newX)]), # removing NA, NaN, Inf
    surface_col = 
      # c('lightyellow', 'lightpink') # nice
      # c('beige', 'lightpink') # nice
      # c('white', 'deeppink') # not good!
      # c('white', 'magenta') # not good!
      c('white', 'lightgreen') # nice
    # c('white', 'darkgreen') # not good!
    # c('white', 'lightgoldenrod') # my R do not recognize
    # c('white', 'lightslateblue') # my R do not recognize
    # c('white', 'yellow') # nice
) {
  
  dots <- list(...)
  if (!all(vapply(dots, FUN = inherits, what = 'gam', FUN.VALUE = NA))) stop('all input needs to be `gam.matrix`')
  
  if (!is.call(formula) || formula[[1L]] != '~' || length(formula) != 2L) stop('`formula` must be one-sided formula')
  xnm <- formula[[2L]] # right-hand-side
  if (!is.symbol(xnm)) stop('Right-hand-side ', xnm |> deparse1() |> col_magenta(), ' must be a symbol')
  
  data_ <- dots |> lapply(FUN = \(i) i$data) |> unique()
  if (length(data_) > 1L) stop('data not same')
  data <- data_[[1L]]
  
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
  
  # plot!!
  # *surface* based on training model
  x_ <- seq.int(from = min(x.), to = max(x.), length.out = n)
  y_ <- seq.int(from = ylim[1L], to = ylim[2L], length.out = n)
  d_xy <- data.frame(
    expand.grid(x = x_, y_), # span `x_` first, then span `y_`
    L = l
  ) |> 
    setNames(nm = c(
      paste(xnm, 'x', sep = '.'),
      paste(xnm, 'y', sep = '.'),
      paste(xnm, 'L', sep = '.')
    ))
  
  zs <- mapply(FUN = \(x) { # (x = dots[[1L]])
    y0 <- predict.gam(x, newdata = d_xy, se.fit = FALSE, type = 'link')
    dim(y0) <- c(n, n)
    t.default(y0) # important!!!
    # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
  }, x = dots, SIMPLIFY = FALSE)
  
  zmin <- zs |> unlist() |> min()
  zmax <- zs |> unlist() |> max()
  
  p <- plot_ly()
  
  for (z_ in zs) {
    p <- p |> 
      add_surface(
        x = x_, y = y_,
        z = z_, cmin = zmin, cmax = zmax, 
        contours = list(
          z = list(
            show = TRUE,
            start = zmin, end = zmax, size = (zmax - zmin)/21,
            usecolormap = TRUE,
            highlightcolor = "#ff0000",
            project = list(z = TRUE)
          )
        ),
        colorscale = list(c(0, 1), surface_col), 
        showscale = FALSE
      )
  }
  
  if (!length(newid)) return(p)
  
  if (!is.integer(newid) || anyNA(newid) || any(newid > nrow(newX))) stop('illegal `newid`')
  
  d <- data.frame(
    x = x.,
    y = newX[newid, , drop = FALSE] |> t.default() |> c(),
    id = rep(newid, each = nx),
    L = l
  )
  if (proj_xy) {
    p <- p |> 
      add_paths(
        x = d$x, y = d$y, z = zmin, name = d$id, color = d$id, 
        showlegend = FALSE,
        line = list(width = 4)
      )
  } # projection on x-y plain
  
  d_ <- d |>
    setNames(nm = c(
      paste(xnm, 'x', sep = '.'),
      paste(xnm, 'y', sep = '.'),
      'id',
      paste(xnm, 'L', sep = '.')
    ))
  z_subj <- mapply(FUN = \(x) {
    predict.gam(x, newdata = d_, se.fit = FALSE, type = 'link')
  }, x = dots, SIMPLIFY = FALSE)
  
  if (proj_xz) {
    # projection on x-z plain, F(p, Q(p)) curve
    # this is only done if (length(dots) == 1L); otherwise too messy
    if (length(dots) == 1L) {
      for (i in seq_along(dots)) {
        p <- p |> 
          add_paths(
            x = d$x, y = ylim[2L], z = z_subj[[i]], name = d$id, color = d$id,
            showlegend = FALSE,
            line = list(width = 4)
          )
      }
    }
  } # projection on x-z plain
  
  for (i in seq_along(dots)) {
    p <- p |> 
      add_paths(
        x = d$x, y = d$y, z = z_subj[[i]], name = d$id, color = d$id,
        showlegend = FALSE,
        line = list(width = 4)
      )
  }
  
  return(p)
  
}











if (FALSE) { # learn projection with plotly
  
    
  # python
  # https://plotly.com/python/3d-line-plots/ # no mention of projection
  # https://community.plotly.com/t/how-to-plot-a-2d-graph-on-the-background-side-wall-of-a-3d-plot/72874/5
  
  # R
  # https://plotly.com/r/3d-line-plots/ # no mention of projection
  # https://stackoverflow.com/questions/53182432/3d-surface-with-a-2d-projection-using-r
  
  plot_ly(z = ~volcano) |> 
    add_surface(
      contours = list(
        z = list(
          show = TRUE,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(z=TRUE)
        ),
        y = list(
          show = TRUE,
          usecolormap = FALSE, # Projection without colormap
          highlightcolor = "#ff0000",
          project = list(y=TRUE)
        ),
        x = list(
          show = TRUE,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(x=TRUE)
        )
      )
    )
  
  # ?plotly::add_trace
  # explanation for parameter `...`
  # Arguments (i.e., attributes) passed along to the trace type. See schema() for a list of acceptable attributes for a given trace type (by going to traces -> type -> attributes).
  
  # ?plotly::schema
}
