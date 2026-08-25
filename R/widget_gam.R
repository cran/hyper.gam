

#' @title Alternative \link[graphics]{persp}ective Plot for \link[mgcv]{gam} Model
#' 
#' @description
#' An interactive \link[graphics]{persp}ective plot for 
#' \link[mgcv]{gam} model, 
#' rendered as an \CRANpkg{htmlwidgets} 
#' using the package \CRANpkg{plotly}.
#' 
#' @param x a \link[mgcv]{gam} model
#' 
#' @param ... parameters of the function \link[mgcv]{vis.gam}
#' 
#' @param colorscale to be passed into the function \link[plotly]{add_surface}
#' 
#' @returns 
#' The function [widget_gam()] returns a pretty \CRANpkg{htmlwidgets} created by **R** package \CRANpkg{plotly}.
#' 
#' @note
#' 
#' The internal utility functions `persp_gam_int()` and `newd_gam_int()` are based on the function \link[mgcv]{vis.gam}.
#' 
#' @examples
#' library(mgcv)
#' colorscale = list(c(0, 1), c('white', 'lightgreen'))
#' 
#' # examples from ?mgcv::te
#' test1 = \(x,z,sx=0.3,sz=0.4) { 
#'  x = x*20
#'  (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
#'   0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
#' }
#' n = 500
#' x = runif(n)/20; z = runif(n);
#' xs = seq(0,1,length=30)/20; zs = seq(0,1,length=30)
#' pr = data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
#' truth = matrix(test1(pr$x,pr$z),30,30)
#' persp(xs,zs,truth); title("truth")
#' f = test1(x,z)
#' y = f + rnorm(n)*0.2
#' b1 = gam(y ~ s(x,z))
#' b2 = gam(y ~ te(x,z))
#' b3 = gam(y ~ ti(x) + ti(z) + ti(x,z))
#' b4 = gam(y ~ ti(x) + ti(x,z,mc=c(0,1))) ## note z constrained!
#' \donttest{
#' widget_gam(b1, colorscale = colorscale)
#' widget_gam(b2, colorscale = colorscale)
#' widget_gam(b3, colorscale = colorscale)
#' widget_gam(b4, colorscale = colorscale)
#' }
#' @importFrom plotly plot_ly add_surface layout
#' @export
widget_gam <- function(x, ..., colorscale) {
  
  p_ <- persp_gam_int(x, ...)
  
  ret <- plot_ly() |> 
    add_surface(
      x = p_$m1, y = p_$m2,
      z = t.default(p_$z), # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
      cmin = p_$min.z, cmax = p_$max.z, 
      contours = list(
        z = list(
          show = TRUE,
          start = p_$min.z, end = p_$max.z, size = (p_$max.z - p_$min.z)/21,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(z = TRUE)
        )
      ),
      colorscale = colorscale,
      showscale = FALSE
    )
  
  if (length(p_$lo.z)) {
    ret <- ret |>
      add_surface(
        x = p_$m1, y = p_$m2,
        z = t.default(p_$lo.z), # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
        cmin = p_$min.z, cmax = p_$max.z,
        opacity = .3, 
        colorscale = colorscale,
        showscale = FALSE
      )
  }
  
  if (length(p_$hi.z)) {
    ret <- ret |>
      add_surface(
        x = p_$m1, y = p_$m2,
        z = t.default(p_$hi.z), # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
        cmin = p_$min.z, cmax = p_$max.z,
        opacity = .3, 
        colorscale = colorscale,
        showscale = FALSE
      )
  }
  
  ret |>
    layout(scene = list(
      xaxis = list(title = p_$xlab), 
      yaxis = list(title = p_$ylab),
      zaxis = list(title = p_$zlab)
    ))
  
}


#' @importFrom mgcv predict.gam
persp_gam_int <- \(
    x, 
    view = NULL, 
    se = -1, 
    type = c('link', 'response'), 
    zlim = NULL, 
    lp = 1, 
    ...
) {
  
  # re-written from ?mgcv::vis.gam
  
  type <- match.arg(type)
  zlab <- switch(type, link = {
    paste("linear predictor")
  }, response = {
    type
  })

  tmp <- newd_gam_int(x = x, view = view, lp = lp, ...)
  newd <- tmp$newd
  ex.tf <- tmp$ex.tf
  m1 <- tmp$m1
  m2 <- tmp$m2

  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (is.matrix(z)) {
    lp <- min(ncol(z), max(1, round(lp)))
    z <- z[, lp]
    fv$fit <- fv$fit[, lp]
    fv$se.fit <- fv$se.fit[, lp]
  }
  if (length(ex.tf)) {
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }

  if (!is.null(zlim)) {
    if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
      stop("Something wrong with zlim")
    min.z <- zlim[1]
    max.z <- zlim[2]
  } else {
    min.z <- min(fv$fit, na.rm = TRUE)
    max.z <- max(fv$fit, na.rm = TRUE)
  }
  if (min.z == max.z) {
    min.z <- min.z - 1
    max.z <- max.z + 1
  }

  dm <- c(length(m1), length(m2))
  # see inside [newd_gam_int()]: `m1` index first, `m2` next
  
  z <- fv$fit |>
    array(dim = dm)
   
  ret <- list(
    m1 = m1, m2 = m2, z = z,
    min.z = min.z, max.z = max.z,
    newd = newd,
    xlab = view[1],
    ylab = view[2],
    zlab = zlab
  )
    
  if (se > 0) {
    
    lo.z <- fv$fit - fv$se.fit * se
    hi.z <- fv$fit + fv$se.fit * se
    
    ret$max.z <- max(hi.z, na.rm = TRUE)
    ret$min.z <- min(lo.z, na.rm = TRUE)
    
    ret$lo.z <- lo.z |>
      array(dim = dm)

    ret$hi.z <- hi.z |>
      array(dim = dm)

  }
    
  return(ret)
    
}


#' @importFrom mgcv exclude.too.far
newd_gam_int <- \(
  x, 
  view = NULL, 
  too.far = 0, 
  cond = list(), 
  n.grid = 501L, # default was `30L` in ?mgcv::vis.gam
  lp = 1, 
  ...
) {
  
  # re-written from ?mgcv::vis.gam
  
  vs <- x$var.summary
  
  v.names <- names(vs)
  
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    #for (i in 1:length(v.names)) {
    for (i in seq_along(v.names)) {
      
      ok <- TRUE
      
      if (is.matrix(vs[[i]])) {
        ok <- FALSE
      } else if (is.factor(vs[[i]])) {
        # tzh wants to deprecate
        if (length(levels(vs[[i]])) <= 1) 
          ok <- FALSE
      } else {
        if (length(unique(vs[[i]])) == 1) 
          ok <- FALSE
      }
      
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      
      if (k == 2) 
        break
    }
    
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
    
  } else {
    
    if (sum(view %in% v.names) != 2) 
      stop(gettextf("view variables must be one of %s", paste(v.names, collapse = ", ")))
    
    for (i in 1:2) {
      if (!inherits(vs[[view[i]]], c("numeric", "factor"))) {
        stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
      }
    }
  }
  
  ok <- TRUE
  
  for (i in 1:2) {
    if (is.factor(vs[[view[i]]])) {
      # tzh wants to deprecate
      if (length(levels(vs[[view[i]]])) <= 1) 
        ok <- FALSE
    } else {
      if (length(unique(vs[[view[i]]])) <= 1) 
        ok <- FALSE
    }
  }
  
  if (!ok) 
    stop(gettextf("View variables must contain more than one value. view = c(%s,%s).", view[1], view[2]))
  
  if (is.factor(vs[[view[1]]])) {
    return(invisible()) # tzh deprecated
    #m1 <- fac.seq(vs[[view[1]]], n.grid)
  } else {
    r1 <- range(vs[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  
  if (is.factor(vs[[view[2]]])) {
    return(invisible()) # tzh deprecated
    #m2 <- fac.seq(vs[[view[2]]], n.grid)
  } else {
    r2 <- range(vs[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(vs)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- vs[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(vs[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(vs[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  
  #if (is.factor(m1)) {
  #  # tzh deprecated
  #  m1 <- as.numeric(m1)
  #  m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  #}
  #if (is.factor(m2)) {
  #  # tzh deprecated
  #  m2 <- as.numeric(m2)
  #  m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  #}
  
  # exclude-too-far
  # was after ?mgcv::predict.gam inside ?mgcv::vis.gam
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
  } else ex.tf <- NULL
  
  return(list(
    newd = newd,
    ex.tf = ex.tf,
    m1 = m1,
    m2 = m2
  ))
  

}



# tzh does not agree (mathematically) with function `fac.seq()` inside ?mgcv::vis.gam








