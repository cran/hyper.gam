

#' @title Sign Adjustment
#' 
#' @param object ..
#' 
#' @param ... parameters of function [cor_xy.gam()]
#' 
#' @keywords internal
#' @name sign_adjust
#' @export
sign_adjust <- function(object, ...) UseMethod(generic = 'sign_adjust')


#' @rdname sign_adjust
#' @returns 
#' Function [sign_adjust.hyper_gam()] returns a \link[base]{numeric} \link[base]{vector}.
#' @export sign_adjust.hyper_gam
#' @export
sign_adjust.hyper_gam <- function(object, ...) {
  
  sgn <- object |>
    cor_xy.gam(...) |>
    sign()
  
  return(sgn * object$linear.predictors)
  
}



