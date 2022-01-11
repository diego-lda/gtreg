#' Calibration Function
#'
#' @description This function calibrates coefficients.
#'
#' @param Y
#' @param X
#' @param qY
#' @param dgp
#' @param link
#'
#' @return
#' @export
#'
#' @examples
calibrate_coef <- function(Y, X, qY = NULL, dgp = "dr", link = "probit"){;

  if(dgp == "loc"){

    coef       <- lm(Y ~ X)$coef
    intercoef  <- c(coef[1],1);      # Calibrated coefs for intercept
    Xcoef	     <- c(coef[2],0);      # Calibrated coefs for Z
    eps        <- Y - (coef[1]+coef[2]*X)
    sigma      <- sd(eps)
  }

  if(dgp == "dr"){

    gridY       <- quantile(Y,qY)     # Grid of Y values for distribution regression estimator
    coefs       <- distrfun(Y = Y, X = X, gridY = gridY, link = link)$coefs0

    intercoef   <- coef( lm(coefs[1,] ~ gridY) );      # Calibrated coefs for intercept
    Xcoef	   	  <- coef( lm(coefs[2,] ~ gridY) );      # Calibrated coefs for X

    sde         <- NULL
  }


  # Store calibrated coefs
  bvec <- c( intercoef[1], Xcoef[1], intercoef[2], Xcoef[2] )

  ans <- list(bvec = bvec, sde = sde)

  return(ans)

}
