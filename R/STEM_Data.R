#' @keywords internal
#' @noRd

`STEM_Data` <-
function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

# STEM_Data components
  comp <- c("z", "coordinates", "covariates")

#Verify if all the required components are given by user
  compInd <- match(comp, names(x))
    if (any(is.na(compInd)))
        stop(paste("Component(s)", paste(comp[is.na(compInd)], collapse=", "),
                   "is (are) missing"))

#Controls over the dimensions of the components
   if (sum(!(apply(x$z,1,is.numeric)))>0)  stop("Component z must be numeric")
   if (sum(!(apply(x$coordinates,1,is.numeric)))>0)  stop("Component coordinates must be numeric")
   if (sum(!(apply(x$covariates,1,is.numeric)))>0)  stop("Component covariates must be numeric")
   r <- ncol(x$covariates)
   n <- nrow(x$z)
   d <- ncol(x$z)

   if(nrow(x$covariates) != n*d) stop("The number or row of covariates must be n*d")
### Missing values are allowed in the response only. The EM algorithm handles
### them as prescribed by Durbin and Koopman (2012, 2nd ed.), Sect. 4.10: the
### measurement equation is restricted to the observed rows at each time point.
### The covariates and the coordinates must be complete, because the design
### matrix enters the closed-form M-step updates directly and the coordinates
### enter the distance matrix; accommodating gaps there would require a
### stochastic E-step.
   if (any( c(is.na(x$coordinates), is.na(x$covariates))))
        stop("Missing values are not allowed in components coordinates, covariates")
   if (all(is.na(x$z)))
        stop("Component z is entirely missing")
   if (any(apply(x$z, 2, function(col) all(is.na(col)))))
        stop("Component z has one or more locations with no observed value at all")
   if(!(nrow(x$coordinates) == d && ncol(x$coordinates) == 2)) stop("The dimension of matrix coordinates must be d*2")

#Definition of the class STEM_Data
  x$r=r
  x$n=n
  x$d=d
  class(x) <- "STEM_Data"
   return(list(z=x$z,
		coordinates=x$coordinates,
		covariates=x$covariates,
		r=x$r,
		n=x$n,
		d=x$d))
}

