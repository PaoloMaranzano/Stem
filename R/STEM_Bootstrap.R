#' Parametric bootstrap
#'
#'
#' @description This functions performs the spatio-temporal parametric bootstrap for computing the parameter standard errors
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param B number of bootstrap iterations.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param precision optional; when given, it replaces \code{em_tol_par} of \code{control} in every refit (see \code{\link{STEM_Estimation}}). Default \code{NULL}.
#' @param control the computational settings of the refits, an object returned by \code{\link{STEM_control}} or a list of some of its settings. Default \code{NULL}, the defaults.
#' @param regularization a non-negative number added to the diagonal of the matrices inverted by the EM algorithm of every refit; see \code{\link{STEM_Estimation}}. Default is 0.
#' @param verbose logical. If TRUE, the progress of each bootstrap iteration is reported through message(). Default is FALSE.
#'
#'
#' @return The function returns a list of elements called \dQuote{boot.output}. \bold{Each} element of the list is an object of class \dQuote{STEM_Model} and so
#' it is composed by the following elements:
#' \itemize{
#' \item{skeleton} a list with components \code{phi}, \code{p}, \code{A} where the \code{phi} vector is given by \cr the \code{StemModel$estimates$phi.hat} vector.
#' \item{data} a list with components \code{z}, \code{coordinates}, \code{covariates}, \code{r}, \code{n} and \code{d} where the \code{z} matrix is
#'  given by the simulated data matrix.
#' \item{estimates} A list of four objects: \code{phi.hat}, \code{y.smoothed}, \code{loglik}, \code{convergence.par} as the output of the \code{\link{STEM_Estimation}} function.
#'}
#'
#'
#' @details The spatio-temporal bootstrap is used for parameter uncertainty assessment. The resampling scheme is based on the estimated model: each bootstrap sample is drawn
#' directly from the Gaussian distributions which define the model with the parameter vector replaced by the corresponding ML estimates. For each of  the \eqn{B} bootstrap
#' samples, the ML estimates are computed (using the procedure of \code{\link{STEM_Estimation}} function). Then the \eqn{B} bootstrap replications are returned in a list.
#'
#' @author Michela Cameletti  < michela.cameletti@unibg.it >
#'
#' @references Amisigo, B.A., Van De Giesen, N.C. (2005) \emph{Using a spatio-temporal dynamic state-space model with the EM algorithm to patch gaps in daily riverflow series}. Hydrology and Earth System Sciences 9, 209--224.
#' Fasso', A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring using heterogeneous networks}. Environmetrics 18, 245--264.
#' Fasso', A., Cameletti, M. (2007) \emph{A general spatio-temporal model for environmental data}. Tech.rep. n.27 \emph{Graspa} - The Italian Group of Environmental Statistics.
#' Fasso, A. and M. Cameletti (2010). \emph{A Unified Statistical Approach for Simulation, Modeling, Analysis and Mapping of Environmental Data}. SIMULATION 86(3): 139-153. <doi: 10.1177/0037549709102150>
#' Xu, K., Wikle, C.K. (2007) \emph{Estimation of parameterized spatio-temporal dynamic models}. Journal of Statistical Inference and Planning 137, 567--588.
#' Mc Lachlan, G.J., Krishnan, T. (1997) \emph{The EM Algorithm and Extensions}. Wiley, New York.
#' Shumway, R.H., Stoffer, D.S. (2006) \emph{Time Series Analysis and Its Applications: with R Examples}. Springer, New York.
#'
#' @examples
#' #load the data
#' data(pm10)
#'
#' #extract the data
#' coordinates <- pm10$coords
#' covariates <- pm10$covariates
#' z <- pm10$z
#'
#' #build the parameter list
#' #(the phi list is used for the algorithm starting values)
#' phi <- list(beta=matrix(c(3.65,0.046,-0.904),3,1),
#'             sigma2eps=0.1,
#'             sigma2omega=0.2,
#'             theta=0.01,
#'             G=matrix(0.77,1,1),
#'             Sigmaeta=matrix(0.3,1,1),
#'             m0=as.matrix(0),
#'             C0=as.matrix(1))
#'
#' A <- matrix(1,ncol(z),1)
#'
#' mod1 <- STEM_Model(z=z,covariates=covariates,
#'                    coordinates=coordinates,phi=phi,A=A)
#' class(mod1)
#'
#' \donttest{
#' #mod1 is given as output by the STEM_Model function
#' mod1.est <- STEM_Estimation(mod1)
#'
#' #it is computer intensive
#' mod1.boot <- STEM_Bootstrap(StemModel=mod1.est, B=3)
#' names(mod1.boot)
#' #the first element of the output list
#' names(mod1.boot$boot.output[[1]])
#'
#'
#' #check if there is no convergence for some bootstrap iteration
#' B <- length(mod1.boot$boot.output)
#' n.null <- sum (unlist (lapply(mod1.boot$boot.output,
#'                               function(x) length(x$StemModel)==1)) )
#' pos.null <- which((unlist (lapply(mod1.boot$boot.output,
#'                                   function(x) length(x$StemModel)==1)) ))
#' message("------------B non null =", B - n.null)
#' if(length(pos.null)>0) boot1.mod = mod1.boot[-pos.null]
#'
#' #put the bootstrap output in a matrix
#' npar <- length(unlist((mod1.boot$boot.output[[1]]$estimates$phi.hat)))-1
#' boot.estimates <- matrix(NA, nrow = (B - n.null), ncol = npar)
#'
#' for(b in 1:(B - n.null)) {
#'   phi.estimated <- mod1.boot$boot.output[[b]]$estimates$phi.hat
#'   boot.estimates[b,] <-  c(phi.estimated$beta,
#'                            phi.estimated$sigma2eps,
#'                            phi.estimated$sigma2omega,
#'                           phi.estimated$theta,
#'                           phi.estimated$G,
#'                           phi.estimated$Sigmaeta,
#'                           phi.estimated$m0)
#' }
#'
#' #compute the parameter standard errors
#' se <- sqrt(diag(var(na.omit(boot.estimates))))
#'
#' #create a summary table with Estimates, Standard Errors (SE) and T-statistics.
#' phi.hat <- mod1.est$estimates$phi.hat
#' MLE <- c(phi.hat$beta, phi.hat$sigma2eps, phi.hat$sigma2omega,
#'          phi.hat$theta, phi.hat$G, phi.hat$Sigmaeta,phi.hat$m0)
#' output1 <- cbind(MLE, se, MLE/se)
#' colnames(output1)<- c("Estimate", "SE", "T-stat.")
#' output1
#'
#' #compute the 95\% confidence intervals
#' IC <- matrix(NA,nrow=npar,ncol=2)
#' for(i in 1 : npar) {
#'   IC[i,] <- c(quantile(boot.estimates[,i],0.025),
#'               quantile(boot.estimates[,i],0.975))
#' }
#'
#'
#' #create a summary table with Estimates, Standard Errors (SE)
#' #and T-statistics and confidence intervals.
#' output2 <- cbind(output1,IC)
#' colnames(output2) <- c("Estimate", "SE", "T-stat.", "IC_inf", "IC_sup")
#' output2
#'
#' }
#'
#' @seealso See Also \code{\link{pm10}}, \code{\link{STEM_Model}} and \code{\link{STEM_Estimation}}
#'
#' @keywords models spatial
#'
#' @export

STEM_Bootstrap<-
  function(StemModel, B,distance='euclidean',precision=NULL,regularization=0, verbose = FALSE, control = NULL) {

    seed.list = list()
    for (i in 1:B) {
      seed.list[[i]] = as.integer(stats::runif(1,min=-1,max=1)*(10^8))
    }

    output = lapply(seed.list, STEM_Bootstrap.fn, StemModel = StemModel, seed.list.out = seed.list,distance=distance,precision=precision,regularization=regularization, verbose = verbose, control = control)
    return(list(boot.output=output))
  }



