#' Simulation of spatio-temporal data
#'
#'
#' @description The function \code{STEM_Simulation} simulates spatio-temporal data. of a hierarchical spatio-temporal model of class \dQuote{STEM_Model}. The estimates are obtained using Kalman filtering and EM algorithm.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#'
#'
#' @return The functions return a \eqn{n \times d}{n * d} matrix of data.
#'
#'
#' @details Note that the values contained in \code{StemModel$skeleton$phi} are used as the true values of the parameters.
#'
#'
#' @author Michela Cameletti  < michela.cameletti@unibg.it >
#'
#' @references Amisigo, B.A., Van De Giesen, N.C. (2005) \emph{Using a spatio-temporal dynamic state-space model with the EM algorithm to patch gaps in daily riverflow series}. Hydrology and Earth System Sciences 9, 209--224.
#' Fasso, A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring using heterogeneous networks}. Environmetrics 18, 245--264.
#' Fasso', A., Cameletti, M. (2007) \emph{A general spatio-temporal model for environmental data}. Tech.rep. n.27 \emph{Graspa} - The Italian Group of Environmental Statistics.
#' Fasso, A. and M. Cameletti (2010). \emph{A Unified Statistical Approach for Simulation, Modeling, Analysis and Mapping of Environmental Data}. SIMULATION 86(3): 139-153. <doi: 10.1177/0037549709102150>
#' Mc Lachlan, G.J., Krishnan, T. (1997) \emph{The EM Algorithm and Extensions}. Wiley, New York.
#' Shumway, R.H., Stoffer, D.S. (2006) \emph{Time Series Analysis and Its Applications: with R Examples}. Springer, New York.
#' Xu, K., Wikle, C.K. (2007) \emph{Estimation of parameterized spatio-temporal dynamic models}. Journal of Statistical Inference and Planning 137,  567--588.
#'
#' @examples
#'
#'
#' data(pm10)
#' names(pm10)
#'
#' #extract the data
#' coordinates <- pm10$coords
#' covariates <- pm10$covariates
#' z <- pm10$z
#'
#' #build the parameter list
#' phi <- list(beta=matrix(c(3.65,0.046,-0.904),3,1),
#'             sigma2eps=0.1,
#'             sigma2omega=0.2,
#'             theta=0.01,
#'             G=matrix(0.77,1,1),
#'             Sigmaeta=matrix(0.3,1,1),
#'             m0=as.matrix(0),
#'             C0=as.matrix(1))
#'
#' A <-matrix(1,ncol(z),1)
#'
#' mod1 <- STEM_Model(z=z,covariates=covariates,
#'                    coordinates=coordinates,phi=phi,A=A)
#'
#' class(mod1)
#'
#' simulateddata = STEM_Simulation(mod1)
#'
#'
#' @seealso See Also \code{\link{STEM_Model}} and \code{\link{pm10}}
#'
#' @keywords models spatial
#'
#'
#' @export


STEM_Simulation <-
function (StemModel,distance='euclidean'){

  n=StemModel$data$n
  d=StemModel$data$d
  r=StemModel$data$r
  p=StemModel$skeleton$p

  z = matrix(NA,nrow = d,ncol = n)
  y = matrix(NA,nrow = p,ncol = n)

  phi 			= StemModel$skeleton$phi
  phi$beta		= matrix(phi$beta,r,1)
  phi$logb 		= log(phi$sigma2eps/phi$sigma2omega)
  phi$logtheta	= log(phi$theta)

  if(distance=='euclidean'){dist     = as.matrix(stats::dist(StemModel$data$coordinates,diag=TRUE))} #distance matrix
  if(distance=='geo'){dist     = as.matrix(geodist::geodist(StemModel$data$coordinates,measure='geodesic'))}
  covariates 	= StemModel$data$covariates
  covariates  = changedimension_covariates(covariates,d=d,r=r,n=n)

  ####################
  ###Matrix definitions
  ####################
  Fmat = StemModel$skeleton$A
  Gmat = phi$G

  #cov.spaz = phi_real$sigma2omega * exp(-phi_real$theta * dist)
  #Vmat     = diag(phi_real$sigma2eps,d) + cov.spaz

  Vmat=phi$sigma2omega * Sigmastar.exp(d=d,logb=phi$logb,logtheta=phi$logtheta,dist=dist)

  Wmat = phi$Sigmaeta

  m0   = phi$m0
  C0   = phi$C0

  #####################
  y0 = matrix(MASS::mvrnorm(n=1, mu=m0, Sigma=C0),nrow=p, ncol=1)

  #t=1
  y[,1] = Gmat %*% y0 + MASS::mvrnorm(n=1, mu=matrix(0,nrow=p,ncol=1), Sigma=Wmat)
  z[,1] = covariates[,,1] %*% phi$beta + Fmat %*% y[,1] + MASS::mvrnorm(n=1, mu=matrix(0,nrow=d,ncol=1),Sigma = Vmat)

  for (t in 2:n) {
      y[,t] = Gmat %*% y[,t-1] + MASS::mvrnorm(n=1, mu=matrix(0,nrow=p,ncol=1), Sigma=Wmat)
      z[,t] = covariates[,,t] %*% phi$beta + Fmat %*% y[,t]   + MASS::mvrnorm(n=1, mu=matrix(0,nrow=d,ncol=1),Sigma=Vmat)
  }


  rownames(z) = rownames(StemModel$data$coordinates)
  return(z=t(z))
}

