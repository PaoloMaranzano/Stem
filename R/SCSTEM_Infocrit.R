#' optimal number of cluster selection for the STEM clustered regression
#'
#'
#' @description The function \code{SCSTEM_Infocrit} computes the likelihood-based information criteria for optimal number of cluster selection
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param mink minimum number of cluster
#' @param maxk maximum number of cluster
#' @param precision  a small positive number used for the EM algorithm convergence in every step of the clustered regression. Default is equal to 0.1. See \code{DETAILS} below.
#' @param precision_full_dataset a small positive number used for the EM algorithm convergence at the inizialization step of the EM algorithm. Default is equal to 0.1. See \code{DETAILS} below.
#' @param max_iter maximum number of iterations for the spatio temporal clustered regression algorithm.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param regularization a small positive number to be added to the digonal of the matrices matrices that need to be inverted . Default is set to 0.01
#' @param crs Integer value. Coordinate reference sySTEM_ something suitable as input to st_crs.command from the sf package (see its documentation for details). Default is set to 4326.
#' @param knn Integer value. The number of nearest neighbour to be taken into account for te spatial penalty, Default is set to 5
#' @param init_method Character. Must be one of: 'AMKM' or 'K-means'. If init_method='AMKM', the Adjacent Matrix K-Means clustering is performed. If method='K-means', K-means clustering is performed.
#' @param nugget_var Logical. If FALSE it returns the Sugasawa clustered regression. If TRUE the STEM spatio temporal clustered regression is performed. Default is set to TRUE
#' @param phi_penalty a small positive number. It is the spatial penalty weight for the clustered regression. Default is set to 1.
#' @param share2conv a small positive number. It is the minimum percentage of observations that must change clusters for the algorithm not to converge. Default is set equal to 0.05
#'
#'
#'
#' @return The function returns a table given by:
#' \itemize{
#' \item \code{BIC} Bayesian IC
#' \item \code{AIC} Akaike's IC
#' \item \code{HQC} Hannan–Quinn IC
#' }
#'
#'
#' @details This function allow an optimal selection for the number k of cluster for the SCSTEM_Estim method.
#'
#'
#' @author Francesco Caccia  < francesco.caccia2000@gmail.com >
#'
#' @references Amisigo, B.A., Van De Giesen, N.C. (2005) \emph{Using a spatio-temporal dynamic state-space model with the EM algorithm to patch gaps in daily riverflow series}. Hydrology and Earth System Sciences 9, 209--224.
#' Fasso, A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring using heterogeneous networks}. Environmetrics 18, 245--264.
#' Fasso', A., Cameletti, M. (2007) \emph{A general spatio-temporal model for environmental data}. Tech.rep. n.27 \emph{Graspa} - The Italian Group of Environmental Statistics.
#' Fassò, A. and M. Cameletti (2010). "A Unified Statistical Approach for Simulation, Modeling, Analysis and Mapping of Environmental Data." SIMULATION 86(3): 139-153. <doi: 10.1177/0037549709102150>
#' Cerqueti, R., Maranzano, P., & Mattera, R. (2025). \emph{Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe}. Journal of Agricultural, Biological and Environmental Statistics, 1-35.
#' Sugasawa, S., & Murakami, D. (2021). \emph{Spatially clustered regression}. Spatial Statistics, 44, 100525.
#'
#' @examples
#' \donttest{
#' #load the data
#' data(pm10)
#'
#' #extract the data
#' coordinates <- pm10$coords*1000
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
#' K <-matrix(1,ncol(z),1)
#'
#' mod1 <- STEM_Model(z=z,covariates=covariates,
#'                    coordinates=coordinates,phi=phi,K=K)
#' class(mod1)
#'
#' SCSTEM_Infocrit(mod1,2,3,distance='euclidean',crs=32632)
#' }
#'
#' @seealso See Also \code{\link{SCSTEM_Estim}} and \code{\link{pm10}}
#'
#' @keywords models spatial
#'
#'
#' @param verbose Logical. If TRUE, convergence information from SCSTEM_Estim is emitted via message(). Default is FALSE.
#'
#' @export

SCSTEM_Infocrit<-function(StemModel,mink,maxk,crs=4326,distance='geo',knn=5,init_method='AMKM',nugget_var=TRUE,precision_full_dataset=0.01,precision=0.1,regularization=0.01,phi_penalty=1,max_iter=10,share2conv=0.05,verbose=FALSE)
  {
  BIC<-list()
  AIC<-list()
  HQC<-list()
  n<-dim(StemModel$skeleton$K)[1]
  for (i in mink:maxk){
    model<-SCSTEM_Estim(StemModel = StemModel, k=i, crs=crs,distance=distance,knn=knn,init_method=init_method,nugget_var=nugget_var,precision_full_dataset=precision_full_dataset,precision=precision,regularization=regularization,phi_penalty=phi_penalty,max_iter=max_iter,share2conv=share2conv,verbose=verbose)
    logvero <- 0  # Inizializza la somma a 0 prima di entrare nel ciclo j
    for (j in 1:i) {
      logvero <- logvero + model$fit_list[[j]]$estimates$loglik
      }
    BIC <- c(BIC, -2 * logvero + i * log(n))
    AIC <- c(AIC,-2 * logvero + i * 2)
    HQC <- c(HQC,-2 * logvero + i * 2 * log(log(n)))

  }
  l<-cbind(BIC,AIC,HQC)
  rownames(l)<-paste(seq(mink,maxk),'cluster')
  return(l)
}
