#' STEM Clustered Regression
#'
#'
#' @description The function \code{SCSTEM_Estim} computes the spatio-temporal clustered regression with the STEM_estimation method for the regression.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param precision  a small positive number used for the EM algorithm convergence in every step of the clustered regression. Default is equal to 0.1. See \code{DETAILS} below.
#' @param precision_full_dataset a small positive number used for the EM algorithm convergence at the inizialization step of the EM algorithm. Default is equal to 0.1. See \code{DETAILS} below.
#' @param max_iter maximum number of iterations for the spatio temporal clustered regression algorithm.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param regularization a small positive number to be added to the digonal of the matrices matrices that need to be inverted . Default is set to 0.01
#' @param crs Integer value. Coordinate reference sySTEM_ something suitable as input to st_crs.command from the sf package (see its documentation for details). Default is set to 4326.
#' @param knn Integer value. The number of nearest neighbour to be taken into account for te spatial penalty, Default is set to 5
#' @param k Integer value. The number of clusters for the clustered regression. Default is set to 3.
#' @param init_method Character. Must be one of: 'AMKM' or 'K-means'. If init_method='AMKM', the Adjacent Matrix K-Means clustering is performed. If method='K-means', K-means clustering is performed.
#' @param nugget_var Logical. If FALSE it returns the Sugasawa clustered regression. If TRUE the STEM spatio temporal clustered regression is performed. Default is set to TRUE
#' @param phi_penalty a small positive number. It is the spatial penalty weight for the clustered regression. Default is set to 1.
#' @param share2conv a small positive number. It is the minimum percentage of observations that must change clusters for the algorithm not to converge. Default is set equal to 0.05
#'
#'
#'
#' @return The function returns a list given by:
#' \itemize{
#' \code{phi.hat} is a matrix with the parameter ML estimates (\code{sigma2omega}, \code{beta}, \code{G}, \code{Sigmaeta}, \code{m0}, \code{C0}, \code{theta}, \code{sigma2eps}) for each cluster.
#' \code{clusters} is a df with 3 columns. Longitude, Latitude and cluster for each observation.
#' \code{fit_list} a list that contains the single cluster model output.
#' }
#'
#'
#' @details This function estimates the spatio temporal clustered regression via the STEM algorithm.
#'
#'
#' @author Francesco Caccia  < francesco.caccia2000@gmail.com >
#'
#' @references Amisigo, B.A., Van De Giesen, N.C. (2005) \emph{Using a spatio-temporal dynamic state-space model with the EM algorithm to patch gaps in daily riverflow series}. Hydrology and Earth System Sciences 9, 209--224.
#'
#' Fasso, A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring using heterogeneous networks}. Environmetrics 18, 245--264. <doi: 10.1002/env.837>
#'
#' Fasso', A., Cameletti, M. (2007) \emph{A general spatio-temporal model for environmental data}. Tech.rep. n.27 \emph{Graspa} - The Italian Group of Environmental Statistics.
#'
#' Fassò, A. and M. Cameletti (2010). "A Unified Statistical Approach for Simulation, Modeling, Analysis and Mapping of Environmental Data." SIMULATION 86(3): 139-153. <doi: 10.1177/0037549709102150>
#'
#' Cerqueti, R., Maranzano, P., & Mattera, R. (2025). \emph{Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe}. Journal of Agricultural, Biological and Environmental Statistics, 1-35.
#'
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
#' SCSTEM_Estim(StemModel = mod1,distance='euclidean',crs=32632)
#' }
#'
#' @seealso See Also \code{\link{STEM_Model}} and \code{\link{pm10}}
#'
#' @keywords models spatial
#'
#'
#' @param verbose Logical. If TRUE, convergence information is emitted via message(). Default is FALSE.
#' @param plot_clusters Logical. If TRUE, the final cluster map is plotted. Default is FALSE.
#'
#' @export

SCSTEM_Estim<-function(StemModel,crs=4326,distance='geo',knn=5,k=3,init_method='AMKM',nugget_var=TRUE,precision_full_dataset=0.01,precision=0.1,regularization=0.01,phi_penalty=1,max_iter=10,share2conv=0.05,verbose=FALSE,plot_clusters=FALSE){

  l<-list()
  z<-StemModel$data$z
  coordinates<-StemModel$data$coordinates
  covariates<-StemModel$data$covariates
  phi<-StemModel$skeleton$phi
  K<-StemModel$skeleton$K
  n<-ncol(z)
  day<-nrow(z)
  ncov<-ncol(covariates)
  mod1 <- STEM_Model(z = z, covariates = covariates,coordinates = coordinates, phi = phi, K = K)
  mod1.est <- STEM_Estimation(mod1, precision = precision_full_dataset,distance=distance,regularization = regularization)
  if(k==1){
    l<-list(unlist(mod1.est$estimates$phi.hat),mod1.est)
    names(l)<-c('phi_hat','fit_list')
    }
  if (k>=2){
    nb <- spdep::knn2nb(spdep::knearneigh(coordinates, k = knn))  #  vicini più prossimi
    W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
    block_size <- day  # Numero di righe per blocco
    covariates3<-as.data.frame(covariates)
    # Aggiungiamo un identificatore di blocco
    covariates3 <- covariates3 %>%
      dplyr::mutate(block = rep(1:(nrow(covariates3) / block_size), each = block_size))
      # Calcoliamo la media per ogni blocco e colonna
    new_dataset <- covariates3 %>%
      dplyr::group_by(.data$block) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)))  # Rimuoviamo la colonna block
    new_dataset<-new_dataset[,-1]
    # Inizializzazione dei cluster spazio-temporali
    dati<-cbind(coordinates,new_dataset[,-c(1)])
    clusters <-SCDA::SC_AMKM(Data_sf=sf::st_as_sf(dati, coords = c(1, 2)),Method = init_method,MinNc =k,MaxNc = k,IndexCol = 0,CRS = crs)  # Cluster iniziali k
    clusters<-clusters$df$cluster
    if (sum((table(clusters)<=1))>0){
      stop('A cluster contains 1 observation, change init_method or try a smaller k')
    }
    G<-max(as.numeric(clusters))

    phi_penalty <-phi_penalty # Penalizzazione spaziale

    max_iter <- max_iter # Numero massimo di iterazioni
    for (iter in 1:max_iter) {
      if (isTRUE(verbose)) message('SCSTEM iteration: ', iter)
      # Step A: Stima dei parametri spazio-temporali per ogni cluster
      phi_hat <- matrix(0, G, 10 )#10 npar
      colnames(phi_hat)<-names(unlist(mod1.est$estimates$phi.hat))
      loglik<-matrix(0,G,1)
      fit_list<-list()
      for (g in 1:G) {
        if (!(g %in% clusters)) next
        indices <- which(clusters == g)  # Seleziona le stazioni appartenenti al cluster g
        # Estraggo le righe temporali per ogni stazione
        idx_list <- unlist(lapply(indices, function(i) {
          start_idx <- (i - 1) * day + 1
          end_idx <- i * day
          start_idx:end_idx
          }))
        K_cluster <- matrix(1, length(indices), 1)

        model_cluster <- STEM_Model(z = z[, indices], covariates = covariates[idx_list, ],
                                    coordinates = coordinates[indices, ], phi = phi, K = K_cluster)

        fit <- STEM_Estimation(model_cluster,precision=precision,regularization = regularization,distance=distance)
        phi_hat[g, ] <- unlist(fit$estimates$phi.hat)
        fit_list[[g]] <-fit
        }

      # Step B: Riassegnazione ai cluster basata su componente spazio-temporale
      new_clusters <- clusters
      likelihoods <- matrix(0, n,G)
      a<-matrix(0,n,G)
      b<-matrix(0,n,G)
      for (i in 1:ncol(z)) {
        for (g in 1:G) {
          #Se il cluster g non ha osservazioni, skip
          if (!(g %in% clusters)) next
          # Parametri stimati per il cluster g
          beta_g <- phi_hat[g, 2:(ncov+1)]  # Coefficienti delle covariate
          sigma2_omega_g <- mod1.est$estimates$phi.hat$sigma2omega
          sigma2eps_g<-mod1.est$estimates$phi.hat$sigma2eps
          # Seleziona la serie temporale della stazione i
          z_i <- z[, i]
          X_i <- covariates[((i-1) * day + 1):(i * day), ]
          y_t <- fit_list[[g]]$estimates$y.smoothed
          l_y_t<-c(phi_hat[g,(ncov+4)],y_t[1:(day-1)])
          # Calcola la previsione condizionata
          f_t <- X_i %*% beta_g + l_y_t
          var_g<-(1+(sigma2eps_g/sigma2_omega_g))

          a[i,g]<-(-day/2)*log((sigma2_omega_g*var_g))
          b[i,g]<--0.5*sum((z_i-f_t)^2/(sigma2_omega_g*var_g))
          # Calcolo della log-verosimiglianza condizionata
          likelihoods[i,g] <-(
            (-day/2)*log((sigma2_omega_g*var_g))
            -0.5*sum((z_i-f_t)^2/(sigma2_omega_g*var_g))
            -0.5*log(phi_hat[g,(ncov+5)])
            -0.5*(y_t[1]-phi_hat[g,(ncov+4)])^2
            -(day/2)*log(phi_hat[g,(ncov+3)])
            -0.5*sum((y_t-phi_hat[g,(ncov+2)]*l_y_t)^2/phi_hat[g,(ncov+4)]))
          }
        # Penalità spaziale basata sul modello di Potts
        penalty <- phi_penalty * colSums(W[i, ] * (clusters == matrix(1:G, ncol(z), G, byrow = TRUE)))
        # Massimizzazione della funzione obiettivo
        if (isTRUE(nugget_var)) {new_clusters[i] <- which.max((a[i,]+b[i,])/n + penalty)}
        if (!isTRUE(nugget_var)) {new_clusters[i] <- which.max((b[i,])/n + penalty)}
        }
      # Controllo di convergenza
      if (sum(new_clusters != clusters) == 0) {
        if (isTRUE(verbose)) message('Exact convergence met in ', iter, ' iteration.')
        break
        }
      if (((sum(new_clusters != clusters))/n) <= share2conv) {
        if (isTRUE(verbose)) message('Convergence reached: less than ', share2conv * 100, '% of the observations change cluster at iteration ', iter, '.')
        break
      }
      if (sum((table(new_clusters)<=1))>0){
        stop('Convergence not reached: A cluster contains 0 or 1 observation, a smaller number of clusters may be required')
        break
      }
      clusters <- new_clusters
      }
    clusters <- new_clusters
    if (isTRUE(plot_clusters)) plot(coordinates, col = clusters, pch = 19, main = "Cluster")
    l<-list(phi_hat,cbind(coordinates,clusters),fit_list)
    names(l)<-c('phi_hat','df','fit_list')
    }

  return(l)
}

