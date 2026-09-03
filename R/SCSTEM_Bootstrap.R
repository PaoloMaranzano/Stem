#' STEM Clustered Regression parametric bootstrap
#'
#'
#' @description The function \code{SCSTEM_Bootstrap} computes the spatio-temporal parametric bootstrap for computing the parameter standard errors of the clustered regression.
#'
#' @param SCSTEM is the \dQuote{SCSTEM_Estim } function output.
#' @param B number of bootstrap iterations.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param precision  a small positive number used for the STEM_Estimation algorithm convergence. Default is equal to 0.01.
#' @param regularization a small positive number used for the STEM_Estimation algorithm. It is the value to be added to the digonal of the hessian matrix to avoid quasi-singularity problem. Default is set to 0.01
#' @param alfa the significance level of the confidence interval. default is set to 0.05
#'
#'
#' @return The function returns a list of datframe (one for each cluster).
#' each dataframe contains the estimates, the standard error, the t statistic value, the p-value and the confidence interval for each parameter.
#'
#'
#' @details This function estimates the spatio temporal clustered regression standard error.
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
#' estim<-SCSTEM_Estim(StemModel = mod1,distance='euclidean',crs=32632)
#' SCSTEM_Bootstrap(estim,10,distance='euclidean')
#' }
#'
#' @seealso See Also \code{\link{STEM_Estimation}} and \code{\link{SCSTEM_Estim}}
#'
#' @keywords models spatial
#'
#'
#' @export
#'

SCSTEM_Bootstrap<-function(SCSTEM,B,distance='geo',precision=0.1,regularization=0.01,alfa=0.05){
  if(length(SCSTEM)==2){
    fit <- SCSTEM$fit_list
    boot<-STEM_Bootstrap(fit, B, distance = distance, precision=precision, regularization=regularization)
    npar <- 9
    boot.estimates <- matrix(NA, nrow = B, ncol = npar)
    for(b in 1:B) {
      phi.estimated <- boot$boot.output[[b]]$estimates$phi.hat
      boot.estimates[b,] <- c(phi.estimated$beta,
                              phi.estimated$sigma2eps,
                              phi.estimated$sigma2omega,
                              phi.estimated$theta,
                              phi.estimated$G,
                              phi.estimated$Sigmaeta,
                              phi.estimated$m0)
    }
    se <- sqrt(diag(stats::var(stats::na.omit(boot.estimates))))
    #create a summary table with Estimates, Standard Errors (SE) and T-statistics.
    phi.hat <- fit$estimates$phi.hat
    MLE <- c(phi.hat$beta, phi.hat$sigma2eps, phi.hat$sigma2omega,
             phi.hat$theta, phi.hat$G, phi.hat$Sigmaeta,phi.hat$m0)
    output1 <- cbind(MLE, se, MLE/se)
    colnames(output1)<- c("Estimate", "SE", "T-stat.")
    IC <- matrix(NA,nrow=npar,ncol=2)
    for(i in 1 : npar) {
      IC[i,] <- c(stats::quantile(boot.estimates[,i],alfa/2),
                  stats::quantile(boot.estimates[,i],1-(alfa/2)))
    }
    IC.df <- data.frame(
      Parameter = c('beta0','beta1','beta2','sigma2eps','sigma2omega','theta','G','sigmaeta','m0'),
      Lower = IC[,1],
      Upper = IC[,2]
    )

    names(IC.df)[2:3]<-c(paste0("Lower_", (1-alfa)*100, "_CI"),paste0("Upper_", (1-alfa)*100, "_CI"))
    params <- c('beta1','beta2','beta3','sigma2eps','sigma2omega','theta','G','sigmaeta','m0')

    # Estrai le colonne da output1
    estimates <- output1[, 1]
    SE <- output1[, 2]
    tvalue <- output1[, 3]


    # Calcola il p-value (approssimando con la normale, in assenza di gradi di libertà esatti)
    pval <- 2 * (1 - stats::pnorm(abs(tvalue)))

    # Assegna i simboli di significatività in stile summary.lm
    sig <- ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01, "**",
                         ifelse(pval < 0.05, "*",
                                ifelse(pval < 0.1, ".", ""))))

    # Crea la tabella finale unendo anche gli intervalli di confidenza e arrotondando a 5 cifre decimali
    final_table <- data.frame(
      Parameter      = params,
      Estimate       = round(estimates, 5),
      `Std. Error`   = round(SE, 5),
      `t value`      = round(tvalue, 5),
      `Pr(>|t|)`     = round(pval, 5),
      Signif.        = sig,
      Upper = round(IC.df[,2], 5),
      Lower = round(IC.df[,3], 5)
    )

    names(final_table)[7:8]<-c(paste0("Lower_", (1-alfa)*100, "_CI"),paste0("Upper_", (1-alfa)*100, "_CI"))
    final_tables<-final_table


  }
  else{
    clusters<-SCSTEM$df[,3]
    G<-max(clusters)
    fit_list<-SCSTEM$fit_list
    final_tables<-list()
    for(g in 1:G){
      fit <- (fit_list[[g]])
      boot<-STEM_Bootstrap(fit, B, distance = distance, precision=precision, regularization=regularization)
      npar <- 9
      boot.estimates <- matrix(NA, nrow = B, ncol = npar)
      for(b in 1:B) {
        phi.estimated <- boot$boot.output[[b]]$estimates$phi.hat
        boot.estimates[b,] <- c(phi.estimated$beta,
                                phi.estimated$sigma2eps,
                                phi.estimated$sigma2omega,
                                phi.estimated$theta,
                                phi.estimated$G,
                                phi.estimated$Sigmaeta,
                                phi.estimated$m0)
      }
      se <- sqrt(diag(stats::var(stats::na.omit(boot.estimates))))
      #create a summary table with Estimates, Standard Errors (SE) and T-statistics.
      phi.hat <- fit$estimates$phi.hat
      MLE <- c(phi.hat$beta, phi.hat$sigma2eps, phi.hat$sigma2omega,
               phi.hat$theta, phi.hat$G, phi.hat$Sigmaeta,phi.hat$m0)
      output1 <- cbind(MLE, se, MLE/se)
      colnames(output1)<- c("Estimate", "SE", "T-stat.")
      IC <- matrix(NA,nrow=npar,ncol=2)
      for(i in 1 : npar) {
        IC[i,] <- c(stats::quantile(boot.estimates[,i],alfa/2),
                    stats::quantile(boot.estimates[,i],1-(alfa/2)))
      }
      IC.df <- data.frame(
        Parameter = c('beta0','beta1','beta2','sigma2eps','sigma2omega','theta','G','sigmaeta','m0'),
        Lower = IC[,1],
        Upper = IC[,2]
      )

      names(IC.df)[2:3]<-c(paste0("Lower_", (1-alfa)*100, "_CI"),paste0("Upper_", (1-alfa)*100, "_CI"))
      params <- c('beta1','beta2','beta3','sigma2eps','sigma2omega','theta','G','sigmaeta','m0')

      # Estrai le colonne da output1
      estimates <- output1[, 1]
      SE <- output1[, 2]
      tvalue <- output1[, 3]


      # Calcola il p-value (approssimando con la normale, in assenza di gradi di libertà esatti)
      pval <- 2 * (1 - stats::pnorm(abs(tvalue)))

      # Assegna i simboli di significatività in stile summary.lm
      sig <- ifelse(pval < 0.001, "***",
                    ifelse(pval < 0.01, "**",
                           ifelse(pval < 0.05, "*",
                                  ifelse(pval < 0.1, ".", ""))))

      # Crea la tabella finale unendo anche gli intervalli di confidenza e arrotondando a 5 cifre decimali
      final_table <- data.frame(
        Parameter      = params,
        Estimate       = round(estimates, 5),
        `Std. Error`   = round(SE, 5),
        `t value`      = round(tvalue, 5),
        `Pr(>|t|)`     = round(pval, 5),
        Signif.        = sig,
        Upper = round(IC.df[,2], 5),
        Lower = round(IC.df[,3], 5)
      )

      names(final_table)[7:8]<-c(paste0("Lower_", (1-alfa)*100, "_CI"),paste0("Upper_", (1-alfa)*100, "_CI"))
      final_tables[[g]]<-final_table
    }
    names(final_tables) <- paste("cluster", 1:G)
  }

  return(final_tables)
}
