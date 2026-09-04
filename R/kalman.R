#' @keywords internal
#' @noRd

`kalman` <-
  function (z, coordinates, p, n, d, r, phi_j, max.iter, precision, covariates, Gdiag, Sigmaetadiag, cov.spat,distance,regularization, verbose = FALSE) {


    zz   = stats::ts(z)
    if(distance=='euclidean'){dist = as.matrix(stats::dist(coordinates,diag=TRUE))} #distance matrix}
    if(distance=='geo'){dist = as.matrix(geodist::geodist(coordinates,measure='geodesic'))} #distance matrix}

    ### Missing observations, following Durbin and Koopman (2012), Sect. 4.10.
    ### The observed-row index of every time point is computed once: the filter
    ### restricts the measurement equation to those rows, and the M-step
    ### completes the sufficient statistics over the missing ones.
    obs_ix = stem_obs_index(z)

    ####################
    ###Model definition
    ###if you want to check the model use phi_j=phi_start
    ####################
    SSmodel  = list(z	= zz,
                    Fmat 	= phi_j$K,
                    Gmat 	= phi_j$G,
                    Vmat 	= phi_j$sigma2omega * cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist),
                    Wmat 	= phi_j$Sigmaeta,
                    m0   		= t(phi_j$m0),
                    C0   		= phi_j$C0,
                    phi  		= phi_j,
                    XXX  		= covariates,
                    beta 	= phi_j$beta,
                    flag.cov 	= TRUE,
                    n		= n,
                    p		= p,
                    d		= d,
                    m 		= NA,
                    C		= NA,
                    loglik		= NA
    )

    ####################
    ###kalman filtering and smoothing
    ####################
    mod1.filter   	= filtering(SSmodel)
    mod1.smoother 	= smoothing(mod1.filter)


    ####calculating  B0 e P_0_n=C_0* using kalman filter, smoother and initial values
    R1    = SSmodel$Gmat %*% SSmodel$C0 %*% t(SSmodel$Gmat) + SSmodel$Wmat  #P_1^0
    B0    = SSmodel$C0 %*% SSmodel$Gmat %*% 	solve(diag(regularization,nrow(R1))+R1)
    P_0_n = SSmodel$C0 + B0 %*% (mod1.smoother$C[[1]] - R1 ) %*% t(B0)

    ########################
    ###Define the elements for B_function (output of mod1$filter)
    ###Use B_function (see functions.R) ---> list of B_t (t=1,...,n)
    ########################
    nobs 	= n
    m    	= mod1.filter$m
    C    	= mod1.filter$C
    B    	= list()

    for (tt in (nobs-1):1) {
      nextstep = B_function(
        m     = matrix(m[tt,],nrow=1),
        C     = C[[tt]],
        Gmatx = SSmodel$Gmat,
        Wtx   = SSmodel$Wmat,
        mx    = matrix(m[tt+1,],nrow=1),
        Cx    = C[[tt+1]])

      m[tt,]  = nextstep$ms #equal to mod1.smoother$m
      C[[tt]] = nextstep$Cs #equal to mod1.smoother$C
      B[[tt]] = nextstep$B  #what I need
    }

    #B_n=C_n * t(Gmat) * solve(R_n+1)
    #where  R_n+1=Gmat * C_n * t(Gmat) +Wmat

    B[[nobs]] = mod1.filter$C[[nobs]] %*% t(SSmodel$Gmat) %*%
      solve(diag(regularization,nrow(SSmodel$Gmat %*% mod1.filter$C[[nobs]] %*%
                                       t(SSmodel$Gmat)+SSmodel$Wmat))+SSmodel$Gmat %*% mod1.filter$C[[nobs]] %*%t(SSmodel$Gmat)+SSmodel$Wmat)

    ################################
    ###Ricursion for LAG ONE COVARIANCE SMOOTHER
    ################################
    CCC = list() #empty list for lag one covariance values

    ###FIRST STEP: Calculate C*_{n,n-1} (start of the iterative procedure --> CCC[[n]])
    ###Pn_n-1=R_n=G_n%*%C_n-1%*%t(G_n)+W_n
    ###A_n=R_n%*%F_n%*%solve(t(F_n)%*%R_n%*%F_n+V_n)
    ### The gain that starts the lag-one recursion has to be the one actually
    ### used at time n, so it is built on the observed rows only (Durbin and
    ### Koopman 2012, Sect. 4.10). With nothing observed at time n the gain is
    ### zero, which is the Z_t = 0 device of that section.
    oi_n            = obs_ix$idx[[nobs]]
    R_n 			= SSmodel$Gmat %*% mod1.filter$C[[nobs-1]] %*%t(SSmodel$Gmat)+ SSmodel$Wmat
    if (length(oi_n) == 0L) {
      CCC[[nobs]]   = SSmodel$Gmat %*% mod1.filter$C[[nobs-1]]
    } else {
      Fmat_n        = SSmodel$Fmat[, oi_n, drop = FALSE]
      Vmat_n        = SSmodel$Vmat[oi_n, oi_n, drop = FALSE]
      Q_n           = t(Fmat_n) %*% R_n %*% Fmat_n + Vmat_n
      A_n 			= R_n %*% Fmat_n %*% solve(diag(regularization,nrow(Q_n)) + Q_n)
      CCC[[nobs]] 	= (diag(p)-A_n %*% t(Fmat_n)) %*% SSmodel$Gmat %*% mod1.filter$C[[nobs-1]]
    }

    ################################
    ###loop for calculating the other elements of CCC using cov_lagone fucntion (see function.R)
    ################################
    for (tt in (nobs):2) {
      passiCCC =cov_lagone(
        C_t_minus_1 	= mod1.filter$C[[tt-1]],
        B_t_minus_1 	= B[[tt-1]],
        Gmat        		= SSmodel$Gmat,
        if (tt>=3)
          B_t_minus_2= B[[tt-2]] else B_t_minus_2 = B0,
        CCCx       		 = CCC[[tt]])

      CCC[[tt-1]] = passiCCC$cov
    }

    ################################################################
    #Parameters estimates and calculate of the elements of Q function
    ################################################################
    ###PARAMETER 1: m0=y_0_n
    y_0_n = t(SSmodel$m0) + B0 %*% (t(matrix(mod1.smoother$m[1,],nrow=1))- SSmodel$Gmat %*% t(SSmodel$m0))
    m0_j  = y_0_n

    ####Q element n.2 (NB: y0_n=m0)
    Q_addendo2 = Q_function_addendo2(C0=phi_j$C0, m_0=m0_j, P_0_n=P_0_n, y_0_n=y_0_n)

    ###S11, S10 and S00
    S11list = list()
    for (tt in 1:nobs) {
      S11list[[tt]] = t(matrix(mod1.smoother$m[tt,],nrow=1)) %*% matrix(mod1.smoother$m[tt,],nrow=1) +  mod1.smoother$C[[tt]]
    }
    S11 = sumMatrices(S11list)

    S00list = list()
    S00list[[1]] = m0_j %*% t(m0_j) + P_0_n
    for (tt in 2:nobs) {
      S00list[[tt]] = t(matrix(mod1.smoother$m[tt-1,],nrow=1)) %*% matrix(mod1.smoother$m[tt-1,],nrow=1) +  mod1.smoother$C[[tt-1]]
    }
    S00 = sumMatrices(S00list)

    S10list = list()
    S10list[[1]] = t(matrix(mod1.smoother$m[1,],nrow=1)) %*% t(m0_j) + CCC[[1]]
    for (tt in 2:nobs) {
      S10list[[tt]] = t(matrix(mod1.smoother$m[tt,],nrow=1)) %*% matrix(mod1.smoother$m[tt-1,],nrow=1) +  CCC[[tt]]
    }
    S10 = sumMatrices(S10list)

    ################################
    ###PARAMETER 2: Sigmaeta
    #Sigmaeta_j = (S11 - S10 %*% t(phi_j$G) - phi_j$G %*% t(S10) + phi_j$G %*% S00 %*% t(phi_j$G))/nobs
    Sigmaeta_j = (S11 - S10 %*% solve(S00) %*% t(S10))/nobs
    if(Gdiag & Sigmaetadiag & p>1) Sigmaeta_j = diag(diag(Sigmaeta_j))

    if (Sigmaetadiag) {
      Sigmaeta_j = matrix(0,p,p)
      num = (S11 - S10 %*% t(phi_j$G) - phi_j$G %*% t(S10) + phi_j$G %*% S00 %*% t(phi_j$G))
      for (i in 1:p) { Sigmaeta_j[i,i] = num[i,i] / nobs }
    }

    ################################
    ###PARAMETER 3: G
    G_j = S10 %*% solve(S00)
    if(Gdiag & Sigmaetadiag & p>1)  G_j = diag(diag(G_j))

    if (Gdiag) {
      G_j = matrix(0,p,p)
      num = solve(diag(regularization,nrow(Sigmaeta_j))+Sigmaeta_j) %*% S10
      den = solve(diag(regularization,nrow(Sigmaeta_j))+Sigmaeta_j) %*% S00
      for (i in 1:p) { G_j[i,i] = num[i,i] / den[i,i]  }
    }


    ###Q element n.3
    Q_addendo3 = Q_function_addendo3(Sigmaeta=Sigmaeta_j, G=G_j, n=n, S11=S11, S00=S00, S10=S10)

    #############################
    ###PARAMETER 4: sigma2omega
    ###sigma2omega=sigma2epsilon if logb=0 (exp(logb)=1)
    #W=eq 13 Fasso'-Cameletti
    #W=BB
    ### E-step for the measurement equation. In the complete case this is
    ###   BB_t = Z C^s_t Z' + r_t r_t' ,
    ### the conditional second moment of the measurement error given the data.
    ### With missing entries the same quantity has to be COMPLETED, because the
    ### EM algorithm maximizes the expected complete-data log-likelihood. The
    ### observed block keeps the form above; the missing block is predicted from
    ### the observed one by Gaussian conditioning on Sigma_e -- the same algebra
    ### as kriging at a fixed time point -- and its conditional variance Omega_t
    ### is added back. Writing S_t for the matrix that lifts a residual defined
    ### on the observed rows to the full vector (identity on the observed rows,
    ### P_t on the missing ones),
    ###
    ###   BB_t = S_t Zo C^s_t Zo' S_t' + Omega_t + (S_t r_o)(S_t r_o)' ,
    ###
    ### which reduces to the complete-data formula when nothing is missing and
    ### to Sigma_e when nothing is observed. Note that the conditional mean of a
    ### missing element is NOT the signal alone: that would be exact only for a
    ### diagonal Sigma_e. The completed observations zhat are kept for the beta
    ### update below.
    blocks  = stem_blocks_cache(SSmodel$Vmat, d, regularization)
    Zmat    = t(SSmodel$Fmat)                     # d x p loading matrix
    BB_list = list()
    zhat    = matrix(NA_real_, nobs, d)

    for (tt in 1:nobs) {
      msm = matrix(mod1.smoother$m[tt,], ncol = 1)          # p x 1
      Csm = mod1.smoother$C[[tt]]
      Xb  = covariates[,,tt] %*% phi_j$beta                 # d x 1
      sig = Zmat %*% msm                                    # d x 1
      oi  = obs_ix$idx[[tt]]

      if (obs_ix$complete[tt]) {
        r_t           = matrix(zz[tt,], ncol = 1) - Xb - sig
        BB_list[[tt]] = Zmat %*% Csm %*% t(Zmat) + r_t %*% t(r_t)
        zhat[tt,]     = zz[tt,]
      } else {
        bl = blocks(obs_ix$key[tt], oi)
        if (length(oi) == 0L) {
          BB_list[[tt]] = bl$Omega
          zhat[tt,]     = Xb + sig
        } else {
          Zo   = Zmat[oi, , drop = FALSE]
          r_o  = matrix(zz[tt, oi], ncol = 1) - Xb[oi, , drop = FALSE] -
                   sig[oi, , drop = FALSE]
          SZo  = bl$S %*% Zo
          rhat = bl$S %*% r_o
          BB_list[[tt]] = SZo %*% Csm %*% t(SZo) + bl$Omega + rhat %*% t(rhat)
          zhat[tt,]     = Xb + sig + rhat
        }
      }
    }
    BB = sumMatrices(BB_list)

    D = solve(diag(regularization,nrow(cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist)))+cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist)) %*% BB
    #sigma2omega_j=tr(sigmaeinersa*W) in Fasso Cameletti 12
    ### The divisor stays n*d, the COMPLETE-data count, and not the number of
    ### observed values: the EM algorithm maximizes the expected complete-data
    ### log-likelihood, and BB above already carries the conditional variance of
    ### whatever was not observed.
    sigma2omega_j = sum(diag(D))/(n*d)

    #############################
    ###PARAMETER 5: beta coefficients
    ############################
    #\sum_t X_t^\prime \Sigma_e^-1 v_t
    #v_t=z_t-K_t y_t
    Sigmae_inversa = solve(diag(regularization,nrow(sigma2omega_j * cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist)))+sigma2omega_j * cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist))

    ### zhat is the observation vector completed by the E-step: it equals z
    ### wherever z was observed and the conditional expectation of the missing
    ### element given everything observed elsewhere. The design matrix and
    ### Sigma_e stay at their full dimension, as they must, since the M-step
    ### maximizes the expected COMPLETE-data log-likelihood.
    vvt_list = list ()
    for (tt in 1:nobs) {
      vt 		=  matrix(zhat[tt,],ncol=1) - t(SSmodel$Fmat) %*% mod1.smoother$m[tt,]
      vvt_list[[tt]] 	= t(covariates[,,tt]) %*%  Sigmae_inversa  %*% vt
    }
    v = sumMatrices(vvt_list)


    MM_list = list()
    for (tt in 1:nobs) {
      MM_list[[tt]] 	= t(covariates[,,tt]) %*% Sigmae_inversa %*% covariates[,,tt]
    }
    MM = sumMatrices(MM_list)

    if(det(MM) != 0) {beta_j = solve(diag(regularization,nrow(MM))+MM) %*% v}
    if(det(MM)  < 10^(-7)) {warning("Error in beta estimation! The matrix can not be inverted!!!!", call. = FALSE)}

    #############################
    ###PARAMETER 6: theta e logb
    ##exponential spatial covariance function
    ##Newton raphson algorithm
    #############################
    n_iter_NR = 1
    n_iter_Hess.list = c()
    convergence_NR = FALSE

    logb     = phi_j$logb
    logtheta = phi_j$logtheta

    while(!convergence_NR && n_iter_NR < 50) {

      Q_addendo1_old = Q_function_addendo1(sigma2omega=sigma2omega_j ,n=n, Sigmastar=do.call(cov.spat, list(d=d, logb=logb, logtheta=logtheta, dist=dist)),B=BB)
      Q_prev = Q_addendo1_old+Q_addendo2+Q_addendo3

      cond.hessiana = FALSE
      n_iter_Hess   = 1
      logtheta.iniz = logtheta

      while(!cond.hessiana && n_iter_Hess < 30) {
        cov.spat.mat = do.call(cov.spat, list(logb=logb,d=d,logtheta=logtheta,dist=dist))

        derivata_prima_logtheta 	= d1_Q(
          n	= n,
          X	= cov.spat.mat,
          d1_X  = d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
          sigma2omega = sigma2omega_j,
          B       = BB)
        derivata_seconda_logtheta 	= d2_Q(
          n	= n,
          X      = cov.spat.mat,
          d1_X	= d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
          d2_X = d2_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
          sigma2omega = sigma2omega_j,
          B      = BB)

        derivata_prima_logb    		= d1_Q(
          n	 = n,
          X       = cov.spat.mat,
          d1_X  = d1_Sigmastar_logb.exp(logb=logb,d=d),
          sigma2omega = sigma2omega_j,
          B       = BB)

        derivata_seconda_logb    	= d2_Q(
          n	= n,
          X 	= cov.spat.mat,
          d1_X	= d1_Sigmastar_logb.exp(logb=logb,d=d),
          d2_X	= d2_Sigmastar_logb.exp(logb=logb,d=d),
          sigma2omega = sigma2omega_j,
          B	= BB)

        derivata_mista        		=  d12_Q(
          n	= n,
          X	= cov.spat.mat,
          d1_X_theta  = d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
          d1_X_logb   = d1_Sigmastar_logb.exp(logb=logb,d=d),
          sigma2omega = sigma2omega_j,
          B	= BB)


        hessiana  = matrix( c(derivata_seconda_logtheta, derivata_mista, derivata_mista, derivata_seconda_logb),2,2)
        ### Robustness: with isTRUE() a non-finite determinant counts as "not
        ### well conditioned" instead of turning the loop condition into NA.
        cond.hessiana = isTRUE(det(hessiana) > 10^(-3))
        ### A Hessian with non-finite entries cannot be repaired by the grid
        ### search below: leave the inner loop and let the guarded Newton step
        ### decide what to do.
        if(!all(is.finite(hessiana))) break

        if(!cond.hessiana) {
          kk1=10
          kk2=10
          ### Guard against an underflowed scale, which would make seq() start
          ### at exactly 0 and produce log(0) = -Inf in the grid.
          theta.scale = max(exp(logtheta), .Machine$double.xmin)
          b.scale     = max(exp(logb),     .Machine$double.xmin)
          logtheta_vec = log(seq((0.01*theta.scale),(10*theta.scale),length=kk1))
          logb_vec  = log(seq((0.01*b.scale),(10*b.scale),length=kk2))
          QQ=matrix(NA,kk1,kk2)

          for(i in 1:kk1){
            for(j in 1:kk2){
              QQ[i,j] = Q_function_addendo1(sigma2omega=sigma2omega_j , n=n ,
                                            Sigmastar=cov.spat(d=d , logb=logb_vec[j] , logtheta=logtheta_vec[i], dist=dist),B=BB)
            }
          }

          val.col.min = apply(QQ,2,min)
          colonna = which.min(val.col.min)
          righe = apply(QQ,2,which.min)
          riga = righe[colonna]

          logtheta = logtheta_vec[riga]
          logb     = logb_vec[colonna]
        }

        n_iter_Hess = n_iter_Hess + 1
      }

      Q_addendo1_prev =	 Q_function_addendo1(sigma2omega=sigma2omega_j, n=n,
                                             Sigmastar=cov.spat(d=d, logb=logb, logtheta=logtheta, dist=dist),B=BB)

      gradiente = matrix( c(derivata_prima_logtheta, derivata_prima_logb),2,1)

      logtheta.logb_old = matrix(c(logtheta , logb),2,1)

      ### Robustness of the Newton step. A singular or non-finite Hessian, or a
      ### step that leaves the finite range, used to propagate NaN into the
      ### convergence test and abort the whole fit with "missing value where
      ### TRUE/FALSE needed". This is harmless for a single pooled fit, where
      ### such subsets are rare, but it makes the clusterwise algorithm brittle,
      ### since every candidate partition produces a new subset of locations.
      ### The step is now validated and, when it is not usable, the loop exits
      ### keeping the last valid iterate.
      delta = try(solve(hessiana + diag(regularization, nrow(hessiana))) %*% gradiente,
                  silent = TRUE)
      if(inherits(delta, "try-error") || !all(is.finite(delta))) {
        n_iter_Hess.list[[n_iter_NR]] = n_iter_Hess - 1
        n_iter_NR = n_iter_NR + 1
        break
      }
      logtheta.logb_new = logtheta.logb_old - delta
      if(!all(is.finite(logtheta.logb_new))) {
        n_iter_Hess.list[[n_iter_NR]] = n_iter_Hess - 1
        n_iter_NR = n_iter_NR + 1
        break
      }

      dist_rel_num = sqrt(t(logtheta.logb_new - logtheta.logb_old) %*% (logtheta.logb_new - logtheta.logb_old))
      ### Floor on the denominator: the relative criterion is undefined when the
      ### current iterate sits exactly at the origin.
      dist_rel_den = max(sqrt(t(logtheta.logb_old) %*% logtheta.logb_old),
                         .Machine$double.eps)
      dist_rel = dist_rel_num / dist_rel_den

      convergence_NR = isTRUE(as.logical(dist_rel < precision))
      logb  = logtheta.logb_new[2,]
      logtheta = logtheta.logb_new[1,]

      n_iter_Hess.list[[n_iter_NR]] = n_iter_Hess - 1
      n_iter_NR = n_iter_NR + 1
      if (isTRUE(verbose)) message("*** NR Algorithm - iteration n. ", n_iter_NR - 1)


    } #end while loop


    logtheta_j 	= logtheta
    logb_j     	= logb


    Q_addendo1 = Q_function_addendo1(sigma2omega=sigma2omega_j, n=n, Sigmastar=cov.spat(d=d, logb=logb_j, logtheta=logtheta_j, dist=dist),B=BB)
    Q_new      	 = Q_addendo1 + Q_addendo2 + Q_addendo3


    ###########################
    ###New parameter vector (k e C0 don't change)
    ############################

    phi_jj = list(
      loglik	        	= mod1.filter$loglik,
      K      	    	= phi_j$K,
      sigma2omega	= sigma2omega_j,
      logtheta      	= logtheta_j,
      logb       		= logb_j,
      beta       		= beta_j,
      G          		= G_j,
      Sigmaeta   	= Sigmaeta_j,
      m0         		= m0_j,
      C0         		= phi_j$C0)

    return(list(phi    = phi_jj,
                Q_prev = Q_prev,
                Q_new  = Q_new,
                n_iter_NR =  n_iter_NR - 1,
                m.smoother = mod1.smoother$m))
    #m.filter   = mod1.filter$m,
    #c.smoother = mod1.smoother$C,
    #c.filter   = mod1.filter$C))
  }

