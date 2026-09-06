#' @keywords internal
#' @noRd

`kalman` <-
  function (z, coordinates, p, n, d, r, phi_j, max.iter, precision, covariates, Gdiag, Sigmaetadiag, cov.spat,distance,regularization, verbose = FALSE, alpha = 0, lambda = 0, penalize = NULL) {



    zz   = stats::ts(z)
    if(distance=='euclidean'){dist = as.matrix(stats::dist(coordinates,diag=TRUE))} #distance matrix}
    if(distance=='geo'){dist = as.matrix(geodist::geodist(coordinates,measure='geodesic'))} #distance matrix}

    ### Missing observations, following Durbin and Koopman (2012), Sect. 4.10.
    ### The observed-row index of every time point is computed once: the filter
    ### restricts the measurement equation to those rows, and the M-step
    ### completes the sufficient statistics over the missing ones.
    obs_ix = stem_obs_index(z)

    ### whether the spatial covariance function can be handed the exponential
    ### kernel it would otherwise rebuild; checked once rather than per
    ### Newton-Raphson iteration
    cov.spat.takes.E = "E" %in% names(formals(cov.spat))

    ### The spatial correlation at the parameter values this E-step conditions
    ### on. It enters the measurement covariance, the update of sigma2omega and
    ### the update of beta, and was being rebuilt for each of them.
    Sigmastar_j = cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist)

    ### Which regression coefficients the elastic net acts on. The default
    ### leaves the intercept alone; see R/stem-penalty.R for why that is not
    ### merely conventional here.
    pen_w = stem_penalized_index(r, penalize)

    ####################
    ###Model definition
    ###if you want to check the model use phi_j=phi_start
    ####################
    SSmodel  = list(z	= zz,
                    Fmat 	= phi_j$K,
                    Gmat 	= phi_j$G,
                    Vmat 	= phi_j$sigma2omega * Sigmastar_j,
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
    ### The three sums of outer products of the smoothed states are crossprods
    ### of the n x p matrix of those states, so each of them is one BLAS call
    ### plus one accumulation of the p x p variances:
    ###
    ###   S11 = sum_t { m_t m_t' + C_t }         = M'M       + sum_t C_t
    ###   S00 = m_0 m_0' + P_0 + sum_{t>=2} ...  = M_-'M_-   + sum_{t<n} C_t + ...
    ###   S10 = sum_t { m_t m_{t-1}' + C*_t }    = M_+'M_-   + sum_t C*_t + ...
    ###
    ### The superseded version built three lists of n matrices and summed them:
    ###
    ###   S11list = list()
    ###   for (tt in 1:nobs) S11list[[tt]] =
    ###     t(matrix(mod1.smoother$m[tt,],nrow=1)) %*%
    ###     matrix(mod1.smoother$m[tt,],nrow=1) + mod1.smoother$C[[tt]]
    ###   S11 = sumMatrices(S11list)
    ###   ... and likewise for S00 and S10.
    Msm = as.matrix(mod1.smoother$m)                     # n x p
    Csm_all = mod1.smoother$C

    S11 = crossprod(Msm) + Reduce(`+`, Csm_all)

    Mlag = Msm[seq_len(nobs - 1L), , drop = FALSE]       # m_1 ... m_{n-1}
    Mlead = Msm[-1L, , drop = FALSE]                     # m_2 ... m_n
    S00 = m0_j %*% t(m0_j) + P_0_n +
          crossprod(Mlag) + Reduce(`+`, Csm_all[seq_len(nobs - 1L)])

    S10 = matrix(Msm[1L, ], ncol = 1) %*% t(m0_j) +
          crossprod(Mlead, Mlag) + Reduce(`+`, CCC)

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
    ### The accumulation is written in matrix form rather than as a loop that
    ### builds a list of T matrices of size d x d. Over the COMPLETE time points
    ### the sum collapses, because Z does not depend on t:
    ###
    ###   sum_t { Z C^s_t Z' + r_t r_t' }  =  Z ( sum_t C^s_t ) Z'  +  R' R ,
    ###
    ### with R the T x d matrix of residuals. The first term needs one p x p
    ### accumulation and two small products; the second is a single crossprod,
    ### which BLAS performs in one call. Only the time points with a gap keep a
    ### per-step treatment, and in a monitoring network there are few of them.
    ###
    ### Besides the arithmetic, this removes the allocation: the previous
    ### version held T matrices of d x d simultaneously, which on a network of
    ### 400 stations observed for 200 periods is 256 MB per EM iteration.
    ###
    ### The superseded loop:
    ###
    ###   BB_list = list()
    ###   for (tt in 1:nobs) {
    ###     msm = matrix(mod1.smoother$m[tt,], ncol = 1)
    ###     Csm = mod1.smoother$C[[tt]]
    ###     Xb  = covariates[,,tt] %*% phi_j$beta
    ###     sig = Zmat %*% msm
    ###     oi  = obs_ix$idx[[tt]]
    ###     if (obs_ix$complete[tt]) {
    ###       r_t           = matrix(zz[tt,], ncol = 1) - Xb - sig
    ###       BB_list[[tt]] = Zmat %*% Csm %*% t(Zmat) + r_t %*% t(r_t)
    ###       zhat[tt,]     = zz[tt,]
    ###     } else { ... the same expression on the observed block ... }
    ###   }
    ###   BB = sumMatrices(BB_list)
    blocks  = stem_blocks_cache(SSmodel$Vmat, d, regularization)
    Zmat    = t(SSmodel$Fmat)                     # d x p loading matrix
    zhat    = matrix(NA_real_, nobs, d)
    msm_all = as.matrix(mod1.smoother$m)          # n x p smoothed states

    ### the regression surface for every location and time, in one pass: the
    ### covariate array is d x r x n, so each slice covariates[, j, ] is d x n
    Xb_dn = matrix(0, d, nobs)
    bvec  = as.numeric(phi_j$beta)
    for (j in seq_len(r)) Xb_dn = Xb_dn + covariates[, j, ] * bvec[j]
    sig_dn  = Zmat %*% t(msm_all)                 # d x n signal
    mu_dn   = Xb_dn + sig_dn
    Res_nd  = t(as.matrix(zz)) - mu_dn            # d x n residuals
    Res_nd  = t(Res_nd)                           # n x d

    cidx = which(obs_ix$complete)
    iidx = which(!obs_ix$complete)

    ### complete time points, in closed form
    if (length(cidx)) {
      Csum = Reduce(`+`, mod1.smoother$C[cidx])
      Rc   = Res_nd[cidx, , drop = FALSE]
      BB   = Zmat %*% Csum %*% t(Zmat) + crossprod(Rc)
      zhat[cidx, ] = as.matrix(zz)[cidx, , drop = FALSE]
    } else {
      BB = matrix(0, d, d)
    }

    ### the time points with a gap, one at a time
    for (tt in iidx) {
      Csm = mod1.smoother$C[[tt]]
      oi  = obs_ix$idx[[tt]]
      bl  = blocks(obs_ix$key[tt], oi)
      if (length(oi) == 0L) {
        BB        = BB + bl$Omega
        zhat[tt,] = mu_dn[, tt]
      } else {
        Zo   = Zmat[oi, , drop = FALSE]
        r_o  = matrix(Res_nd[tt, oi], ncol = 1)
        SZo  = bl$S %*% Zo
        rhat = bl$S %*% r_o
        BB   = BB + SZo %*% Csm %*% t(SZo) + bl$Omega + tcrossprod(rhat)
        zhat[tt,] = mu_dn[, tt] + rhat
      }
    }

    ### The spatial correlation at the CURRENT parameter values is the same
    ### matrix that built Vmat at the top of this function and that the update
    ### of beta needs below. It was being rebuilt five times in one pass -- and
    ### twice within each of the two expressions, once only to read nrow(), which
    ### is d. It is now built once, at line 31, and reused.
    ###
    ### Superseded:
    ###   D = solve(diag(regularization, nrow(cov.spat(d=d, logb=phi_j$logb,
    ###         logtheta=phi_j$logtheta, dist=dist))) +
    ###       cov.spat(d=d, logb=phi_j$logb, logtheta=phi_j$logtheta, dist=dist)) %*% BB
    D = solve(diag(regularization, d) + Sigmastar_j) %*% BB
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
    ### Superseded:
    ###   Sigmae_inversa = solve(diag(regularization,
    ###     nrow(sigma2omega_j * cov.spat(d=d, logb=phi_j$logb,
    ###          logtheta=phi_j$logtheta, dist=dist))) +
    ###     sigma2omega_j * cov.spat(d=d, logb=phi_j$logb,
    ###                              logtheta=phi_j$logtheta, dist=dist))
    Sigmae_inversa = solve(diag(regularization, d) + sigma2omega_j * Sigmastar_j)

    ### zhat is the observation vector completed by the E-step: it equals z
    ### wherever z was observed and the conditional expectation of the missing
    ### element given everything observed elsewhere. The design matrix and
    ### Sigma_e stay at their full dimension, as they must, since the M-step
    ### maximizes the expected COMPLETE-data log-likelihood.
    ### Both accumulations are sums over t of X_t' S X_t and X_t' S v_t with S
    ### constant, and a sum of that shape IS a crossprod once the per-period
    ### blocks are stacked: writing Xbig for the (n d) x r matrix whose rows are
    ### the blocks X_1, ..., X_n laid one under the other,
    ###
    ###   sum_t X_t' S X_t = Xbig' (S X)big   and   sum_t X_t' S v_t = Xbig' (S v)big ,
    ###
    ### so the whole thing is one multiplication of S against every period at
    ### once followed by two crossprods, instead of n iterations each doing two
    ### d x d products.
    ###
    ### The superseded loops:
    ###
    ###   vvt_list = list()
    ###   for (tt in 1:nobs) {
    ###     vt = matrix(zhat[tt,],ncol=1) - t(SSmodel$Fmat) %*% mod1.smoother$m[tt,]
    ###     vvt_list[[tt]] = t(covariates[,,tt]) %*% Sigmae_inversa %*% vt
    ###   }
    ###   v = sumMatrices(vvt_list)
    ###   MM_list = list()
    ###   for (tt in 1:nobs)
    ###     MM_list[[tt]] = t(covariates[,,tt]) %*% Sigmae_inversa %*% covariates[,,tt]
    ###   MM = sumMatrices(MM_list)
    ###
    ### `covariates` is d x r x n, so matrix(covariates, nrow = d) lays the
    ### periods side by side and aperm(., c(1,3,2)) stacks them by row, with the
    ### location index varying fastest -- the same ordering the residual vector
    ### below uses, which is what makes the two crossprods conformable.
    SXcat = Sigmae_inversa %*% matrix(covariates, nrow = d)      # d x (r n)
    Xbig  = matrix(aperm(covariates, c(1, 3, 2)), d * nobs, r)
    SXbig = matrix(aperm(array(SXcat, c(d, r, nobs)), c(1, 3, 2)), d * nobs, r)
    MM    = crossprod(Xbig, SXbig)

    vt_dn = t(zhat) - sig_dn                                     # d x n
    v     = crossprod(Xbig, as.vector(Sigmae_inversa %*% vt_dn))

    ### The update of beta, ordinary or regularized. With lambda = 0 this is
    ### exactly solve(diag(regularization, r) + MM) %*% v, the estimator the
    ### package has always computed; with lambda > 0 it is the elastic net of
    ### R/stem-penalty.R, closed form for ridge and coordinate descent otherwise.
    ### `beta_df` is the effective number of coefficients, which the information
    ### criteria need in place of r once a penalty is in force.
    if(det(MM)  < 10^(-7)) {warning("Error in beta estimation! The matrix can not be inverted!!!!", call. = FALSE)}
    beta_upd = stem_beta_update(M = MM, v = v, alpha = alpha, lambda = lambda,
                                w = pen_w, beta0 = phi_j$beta,
                                ridge_reg = regularization)
    beta_j  = matrix(beta_upd$beta, ncol = 1)
    beta_df = beta_upd$df

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

        ### The covariance and both of its derivatives with respect to
        ### log(theta) are built on the same exponential kernel exp(-theta h),
        ### which is a d x d elementwise exp(): it is evaluated once here and
        ### handed to all three, instead of three times over.
        Ker = exp(-exp(logtheta) * dist)
        cs_args = list(logb=logb,d=d,logtheta=logtheta,dist=dist)
        ### a covariance function that does not take the kernel builds it
        ### itself, so a user-supplied cov.spat keeps working unchanged
        if (cov.spat.takes.E) cs_args$E = Ker
        cov.spat.mat = do.call(cov.spat, cs_args)

        ### The five derivative evaluations below all need the inverse of the
        ### same matrix and its product with BB. They are formed once here and
        ### passed down, instead of each function taking its own inverse: the
        ### superseded d2_Q() alone called solve(X) seven times.
        Xi_cur  = solve(cov.spat.mat)
        XiB_cur = Xi_cur %*% BB
        d1theta = d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist,E=Ker)
        d2theta = d2_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist,E=Ker)
        ### these two are exp(logb) times the identity, and are carried as that
        ### single number: the products against them are scalings, not matrix
        ### products (see d1_Sigmastar_logb.exp)
        d1logb  = d1_Sigmastar_logb.exp(logb=logb,d=d)
        d2logb  = d2_Sigmastar_logb.exp(logb=logb,d=d)

        ### X^{-1} dSigma/dparameter enters three of the five evaluations for
        ### log(theta) and three for log(b). Formed once each.
        Pt_cur  = Xi_cur %*% d1theta
        P2t_cur = Xi_cur %*% d2theta
        Pb_cur  = d1logb * Xi_cur

        derivata_prima_logtheta 	= d1_Q(
          n	= n,
          X	= cov.spat.mat,
          d1_X  = d1theta,
          sigma2omega = sigma2omega_j,
          B       = BB, Xi = Xi_cur, XiB = XiB_cur, P = Pt_cur)
        derivata_seconda_logtheta 	= d2_Q(
          n	= n,
          X      = cov.spat.mat,
          d1_X	= d1theta,
          d2_X = d2theta,
          sigma2omega = sigma2omega_j,
          B      = BB, Xi = Xi_cur, XiB = XiB_cur, P = Pt_cur, P2 = P2t_cur)

        derivata_prima_logb    		= d1_Q(
          n	 = n,
          X       = cov.spat.mat,
          d1_X  = d1logb,
          sigma2omega = sigma2omega_j,
          B       = BB, Xi = Xi_cur, XiB = XiB_cur, P = Pb_cur)

        derivata_seconda_logb    	= d2_Q(
          n	= n,
          X 	= cov.spat.mat,
          d1_X	= d1logb,
          d2_X	= d2logb,
          sigma2omega = sigma2omega_j,
          B	= BB, Xi = Xi_cur, XiB = XiB_cur, P = Pb_cur, P2 = Pb_cur)

        derivata_mista        		=  d12_Q(
          n	= n,
          X	= cov.spat.mat,
          d1_X_theta  = d1theta,
          d1_X_logb   = d1logb,
          sigma2omega = sigma2omega_j,
          B	= BB, Xi = Xi_cur, XiB = XiB_cur, Pt = Pt_cur, Pb = Pb_cur)


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
                beta_df = beta_df,
                m.smoother = mod1.smoother$m))
    #m.filter   = mod1.filter$m,
    #c.smoother = mod1.smoother$C,
    #c.filter   = mod1.filter$C))
  }
