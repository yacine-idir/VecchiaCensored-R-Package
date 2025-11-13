# ============================================================================
# FONCTIONS   BAYESIENNES
# ============================================================================


fit_SVC_nocensored_bayesian_Vecchia=function(Y_obs,locs_obs,X_obs ,M,svc_indices=NULL,censored_indices,
                                             chains,iter_warmup,iter_sampling,parallel_chains,FIXED){

  n=length(Y_obs)

  if(sum(censored_indices)>0){
    censored_indices=which(censored_indices==1)
    non_censored_indices=setdiff(1:n,censored_indices)

    locs_cen=locs_obs[censored_indices,,drop=FALSE]
    locs_n_cen=locs_obs[non_censored_indices,,drop=FALSE]
    Y_cen=Y_obs[censored_indices]
    Y_n_cen=Y_obs[non_censored_indices]
    X_cen=X_obs[censored_indices,,drop=FALSE]
    X_n_cen=X_obs[non_censored_indices,,drop=FALSE]



    ord=GpGp::order_maxmin(locs_n_cen)
    locs_n_cen=locs_n_cen[ord,]
    Y_n_cen= Y_n_cen[ord]
    X_n_cen=X_n_cen[ord,]


    all_locs=rbind(locs_n_cen,locs_cen)
    all_X=rbind(X_n_cen,X_cen)
    all_Y=c(Y_n_cen,Y_cen)

    censored_indices=c(rep(0,length(non_censored_indices)) ,rep(1,length(censored_indices)))
    ## NN modified
    M=min(M,n-1)

    NN=GpGp::find_ordered_nn(all_locs[,1:2],M)

    dmat <- Rfast::Dist(all_locs)

    dmat=dmat[(length(Y_n_cen)+1):n,1:length(Y_n_cen),drop=FALSE]  # keep the distance between censored and non censored
    for(i in 1:nrow(dmat)){

      NN[i+length(non_censored_indices),2:(M+1)]=order(dmat[i,])[1:M] # change only rows of censored data
    }

    NN[which(is.na(NN))]=0

    nb_NN=apply(NN,MARGIN = 1,function(x) sum(x!=0))






  }else{

    ord=GpGp::order_maxmin(locs_obs)


    all_locs=locs_obs[ord,]
    all_X=X_obs[ord,]
    all_Y=Y_obs[ord]

    M=min(M,n-1)

    NN=GpGp::find_ordered_nn(all_locs,M)

    NN[which(is.na(NN))]=0

    nb_NN=apply(NN,MARGIN = 1,function(x) sum(x!=0))


  }

  gpgp_restricted=GpGp::fit_model(y=all_Y,locs =all_locs,X = all_X,covfun_name = "exponential_isotropic",silent=TRUE,m_seq = c(10) )


  if(is.null(svc_indices)){
    p=1
    A=array(data = 1,dim = c(p,n,n))
  }else{

    p=length(svc_indices)

    A=array(data = 0,dim = c(p,n,n))
    X_svc=all_X[,svc_indices,drop=FALSE]
    for(i in 1:p){

      A[i,,]=X_svc[,i,drop=FALSE]%*%t(X_svc[,i,drop=FALSE])
    }

  }





  data_stan=  list(N=n,
                   Z=all_Y,
                   M=M+1,
                   dist= Rfast::Dist(all_locs),
                   p=p,
                   pX=ncol(all_X),
                   X=all_X,
                   uncensored_idx=censored_indices,
                   NN=NN,
                   nb_NN=nb_NN,
                   means_prior_means=as.vector(gpgp_restricted$betahat),
                   var_prior_means=as.vector(gpgp_restricted$sebeta),
                   mean_prior_sigma=gpgp_restricted$covparms[1],
                   var_prior_sigma=10,
                   mean_prior_phi=gpgp_restricted$covparms[2],
                   var_prior_phi=10,
                   mean_prior_tau=gpgp_restricted$covparms[1]*gpgp_restricted$covparms[3],
                   var_prior_tau=5,
                   A=A
  )


  # model <- cmdstan_model("reduce_sum_claude.stan", cpp_options = list(stan_threads = TRUE))
  #MCMC <- model$sample(
  #          data = data_stan,
  #chains = chains,  init = init_fun,
  #parallel_chains = parallel_chains,
  # threads_per_chain = threads_per_chain,iter_warmup = iter_warmup,iter_sampling = iter_sampling  # each chain can use 4 threads
  # )



  if(FIXED==TRUE){
    print("Using Fixed range across SVC")
    stan_file <- system.file("stan", "SVC_GP_vecchia_onthefly_fixed.stan",
                             package = "VecchiaCensored")
    if(stan_file == "") {
      stop("Stan file not found. Make sure SVC_GP_vecchia_onthefly_fixed.stan is in inst/stan/")
    }
    model <- cmdstan_model(stan_file)

  }else{
    print("NOT Using Fixed range across SVC")
    stan_file <- system.file("stan", "SVC_GP_cen_vecchia_onthefly.stan",
                             package = "VecchiaCensored")
    if(stan_file == "") {
      stop("Stan file not found. Make sure SVC_GP_cen_vecchia_onthefly.stan is in inst/stan/")
    }
    model <- cmdstan_model(stan_file)
  }
  MCMC <- model$sample(
    data = data_stan,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,iter_sampling = iter_sampling  # each chain can use 4 threads
  )



  #  model <- stan_model("SVC_GP_cen_vecchia_onthefly.stan")
  # MCMC <- sampling(model, data = data_stan, chains = 5, iter = 500,cores=5)#,init=init_fun


  return(list(MCMC,data_stan))

}


# ============================================================================
# FONCTIONS  FREQUENTISTES
# ============================================================================


fit_SVC_nocensored_freq_Vecchia <- function(Y_obs, locs_obs, X_obs, M, svc_indices, fixed_range ) {

  ord <- GpGp::order_maxmin(locs_obs)
  Y_obs <- Y_obs[ord]
  locs_obs <- locs_obs[ord, , drop = FALSE]
  X_obs <- X_obs[ord, , drop = FALSE]
  NN <- GpGp::find_ordered_nn(locs_obs, M)

  lm_fit <- lm(Y_obs ~ X_obs - 1)
  mean_start <- as.vector(lm_fit$coefficients)

  p <- ncol(X_obs)
  sigma0 <- var(Y_obs) / (20 * p)

  dmat <- Rfast::Dist(locs_obs)
  phi0 <- 0.5 * mean(dmat[lower.tri(dmat)])

  n_svc <- length(svc_indices)

  # Adapt parameter structure based on fixed_range
  if (fixed_range && n_svc > 1) {
    # Structure: [sigma1, sigma2, ..., sigmaN, phi_shared, nugget]
    cov_start <- c(rep(sigma0, n_svc), phi0, sigma0 * 1e-2)
    lower <- c(rep(-Inf, length(mean_start)), rep(1e-6, length(cov_start)))
    upper <- c(rep(Inf, length(mean_start)), rep(Inf, length(cov_start)))
  } else {
    # Standard structure: [sigma1, phi1, sigma2, phi2, ..., nugget]
    cov_start <- c(rep(c(sigma0, phi0), n_svc), sigma0 * 1e-2)
    lower <- c(rep(-Inf, length(mean_start)), rep(1e-6, length(cov_start)))
    upper <- c(rep(Inf, length(mean_start)), rep(Inf, length(cov_start)))
  }

  par_start <- c(mean_start, cov_start)
  # Negative log-likelihood to minimize
  objective <- function(par) {
    mean_par <- par[1:length(mean_start)]
    cov_par <- par[(length(mean_start) + 1):length(par)]

    # Construct list of covariance blocks
    cov_list <- vector("list", n_svc + 1)

    if (fixed_range && n_svc > 1) {
      # Extract shared phi (second to last parameter)
      phi_shared <- cov_par[n_svc + 1]

      for (i in seq_len(n_svc)) {
        cov_list[[i]] <- c(sigma_sq = cov_par[i], phi = phi_shared)
      }
      cov_list$nugget <- c(nugget = cov_par[length(cov_par)])
    } else {
      # Standard: every two elements correspond to one svc component
      for (i in seq_len(n_svc)) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        cov_list[[i]] <- c(sigma_sq = cov_par[idx[1]], phi = cov_par[idx[2]])
      }
      cov_list$nugget <- c(nugget = cov_par[length(cov_par)])
    }

    ll <- Like_SVC_nocensored_freq_Vecchia(mean_par, cov_list, Y_obs, locs_obs, X_obs, NN, svc_indices)
    return(-ll)
  }

  objective(par = par_start)
  opt <- optim(par_start, objective, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(trace = 1, factr = 1e7))

  # Reconstruct parameters in standard format for output
  if (fixed_range && n_svc > 1) {
    mean_par <- opt$par[1:length(mean_start)]
    cov_par <- opt$par[(length(mean_start) + 1):length(opt$par)]

    phi_shared <- cov_par[n_svc + 1]
    nugget <- cov_par[length(cov_par)]

    # Reconstruct as [sigma1, phi1, sigma2, phi2, ..., nugget]
    reconstructed_par <- numeric(length(mean_start) + 2 * n_svc + 1)
    reconstructed_par[1:length(mean_start)] <- mean_par

    for (i in seq_len(n_svc)) {
      idx <- length(mean_start) + (2 * (i - 1) + 1):(2 * i)
      reconstructed_par[idx[1]] <- cov_par[i]        # sigma
      reconstructed_par[idx[2]] <- phi_shared         # phi (shared)
    }
    reconstructed_par[length(reconstructed_par)] <- nugget

    opt$par <- reconstructed_par
  }

  return(list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message
  ))
}


Like_SVC_nocensored_freq_Vecchia <- function(mean_parameters,
                                             cov_parameters,
                                             Y_obs,
                                             locs_obs,
                                             X_obs,
                                             NN,
                                             svc_indices) {

  # extract the covariate block used in the SVC‐covariance
  X_svc <- X_obs[, svc_indices, drop = FALSE]
  # build Cholesky of the precision: Λ = L %*% t(L)
  L <- precision_decomp_sparse(cov_parameters, X_svc, NN, locs_obs)

  # compute residual vector
  # (assumes mean_parameters corresponds to columns of X_obs)
  mu    <- X_obs %*% mean_parameters
  resid <- drop(Y_obs - mu)         # numeric vector of length n

  # quadratic form: resid' Λ resid = || L %*% resid ||^2
  z    <- L %*% resid               # sparse-dense matrix‐vector multiply
  quad <- sum(z * z)                # sum of squares

  # log-determinant of precision: det(Λ) = (prod diag(L))^2
  # FIXED: Use Matrix::diag for sparse matrices
  logdet_prec <- 2 * sum(log(Matrix::diag(L)))

  # dimension
  n <- length(resid)

  # constant term
  c0 <- n * log(2 * pi)

  # final log-likelihood
  ll <- 0.5 * logdet_prec - 0.5 * (c0 + quad)

  return(ll)
}




# ============================================================================
# PREDICTION  HELPER
# ============================================================================



prediction_freq=function(beta,sigma,phi,tau,M,Y_obs,X_obs,locs_obs,svc_indices,cen_indices,X_pred,locs_pred,simulation,pred_coef=FALSE){

  ##  prediction aux points censurés
  # if censored, it will change y_obs
  if(sum(cen_indices)>0){


    cen_indx=which(cen_indices==1)
    n_cen_indx=which(cen_indices==0)

    n_cen=length(cen_indx)
    n_n_cen=length(n_cen_indx)
    n_obs=length(Y_obs)


    Y_obs_n_cen=Y_obs[n_cen_indx]
    Y_obs_cen=Y_obs[cen_indx]

    X_obs_n_cen=X_obs[n_cen_indx,,drop=FALSE]
    X_obs_cen=X_obs[cen_indx,,drop=FALSE]

    locs_obs_n_cen=locs_obs[n_cen_indx,,drop=FALSE]
    locs_obs_cen=locs_obs[cen_indx,,drop=FALSE]


    ord=GpGp::order_maxmin(locs_obs_n_cen)

    ALL_X_obs_ordered=rbind(X_obs_n_cen[ord,],X_obs_cen)
    ALL_Y_obs_ordered=c(Y_obs_n_cen[ord],Y_obs_cen)
    ALL_locs_obs_ordered=rbind(locs_obs_n_cen[ord,,drop=FALSE],locs_obs_cen)

    if(is.null(svc_indices)){
      X_svc=matrix(1,nrow = n_obs,ncol = 1)
      n_svc=1
    }else{
      n_svc=length(svc_indices)
      X_svc <- ALL_X_obs_ordered[, svc_indices, drop = FALSE]
    }




    cov_list <- vector("list", n_svc)
    for (i in seq_len(n_svc)) {
      cov_list[[i]] <- c( sigma_sq=sigma[i],phi=phi[i] )
    }
    cov_list$nugget <- as.numeric(tau)



    NN=GpGp::find_ordered_nn(ALL_locs_obs_ordered,m = M)

    dmat <- Rfast::Dist(ALL_locs_obs_ordered)

    dmat=dmat[(n_n_cen+1):n_obs,1:n_n_cen,drop=FALSE]  # keep the distance between censored and non censored
    for(i in 1:nrow(dmat)){

      NN[i+n_n_cen,2:(M+1)]=order(dmat[i,])[1:M] # change only rows of censored data
    }


    # build Cholesky of the precision: Λ = L %*% t(L)
    # FIXED: Use Matrix::t for sparse matrices
    L <- Matrix::t(precision_decomp_sparse(cov_list, X_svc, NN, ALL_locs_obs_ordered,n_n_cen+1))

    y_cens_pred <- X_obs_cen %*% beta - forwardsolve(
      Matrix::t(L[(n_n_cen+1):n_obs, (n_n_cen+1):n_obs]),
      Matrix::t(L[1:n_n_cen, (n_n_cen+1):n_obs]) %*% (Y_obs_n_cen[ord] - X_obs_n_cen[ord,] %*% beta ))

    # FIXED: Extract diagonal efficiently
    diag_indices <- (n_n_cen+1):n_obs
    y_cens_var <- 1/sqrt(L[cbind(diag_indices, diag_indices)])


    truncated_values <- numeric(n_cen)

    for (j in 1:n_cen) {
      if (y_cens_var[j] > 0) {
        z <- (Y_obs_cen[j] - y_cens_pred[j]) / sqrt(y_cens_var[j])
        # Handle numerical issues
        if (is.finite(z)) {
          pnorm_z <- pnorm(z)
          if (pnorm_z > 1e-10) {  # Avoid division by very small numbers
            mills_ratio <- dnorm(z) / pnorm_z
            truncated_values[j] <- y_cens_pred[j] - sqrt(y_cens_var[j]) * mills_ratio
          } else {
            truncated_values[j] <- Y_obs_cen[j]  # Use censoring threshold
          }
        } else {
          truncated_values[j] <- Y_obs_cen[j]
        }
      } else {
        truncated_values[j] <- y_cens_pred[j]
      }
    }

    Y_obs[cen_indx]=truncated_values



  }

  n_pred=nrow(X_pred)
  n_obs= length(Y_obs)
  n=n_pred+n_obs

  ord=GpGp::order_maxmin(locs_obs)
  ord_pred=GpGp::order_maxmin(locs_pred)


  ALL_X_ordered=rbind(X_obs[ord,],X_pred[ord_pred,])
  Y_obs_ordered=Y_obs[ord]
  ALL_locs_ordered=rbind(locs_obs[ord,,drop=FALSE],locs_pred[ord_pred,,drop=FALSE])


  if(is.null(svc_indices)){
    X_svc=matrix(1,nrow = n,ncol = 1)
    n_svc=1
  }else{
    n_svc=length(svc_indices)
    X_svc <- ALL_X_ordered[, svc_indices, drop = FALSE]
  }




  cov_list <- vector("list", n_svc)
  for (i in seq_len(n_svc)) {
    cov_list[[i]] <- c( sigma_sq=sigma[i],phi=phi[i] )
  }
  cov_list$nugget <- as.numeric(tau)



  NN=GpGp::find_ordered_nn(ALL_locs_ordered,m = M)

  dmat <- Rfast::Dist(ALL_locs_ordered)

  dmat=dmat[(n_obs+1):n,1:n_obs,drop=FALSE]  # keep the distance between censored and non censored
  for(i in 1:nrow(dmat)){

    NN[i+n_obs,2:(M+1)]=order(dmat[i,])[1:M] # change only rows of censored data
  }


  # FIXED: Use Matrix::t for sparse matrices
  L <- Matrix::t(precision_decomp_sparse(cov_list, X_svc, NN, ALL_locs_ordered,n_obs+1))

  prediction <- X_pred[ord_pred,] %*% beta - forwardsolve(
    Matrix::t(L[(n_obs+1):n, (n_obs+1):n]),
    Matrix::t(L[1:n_obs, (n_obs+1):n]) %*% (Y_obs[ord] - X_obs[ord,] %*% beta ))


  if(pred_coef ){

    # prediction du coef i
    NN=GpGp::find_ordered_nn(locs_obs[ord,],M)
    L_i=Matrix::t(precision_decomp_sparse(cov_list, X_obs[ord,svc_indices,drop=FALSE], NN, locs_obs[ord,,drop=FALSE],1))
    coef_prediction=matrix(NA,nrow = n_pred,ncol = length(svc_indices))

    for(i in 1:length(svc_indices)){

      cov_Y_i =cov_list[[i]][1]*exp(-dmat/cov_list[[i]][2])  # 4*12
      cov_Y_i=sweep(cov_Y_i, 2,X_obs[ord,svc_indices,drop=FALSE][,i] , `*`)

      coef_prediction[ord_pred,i]= beta[svc_indices[i]]+ as.vector((cov_Y_i%*% (L_i)%*%Matrix::t(L_i) %*% (Y_obs[ord] - X_obs[ord,] %*% beta )))

    }
    ### marche pour le premier, pas pour le second, surement la covariance cov_Y_i
  }


  if(simulation){

    iid <- rnorm(n_pred)
    sim <- prediction + forwardsolve(
      Matrix::t(L[(n_obs+1):n,(n_obs+1):n]),
      iid
    )

    sim[ord_pred]=sim
    prediction[ord_pred]=prediction
    final_list=list(prediction=prediction,simulation=sim)

  }else{

    prediction[ord_pred]=prediction  # inverse de ord_pred

    final_list=list(prediction=prediction)
  }

  if(pred_coef){
    final_list <- append(final_list, list(coef_prediction))

  }

  return(final_list)
}




# ============================================================================
# ESTIMATION  HELPER
# ============================================================================


fit_censored_freq_Vecchia <- function(Y_obs, locs_obs, X_obs, M, svc_indices,
                                      censored_indices, fixed_range ) {

  n <- length(Y_obs)
  censored_indices_pos <- which(censored_indices == 1)
  non_censored_indices <- setdiff(1:n, censored_indices_pos)

  locs_cen <- locs_obs[censored_indices_pos, , drop = FALSE]
  locs_n_cen <- locs_obs[non_censored_indices, , drop = FALSE]
  Y_cen <- Y_obs[censored_indices_pos]
  Y_n_cen <- Y_obs[non_censored_indices]
  X_cen <- X_obs[censored_indices_pos, , drop = FALSE]
  X_n_cen <- X_obs[non_censored_indices, , drop = FALSE]

  ord <- GpGp::order_maxmin(locs_n_cen)
  locs_n_cen <- locs_n_cen[ord,  , drop = FALSE]
  Y_n_cen <- Y_n_cen[ord]
  X_n_cen <- X_n_cen[ord,  , drop = FALSE]

  all_locs <- rbind(locs_n_cen, locs_cen)
  all_X <- rbind(X_n_cen, X_cen)
  all_Y <- c(Y_n_cen, Y_cen)

  ## NN modified
  M <- min(M, n - 1)
  NN <- GpGp::find_ordered_nn(all_locs, M)

  dmat <- Rfast::Dist(all_locs)
  phi0 <- 0.5 * mean(dmat[lower.tri(dmat)])

  if (length(non_censored_indices) != n) {
    dmat <- dmat[(length(Y_n_cen) + 1):n, 1:length(Y_n_cen), drop = FALSE]
    for (i in 1:nrow(dmat)) {
      NN[i + length(non_censored_indices), 2:(M + 1)] <- order(dmat[i, ])[1:M]
    }
  }
  lm_fit <- lm(Y_obs ~  as.matrix(X_obs) - 1)
  mean_start <- as.vector(lm_fit$coefficients)

  p <- ncol(X_obs)
  sigma0 <- var(Y_obs) / (20 * p)

  if (is.null(svc_indices)) {
    n_svc <- 1
  } else {
    n_svc <- length(svc_indices)
  }

  # Adapt parameter structure based on fixed_range
  if (fixed_range && n_svc > 1) {
    # Structure: [sigma1, sigma2, ..., sigmaN, phi_shared, nugget]
    cov_start <- c(rep(sigma0, n_svc), phi0, sigma0 * 1e-2)
    lower <- c(rep(-Inf, length(mean_start)), rep(1e-6, length(cov_start)))
    upper <- c(rep(Inf, length(mean_start)), rep(Inf, length(cov_start)))
  } else {
    # Standard structure: [sigma1, phi1, sigma2, phi2, ..., nugget]
    cov_start <- c(rep(c(sigma0, phi0), n_svc), sigma0 * 1e-2)
    lower <- c(rep(-Inf, length(mean_start)), rep(1e-6, length(cov_start)))
    upper <- c(rep(Inf, length(mean_start)), rep(Inf, length(cov_start)))
  }

  par_start <- c(mean_start, cov_start)

  # Negative log-likelihood to minimize
  objective <- function(par) {
    mean_par <- par[1:length(mean_start)]
    cov_par <- par[(length(mean_start) + 1):length(par)]

    # Construct list of covariance blocks
    cov_list <- vector("list", n_svc + 1)

    # In the objective function:
    if (fixed_range && n_svc > 1) {
      phi_shared <- cov_par[n_svc + 1]

      for (i in seq_len(n_svc)) {
        cov_list[[i]] <- c(sigma_sq = cov_par[i], phi = phi_shared)
      }
      cov_list$nugget <- cov_par[length(cov_par)]  # Remove the c(nugget = ...)
    } else {
      for (i in seq_len(n_svc)) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        cov_list[[i]] <- c(sigma_sq = cov_par[idx[1]], phi = cov_par[idx[2]])
      }
      cov_list$nugget <- cov_par[length(cov_par)]  # Remove the c(nugget = ...)
    }

    ll <- Like_SVC_censored_freq_Vecchia(mean_par, cov_list, all_Y, all_locs, all_X,
                                         NN, svc_indices, n_no_cen = length(non_censored_indices))
    return(-ll)
  }

  objective(par = par_start)
  opt <- optim(par_start, objective, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(trace = 1, factr = 1e7))

  # Reconstruct parameters in standard format for output
  if (fixed_range && n_svc > 1) {
    mean_par <- opt$par[1:length(mean_start)]
    cov_par <- opt$par[(length(mean_start) + 1):length(opt$par)]

    phi_shared <- cov_par[n_svc + 1]
    nugget <- cov_par[length(cov_par)]

    # Reconstruct as [sigma1, phi1, sigma2, phi2, ..., nugget]
    reconstructed_par <- numeric(length(mean_start) + 2 * n_svc + 1)
    reconstructed_par[1:length(mean_start)] <- mean_par

    for (i in seq_len(n_svc)) {
      idx <- length(mean_start) + (2 * (i - 1) + 1):(2 * i)
      reconstructed_par[idx[1]] <- cov_par[i]        # sigma
      reconstructed_par[idx[2]] <- phi_shared         # phi (shared)
    }
    reconstructed_par[length(reconstructed_par)] <- nugget

    opt$par <- reconstructed_par
  }

  return(list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message
  ))
}

Like_SVC_censored_freq_Vecchia <- function(mean_parameters, cov_parameters, Y_obs,
                                           locs_obs, X_obs, NN, svc_indices, n_no_cen) {
  n <- length(Y_obs)

  if (is.null(svc_indices)) {
    X_svc <- matrix(1, nrow = n, ncol = 1)
  } else {
    X_svc <- X_obs[, svc_indices, drop = FALSE]
  }

  L <- precision_decomp_sparse(cov_parameters, X_svc, NN, locs_obs)
  mu <- X_obs %*% mean_parameters
  resid <- drop(Y_obs - mu)

  # Non-censored likelihood
  # FIXED: Extract diagonal more efficiently
  L_sub <- L[1:n_no_cen, 1:n_no_cen]
  z <- L_sub %*% resid[1:n_no_cen]
  quad <- sum(z * z)

  # Extract diagonal values directly
  diag_indices <- 1:n_no_cen
  diag_vals <- L[cbind(diag_indices, diag_indices)]
  logdet_prec <- 2 * sum(log(diag_vals))

  c0 <- n_no_cen * log(2 * pi)
  ll_no_censored <- 0.5 * logdet_prec - 0.5 * (c0 +quad)

  # Censored likelihood
  # FIXED: Use Matrix::t for sparse matrices
  L <- Matrix::t(L)

  cond_exp <- mu[(n_no_cen + 1):n] -
    forwardsolve(
      Matrix::t(L[(n_no_cen + 1):n, (n_no_cen + 1):n]),
      Matrix::t(L[1:n_no_cen, (n_no_cen + 1):n]) %*% (resid[1:n_no_cen])
    )

  # FIXED: Extract diagonal efficiently for conditional variance too
  diag_indices_cen <- (n_no_cen + 1):n
  cond_diag_vals <- L[cbind(diag_indices_cen, diag_indices_cen)]
  cond_sd <- 1 / cond_diag_vals

  ll_cen <- pnorm(q = Y_obs[(n_no_cen + 1):n], mean = cond_exp, sd = cond_sd, log.p = TRUE)

  ll <- ll_no_censored + sum(ll_cen)
  return(ll)
}
