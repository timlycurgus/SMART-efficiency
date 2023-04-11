### This script will analyze a SMART for the scenario with estimated weights

analyze_SMART2 <- function(d,
                           design = 2,
                           rounding = "up",
                           conservative = TRUE,
                           corstr = c("identity", "exchangeable", "ar1"),
                           L = NULL,
                           corstr.estimate = corstr,
                           pool.time = TRUE,
                           pool.dtr = TRUE,
                           niter = 5000,
                           tol = 1e-8,
                           maxiter.solver = 1000,
                           wts_meth,
                           ...
){
  
  wts = calc_weights(d$data, long = FALSE, method = wts_meth)
  d$data$weight = wts[[1]]
  S = cbind(wts[[2]], wts[[3]])
  
  d1 <- reshape(d$data, varying = list(grep("Y", names(d$data))),
                ids = d$data$id, 
                times = times, direction = "long", v.names = "Y")
  d1 <- d1[order(d1$id, d1$time), ]
  d1 <- d1[(d1$time == 2),]
  
  param.hat <- try(estimate.params0(d1, diag(1), design, rep(0, 4),
                                    maxiter.solver, tol))
  
  sigma2.hat <- estimate.sigma2_0(d1, design, param.hat,
                                  pool.time = TRUE, pool.dtr = TRUE)
  
  condVars <- lapply(split.SMART(d1), function(x) {
    var(subset(x, select = grep("Y", names(x), value = TRUE)))
  })
  
  # Iterate parameter estimation
  outcome.var <- varmat_stand(sigma2.hat)
  
  param.hat <- estimate.params0(d1, outcome.var, design, param.hat,
                                maxiter.solver, tol)
  
  # Iterate until estimates of gammas and rho converge
  for (j in 1:maxiter.solver) {
    sigma2.new <- estimate.sigma2_0(d1, design, param.hat,
                                    pool.time = TRUE, pool.dtr = TRUE)
    
    outcomeVar.new <- varmat_stand(sigma2.new)
    
    param.new <- estimate.params0(d1, outcomeVar.new, design, param.hat,
                                  maxiter.solver, tol)
    
    if (norm(param.new - param.hat, type = "F") <= tol &
        norm(as.matrix(sigma2.new) - as.matrix(sigma2.hat), 
             type = "F") <= tol) {
      param.hat <- param.new
      sigma2.hat <- sigma2.new
      iter <- j
      break
    } else {
      param.hat <- param.new
      sigma2.hat <- sigma2.new
    }
  }
  
  param.var <- estimate.paramvar_2(d1, varmat_stand(sigma2.hat), design,
                                   gammas = param.hat, S)
  
  confLB <- 
    L %*% param.hat - sqrt(L %*% param.var %*% L) * qnorm(.975)
  confUB <-
    L %*% param.hat + sqrt(L %*% param.var %*% L) * qnorm(.975)
  
  CI = c(confLB, confUB)
  
  return(list(CI, param.hat, param.var))
}

########################### Unique Functions Needed ############################
### Note: these are mostly functions to compute weights and score functions. 

calc_weights <- function(d, long, method){
  
  if(method == 'estimated'){
    multiplier.A1 <- mean(d$A1 == 1)
    multiplier.A2.1 <- sum(d$A2NR == 1 & d$A1 == 1 & d$R == 0) / sum(d$A1 == 1 & d$A2NR != 0) 
    multiplier.A2.0 <- sum(d$A2NR == 1 & d$A1 == -1 & d$R == 0) / sum(d$A1 == -1 & d$A2NR != 0)  
    
    d$weight <- (1 / multiplier.A1) * (d$A1 == 1) + 
      (1 / (1 - multiplier.A1)) * (d$A1 == -1)
    
    d$weight[d$R == 0 & d$A1 == 1] <- d$weight[d$R == 0 & d$A1 == 1] *  
      ((1 / multiplier.A2.1) * (d$A2NR[d$R == 0 & d$A1 == 1] == 1) + 
         (1 / (1 - multiplier.A2.1)) * (d$A2NR[d$R == 0 & d$A1 == 1] == -1))
    d$weight[d$R == 0 & d$A1 == -1] <- d$weight[d$R == 0 & d$A1 == -1] *  
      ((1 / multiplier.A2.0) * (d$A2NR[d$R == 0 & d$A1 == -1] == 1) + 
         (1 / (1 - multiplier.A2.0)) * (d$A2NR[d$R == 0 & d$A1 == -1] == -1))
    scores = calc_score(d, long, method, multiplier.A1, 
                        multiplier.A2.1, multiplier.A2.0)
    
  } else {
    
    d$binA1 <- as.numeric(d$A1 == 1)
    d$binA2 <- as.numeric(d$A2NR == 1)
    
    if(long == TRUE){
      mod_stage1 <- glm(binA1 ~ Y0 + X, family = binomial, data = d)
      multiplier.A1 <- predict(mod_stage1, type = "response")
      mod_stage2 <- glm(binA2 ~ Y0 + 
                          A1 + Y1 + L, 
                        data = subset(d, R == 0),
                        family = binomial)
      multiplier.A2 <- predict(mod_stage2, type = "response")
    } else{
      mod_stage1 <- glm(binA1 ~ 1 + X, family = binomial, data = d)
      multiplier.A1 <- predict(mod_stage1, type = "response")
      mod_stage2 <- glm(binA2 ~ A1 + L, 
                        data = subset(d, R == 0),
                        family = binomial)
      multiplier.A2 <- predict(mod_stage2, type = "response")
    }
    d$weight <- (1 / multiplier.A1) * (d$A1 == 1) + 
      (1 / (1 - multiplier.A1)) * (d$A1 == -1)
    
    d$weight[d$R == 0] <- d$weight[d$R == 0] *  ((1 / multiplier.A2) * (d$A2NR[d$R == 0] == 1) + 
                                                   (1 / (1 - multiplier.A2)) * (d$A2NR[d$R == 0] == -1))
    
    scores = calc_score(d, long, method, multiplier.A1, 
                        multiplier.A2, multiplier.A2)
  }
  
  sa1 = scores[[1]]
  sa2 = scores[[2]]
  
  return(list(d$weight, sa1, sa2))
}

calc_score <- function(d, long, method, ma1, ma21, ma20){
  
  if(method == 'estimated'){
    score_alpha1 = as.numeric(d$A1 == 1) - ma1
    score_alpha2 = as.numeric(d$A2NR == 1) - ma21*(d$A1 == 1 & d$R == 0) -
      ma20*(d$A1 == -1 & d$R == 0)
  } else{
    
    d$binA1 <- as.numeric(d$A1 == 1)
    d$binA2 <- as.numeric(d$A2NR == 1)
    
    if(long == TRUE){
      score_alpha1 = cbind(1, d[,c('Y0','X')]) *  
        matrix(rep(d$binA1 - ma1, 3), nrow = nrow(d))
      
      d$dev = 0
      d$dev[d$R == 0] = d$binA2[d$R == 0] - ma21
      score_alpha2 = cbind(1, d[,c('Y0', 'Y1','L')]) * 
        matrix(rep(d$dev,4), nrow = nrow(d))
    } else{
      score_alpha1 = cbind(1, d[,c('X')]) *  
        matrix(rep(d$binA1 - ma1, 2), nrow = nrow(d))
      
      d$dev = 0
      d$dev[d$R == 0] = d$binA2[d$R == 0] - ma21
      score_alpha2 = cbind(1, d[,c('L')]) * 
        matrix(rep(d$dev,2), nrow = nrow(d))
    }
  }
  return(list(score_alpha1, score_alpha2))
}

calc_score_new <- function(d, long, method, ma1, ma21, ma20){
  
  if(method == 'estimated'){
    d$binA1 <- as.numeric(d$A1 == 1)
    d$binA2 <- as.numeric(d$A2NR == 1)
    
    score_alpha1 = matrix(d$binA1 - ma1, nrow = nrow(d))
    
    d$dev = 0
    d$dev[d$R == 0] = d$binA2[d$R == 0] - ma21
    score_alpha2 = matrix(d$dev, nrow = nrow(d))
    
    score_alpha1 = as.numeric(d$A1 == 1) - ma1
    score_alpha2 = as.numeric(d$A2NR == 1) - ma21*(d$A1 == 1 & d$R == 0) -
      ma20*(d$A1 == -1 & d$R == 0)
  } else{
    
    d$binA1 <- as.numeric(d$A1 == 1)
    d$binA2 <- as.numeric(d$A2NR == 1)
    
    if(long == TRUE){
      score_alpha1 = cbind(1, d[,c('Y0','X')]) *  
        matrix(rep(d$binA1 - ma1, 3), nrow = nrow(d))
      
      d$dev = 0
      d$dev[d$R == 0] = d$binA2[d$R == 0] - ma21
      score_alpha2 = cbind(1, d[,c('Y0', 'Y1','L')]) * 
        matrix(rep(d$dev,4), nrow = nrow(d))
    } else{
      score_alpha1 = cbind(1, d[,c('X')]) *  
        matrix(rep(d$binA1 - ma1, 2), nrow = nrow(d))
      
      d$dev = 0
      d$dev[d$R == 0] = d$binA2[d$R == 0] - ma21
      score_alpha2 = cbind(1, d[,c('L')]) * 
        matrix(rep(d$dev,2), nrow = nrow(d))
    }
  }
  return(list(score_alpha1, score_alpha2))
}

estimate.paramvar_2 <- function(d, V, design, gammas, S) {
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  n <- length(unique(d$id))
  J <- esteqn.jacobian0(d, V, design)
  
  solve(-J) %*%  meat.compute_2(d, V, design, gammas, S) %*%
    solve(-J) / n
}

meat.compute_2 <- function(d, V, design, gammas, S) {
  
  # d should be LONG
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- matrix(meanvec_stand(design, gammas)[3,], nrow = 1)
  dmat <- mod.derivs_stand(design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE) 
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR), 
                  X = d$X)
  
  # sum over IDs
  mmt <- Reduce("+", lapply(split.data.frame(resids, resids$id), function(obs) {
    # compute weighted residual matrix for individual
    # resid <- obs$weight * obs[, grepl("dtr", names(obs))] * 
    # (obs$Y - meanvec(times, spltime, design, gammas))
    # Sum over DTRs
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      matrix(c(dmat[[dtr]]), ncol = 1) %*% solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% t(m)
  })) / n
  
  resids1 = cbind(resids, S)
  nvars = ncol(S)
  
  mst <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      matrix(c(dmat[[dtr]]), ncol = 1) %*% solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])
  })) / n
  
  smt <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      matrix(c(dmat[[dtr]]), ncol = 1) %*% solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    t(as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])) %*%
      t(m)
  })) / n
  
  sst <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    t(as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])) %*% 
      as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])
  })) / n
  
  mmt - mst%*%solve(sst)%*%smt
  
}


########################### Generic Functions Needed ###########################
### Note: there are more because this solely differs from Scenario 0 through the
### standard error estimation. 

estimate.params0 <- function(d, V, design, 
                             start, maxiter.solver, tol) {
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i)
      V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  params.hat <- start
  j <- 0
  epsilon <- 10
  J <- esteqn.jacobian0(d, V, design)
  
  for (j in 1:maxiter.solver) {
    params.hat.new <-
      params.hat - solve(J) %*% esteqn.compute_0(d, V, design, 
                                                 gammas = params.hat)
    epsilon <- norm(params.hat.new - params.hat, type = "F")
    params.hat <- params.hat.new
    if (epsilon <= tol)
      break
  }
  params.hat
}

esteqn.jacobian0 <- function(d, V, design) {
  
  nDTR <- switch(design, 8, 4, 3)
  deriv <- mod.derivs_stand(design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR) {
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  }
  
  Reduce("+", lapply(split.data.frame(d, d$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      unique(x[[paste0("dtr", dtr)]]) * unique(x$weight) *
        crossprod(deriv[[dtr]],solve(V[[dtr]])) %*% deriv[[dtr]]
    }))
  })) / -length(unique(d$id))
}

esteqn.compute_0 <- function(d, V, design, gammas) {
  # d is an unreplicated data set in LONG format
  # V is a working covariance matrix
  
  n <- length(unique(d$id))
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  mvec <- matrix(meanvec_stand(design, gammas)[3,], nrow = 1)
  deriv <- mod.derivs_stand(design)
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR))
  
  DVinv <- lapply(1:nDTR, function(i) crossprod(deriv[[i]], solve(V[[i]])))
  
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      DVinv[[dtr]] %*% matrix(x[[paste0("dtr", dtr)]], ncol = 1)
    }))
  })
  ) / n
}

estimate.sigma2_0 <- function(d, design, gammas,
                              pool.time = FALSE, pool.dtr = FALSE) {
  
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- matrix(meanvec_stand(design, gammas)[3,], nrow = 1)
  Dmat <- mod.derivs_stand(design)
  
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format")
  
  resids <- (matrix(rep(d$Y, nDTR), ncol = nDTR) -
               matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE))^2
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR))
  
  weightmat <- cbind("id" = d$id, "time" = d$time,
                     d$weight * d[, grep("dtr", names(d))])
  
  if (pool.time & pool.dtr) {
    numerator <-
      sum(subset(resids, select = grep("dtr", names(resids))))
    denominator <-
      (sum(subset(weightmat, select = grep(
        "dtr", names(weightmat)
      ))) - length(gammas))
  } else if (pool.time & !pool.dtr) {
    numerator <-
      apply(subset(resids, select = grep("dtr", names(resids))), 2, sum)
    denominator <-
      apply(subset(weightmat, select = grep("dtr", names(weightmat))), 2, sum) -
      length(gammas)
  } else if (!pool.time & pool.dtr) {
    numerator <- do.call("c", lapply(times, function(x) {
      sum(subset(resids, time == x, select = grep("dtr", names(resids))))
    }))
    denominator <- do.call("c", lapply(times, function(x) {
      sum(subset(weightmat, time == x, select = grep("dtr", names(weightmat))))
    }))
    names(numerator) <- names(denominator) <- paste0("time", times)
  } else {
    numerator <- Reduce(function(...)
      merge(..., by = "time"),
      lapply(1:nDTR, function(dtr) {
        x <- list(resids[[paste0("dtr", dtr)]])
        sumsqrs <-
          aggregate(x = setNames(x, paste0("dtr", dtr)),
                    by = list("time" = resids$time),
                    sum)
        x <- list(weightmat[[paste0("dtr", dtr)]])
        sumwts <-
          aggregate(x = setNames(x, paste0("dtr", dtr)),
                    by = list("time" = resids$time),
                    sum)
        sumsqrs[, 2] <- sumsqrs[, 2] / sumwts[, 2]
        sumsqrs
      }))
    rownames(numerator) <- numerator$time
    numerator <- as.matrix(numerator[,-1])
    denominator <- 1
  }
  as.matrix(numerator / denominator)
}

mod.derivs_stand <- function(design) {
  
  nDTR <- switch(design, 8, 4, 3)
  
  if(design == 2){
    dtr11 = matrix(c(1, 1, 1, 1), nrow = 1)
    dtr10 = matrix(c(1, 1, -1, -1), nrow = 1)
    dtr01 = matrix(c(1, -1, 1, -1), nrow = 1)
    dtr00 = matrix(c(1, -1, -1, 1), nrow = 1)
    dmat = list(dtr11, dtr10, dtr01, dtr00)
  }
  return(dmat)
}

meanvec_stand <- function(design, gammas) {
  
  A <- dtrIndex(design)
  do.call(cbind, lapply(1:length(A$a1), function(dtr) {
    do.call(rbind, lapply(times, function(obsTime) {
      marginal.model_stand(A$a1[dtr], A$a2R[dtr], A$a2NR[dtr], 
                           design, gammas)
    }))
  }))
}

marginal.model_stand <- function(a1, a2R, a2NR, design, gammas) {
  
  if(design == 2){
    mu <- gammas[1] + gammas[2]*a1 + gammas[3]*a2NR + gammas[4]*a1*a2NR
  }
  mu
}

varmat_stand <- function(sigma2) {
  
  nDTR <- switch(design, 8, 4, 3)
  
  lapply(1:nDTR, function(dtr) {
    diag(sqrt(sigma2)) %*% diag(sqrt(sigma2))
  })
}

