### This script will analyze a SMART using all four techniques at once.

analyze_SMARTfull <- function(d, 
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
  
  if(wts_meth != 'known'){
    
    wts = calc_weights(d$data, long = TRUE, method = wts_meth)
    d$data$weight = wts[[1]]
    S = data.frame(cbind(wts[[2]], wts[[3]]))
    S$nrep = 3
    S <- as.data.frame(lapply(S, rep, S$nrep))
    #colnames(S) = c('S1', 'SY0', 'SX', 'S11', 'S01', 'S1', 'SL')
    
    d1 <- reshape(d$data, varying = list(grep("Y", names(d$data))),
                  ids = d$data$id, 
                  times = times, direction = "long", v.names = "Y")
    d1 <- d1[order(d1$id, d1$time), ]
    #S <- d1[ , grepl( "S" , names( d1 ) ) ]
    
  } else{
    d1 <- reshape(d$data, varying = list(grep("Y", names(d$data))),
                  ids = d$data$id, 
                  times = times, direction = "long", v.names = "Y")
    d1 <- d1[order(d1$id, d1$time), ]
    S = matrix(rep(0, 2*nrow(d1)), ncol = 2)
  }
  
  param.hat <- try(estimate.paramsfull(d1,
                                       diag(rep(1, length(times))),
                                       times, spltime,
                                       design, rep(0, 8),
                                       maxiter.solver, tol))
  
  sigma2.hat <- estimate.sigma2_full(d1, times, spltime, design, param.hat,
                                     pool.time = FALSE, pool.dtr = FALSE)
  
  rho.hat <- estimate.rho_full(d1, times, spltime, design, sqrt(sigma2.hat),
                               param.hat, corstr = corstr.estimate,
                               pool.dtr = FALSE)
  
  # Compute variance matrices for all conditional cells
  condVars <- lapply(split.SMART(d$data), function(x) {
    var(subset(x, select = grep("Y", names(x), value = TRUE)))
  })
  
  # Iterate parameter estimation
  outcome.var <- varmat_adj(sigma2.hat, rho.hat, times, design,
                            corstr = corstr.estimate)
  
  param.hat <- estimate.paramsfull(d1, outcome.var, times, spltime,
                                   design, param.hat,
                                   maxiter.solver, tol)
  
  # Iterate until estimates of gammas and rho converge
  for (j in 1:maxiter.solver) {
    sigma2.new <- estimate.sigma2_full(d1, times, spltime, design, param.hat,
                                       pool.time = FALSE, pool.dtr = FALSE)
    rho.new <- estimate.rho_full(d1, times, spltime, design, sqrt(sigma2.hat),
                                 param.hat, corstr = corstr.estimate, 
                                 pool.dtr = FALSE)
    
    outcomeVar.new <- varmat_adj(sigma2.new, rho.new, times, design, 
                                 corstr = corstr.estimate)
    
    param.new <- estimate.paramsfull(d1, outcomeVar.new, times, spltime,
                                     design, start = param.hat,
                                     maxiter.solver, tol)
    
    if (norm(param.new - param.hat, type = "F") <= tol &
        norm(as.matrix(sigma2.new) - as.matrix(sigma2.hat), 
             type = "F") <= tol &
        norm(as.matrix(rho.new) - as.matrix(rho.hat),
             type = "F") <= tol) {
      param.hat <- param.new
      sigma2.hat <- sigma2.new
      rho.hat <- rho.new
      iter <- j
      break
    } else {
      param.hat <- param.new
      sigma2.hat <- sigma2.new
      rho.hat <- rho.new
    }
  }
  
  param.var <- estimate.paramvar_full(d1,
                                      varmat_adj(sigma2.hat, rho.hat, times, 
                                                 design, corstr.estimate),
                                      times, spltime, design,
                                      gammas = param.hat, S)
  
  confLB <- 
    L %*% param.hat - sqrt(L %*% param.var %*% L) * qnorm(.975)
  confUB <-
    L %*% param.hat + sqrt(L %*% param.var %*% L) * qnorm(.975)
  
  CI = c(confLB, confUB)
  
  return(list(CI, param.hat, param.var, rho.hat))
}

########################### Unique Functions Needed ############################
### Note: these functions largely combine the functions from Scenario 3 with
### that from Scenario 0 and 2.

estimate.paramsfull <- function(d, V, times, spltime, design,
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
    J <- esteqn.jacobianfull(d, V, times, spltime, design)
    
    for (j in 1:maxiter.solver) {
      params.hat.new <-
        params.hat - solve(J) %*% esteqn.compute_full(d, V, times, spltime, 
                                                      design, gammas = params.hat)
      epsilon <- norm(params.hat.new - params.hat, type = "F")
      params.hat <- params.hat.new
      if (epsilon <= tol)
        break
    }
    params.hat
  }

esteqn.jacobianfull <- function(d, V, times, spltime, design) {
  
  nDTR <- switch(design, 8, 4, 3)
  deriv <- mod.derivs(times, spltime, design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR) {
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  }
  
  Reduce("+", lapply(split.data.frame(d, d$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      unique(x[[paste0("dtr", dtr)]]) * unique(x$weight) *
        crossprod(matrix(c(deriv[[dtr]], x$X), nrow = length(times)), 
                  solve(V[[dtr]])) %*% 
        matrix(c(deriv[[dtr]], x$X), nrow = length(times))
    }))
  })) / -length(unique(d$id))
}

esteqn.compute_full <- function(d, V, times, spltime, design, gammas) {
  # d is an unreplicated data set in LONG format
  # V is a working covariance matrix
  
  n <- length(unique(d$id))
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  mvec <- meanvec_full(times, spltime, design, gammas)
  deriv <- mod.derivs(times, spltime, design)
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    (matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE) +  
       gammas[8]*matrix(rep(d$X, nDTR), ncol = nDTR))
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR), X = d$X)
  
  DVinv <- lapply(1:nDTR, function(i) crossprod(deriv[[i]], solve(V[[i]])))
  
  Reduce("+", lapply(split.data.frame(resids, resids$id), function(x) {
    Reduce("+", lapply(1:nDTR, function(dtr) {
      crossprod(matrix(c(deriv[[dtr]], rep(x$X,times)), nrow = length(times)), 
                solve(V[[dtr]])) %*%
        matrix(x[[paste0("dtr", dtr)]], ncol = 1)
      #matrix(c(DVinv[[dtr]],x$X*solve(V[[dtr]])), ncol = 1) %*% matrix(x[[paste0("dtr", dtr)]], ncol = 1)
    }))
  })
  ) / n
}

meanvec_full <- function(times, spltime, design, gammas) {
  
  A <- dtrIndex(design)
  do.call(cbind, lapply(1:length(A$a1), function(dtr) {
    do.call(rbind, lapply(times, function(obsTime) {
      marginal.model_full(A$a1[dtr], A$a2R[dtr], A$a2NR[dtr], obsTime,
                          spltime, design, gammas)
    }))
  }))
}

marginal.model_full <- function(a1, a2R, a2NR, obsTime, spltime, design, gammas) {
  if (obsTime <= spltime) {
    mu <- gammas[1] + obsTime * (gammas[2] + gammas[3] * a1)
  } else {
    if (design == 1) {
      if (length(gammas) != 9)
        stop("for design 1 gammas must be length 9.")
      mu <-
        gammas[4] + gammas[5] * a1 + gammas[6] * a2R + gammas[7] * a2NR +
        gammas[8] * a1 * a2R + gammas[9] * a1 * a2NR
    } else if (design == 2) {
      mu <-
        gammas[4] + gammas[5] * a1 + gammas[6] * a2NR + gammas[7] * a1 * a2NR
    } else if (design == 3) {
      if (length(gammas) != 6)
        stop("for design 3 gammas must be length 6.")
      mu <-
        gammas[4] + gammas[5] * a1 + gammas[6] * (a1 == 1) * a2NR
    }
    mu <- gammas[1] + spltime * (gammas[2] + gammas[3] * a1) +
      (obsTime - spltime) * mu
  }
  mu
}

estimate.sigma2_full <- function(d, times, spltime, design, gammas,
                                 pool.time = FALSE, pool.dtr = FALSE) {
    
    n <- length(unique(d$id))
    nDTR <- switch(design, 8, 4, 3)
    mvec <- meanvec_full(times, spltime, design, gammas)
    Dmat <- mod.derivs(times, spltime, design)
    
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

estimate.rho_full <- function(d, times, spltime, design, sigma, gammas,
                              corstr = c("exchangeable", "ar1", "unstructured"),
                              pool.dtr = TRUE) {
  ## FIXME: Allow for different rhos across DTRs
  
  corstr <- match.arg(corstr)
  if (any(is(d$Y) == "NULL")) stop("d has to be in long format.")
  
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec_full(times, spltime, design, gammas)
  Dmat <- mod.derivs(times, spltime, design)
  
  sigma <- reshapeSigma(sigma, times, design)
  
  # if (sum(rownames(sigma) == times) == length(times) & ncol(sigma) != nDTR) {
  #   sigma <- matrix(rep(sigma, nDTR), ncol = nDTR)
  # } else if (sum(grepl("dtr", rownames(sigma))) != 0) {
  #   sigma <- t(matrix(rep(sigma, length(times)), ncol = length(times)))
  # } else if (nrow(sigma) == ncol(sigma) & ncol(sigma) == 1) {
  #   sigma <- matrix(rep(sigma, nDTR * length(times)), ncol = nDTR)
  # } else if (nrow(sigma) == length(times) & ncol(sigma) == 1) {
  #   sigma <- matrix(rep(sigma, nDTR), ncol = nDTR)
  # }
  
  colnames(sigma) <- paste0("dtr", 1:nDTR)
  rownames(sigma) <- paste0("time", times)
  
  # Compute residuals (Y_{it} - mu_{t}(a1, a2))
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE)
  
  # Create data.frame from resids, indexed by time and ID
  # Multiply residuals by DTR indicators (dtr1[i] = 1 iff i is consistent with
  # DTR 1)
  resids <- cbind("id" = d$id, "time" = d$time,
                  resids * d[, grep("dtr", names(d))])
  
  weights <- aggregate(weight ~ id, d, unique)
  # Create matrix of weights per person-time per DTR
  weightmat <- cbind("id" = d$id, "time" = d$time, 
                     d$weight * d[, grep("dtr", names(d))])
  
  # Sum weights over individuals 
  # only use one weight per person -- weightmat has duplicated rows: 1 per time
  sumweights.time <- Reduce(function(...)
    merge(..., by = "time"),
    lapply(1:sum(grepl("dtr", names(
      d
    ))), function(dtr) {
      x <- list(weightmat[[paste0("dtr", dtr)]])
      sumwts <-
        aggregate(x = setNames(x, paste0("dtr", dtr)),
                  by = list("time" = resids$time),
                  sum)
    }))[,-1]
  rownames(sumweights.time) <- times
  sumweights <- apply(sumweights.time, 2, unique)
  
  if (corstr == "exchangeable") {
    # For every id and every dtr, compute 
    # (\sum_{s<t} r_{s} * r_{t} / sigma_{s} sigma_{t}
    r <-
      do.call(
        rbind,
        lapply(split.data.frame(resids, resids$id),
               function(subdat) {
                 sapply(1:sum(grepl("dtr", names(subdat))),
                        function(dtr) {
                          sum(sapply(2:length(subdat$time), function(ti) {
                            sum((subdat[ti, paste0("dtr", dtr)] /
                                   sigma[ti, dtr]) *
                                  (subdat[1:(ti - 1), paste0("dtr", dtr)] /
                                     sigma[1:(ti - 1), dtr]))
                          }))
                        }) * weights$weight[weights$id == unique(subdat$id)]
               }))
    colnames(r) <- grep("dtr", names(resids), value = T)
    denominator <- sumweights * (length(times) * (length(times) - 1) / 2) - 
      length(gammas)
  } else if (corstr == "ar1") {
    r <-
      do.call(
        rbind,
        lapply(split.data.frame(resids, resids$id),
               function(subdat) {
                 sapply(1:sum(grepl("dtr", names(subdat))), function(dtr) {
                   sum(sapply(2:length(subdat$time), function(time) {
                     sum((subdat[time, paste0("dtr", dtr)] / sigma[time, dtr]) *
                           (subdat[time - 1, paste0("dtr", dtr)] /
                              sigma[time - 1, dtr]))
                   }))
                 }) * weights$weight[weights$id == unique(subdat$id)]
               }))
    colnames(r) <- grep("dtr", names(resids), value = T)
    denominator <- sumweights * (length(times) - 1) - length(gammas)
  } else if (corstr == "unstructured") {
    m <- combn(times, 2)
    r <- 
      lapply(
        split.data.frame(resids, resids$id),
        function(subdat) {
          x <- do.call(rbind, lapply(1:ncol(m), function(tpair) {
            # data.frame("tpair" = paste(m[, tpair], collapse = ""),
            sapply(1:sum(grepl("dtr", names(subdat))), function(dtr) {
              prod(subdat[subdat$time %in% m[, tpair], paste0("dtr", dtr)]) /
                prod(sigma[rownames(sigma) %in% m[, tpair], paste0("dtr", dtr)])
            }) * weights$weight[weights$id == unique(subdat$id)]
          }))
          rownames(x) <-
            sapply(1:ncol(m), function(j)
              paste(m[, j], collapse = ""))
          x
        })
    r <- lapply(1:ncol(m), function(tpair) {
      do.call(rbind,
              lapply(r, function(x) x[paste(m[, tpair], collapse = ""), ]))
    })
    names(r) <- sapply(1:ncol(m), function(j) paste(m[, j], collapse = ""))
  }
  
  if (is.list(r)) {
    numerator <- lapply(r, colSums)
    denominator <- sumweights - length(gammas)
    rhoHat <- do.call(rbind, lapply(numerator, function(x) x / denominator))
    if (pool.dtr) 
      rhoHat <- apply(rhoHat, 1, mean)
  } else {
    numerator <- colSums(r)
    rhoHat <- numerator / denominator
    if (pool.dtr)
      rhoHat <- mean(rhoHat)
  }
  
  return(rhoHat)
}

estimate.paramvar_full <- function(d, V, times, spltime, design, gammas, S) {
  
  nDTR <- switch(design, 8, 4, 3)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  n <- length(unique(d$id))
  J <- esteqn.jacobianfull(d, V, times, spltime, design)
  solve(-J) %*%  meat.compute_full(d, V, times, spltime, design, gammas, S) %*%
    solve(-J) / n
}

meat.compute_full <- function(d, V, times, spltime, design, gammas, S) {
  # d should be LONG
  n <- length(unique(d$id))
  nDTR <- switch(design, 8, 4, 3)
  mvec <- meanvec_full(times, spltime, design, gammas)
  dmat <- mod.derivs(times, spltime, design)
  
  if (!is.list(V) | (is.list(V) & length(V) == 1)) {
    V <- lapply(1:nDTR, function(i) V)
  } else if (length(V) != nDTR)
    stop(paste("V must be a single matrix or a list of", nDTR, "matrices."))
  
  resids <- matrix(rep(d$Y, nDTR), ncol = nDTR) -
    (matrix(rep(t(mvec), n), ncol = ncol(mvec), byrow = TRUE) +  
       gammas[8]*matrix(rep(d$X, nDTR), ncol = nDTR))
  resids <- cbind("id" = d$id, resids * d[, grep("dtr", names(d))] *
                    matrix(rep(d$weight, nDTR), ncol = nDTR), X = d$X)
  
  # sum over IDs
  mmt <- Reduce("+", lapply(split.data.frame(resids, resids$id), function(obs) {
    # compute weighted residual matrix for individual
    # resid <- obs$weight * obs[, grepl("dtr", names(obs))] * 
    # (obs$Y - meanvec(times, spltime, design, gammas))
    # Sum over DTRs
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      t(matrix(c(dmat[[dtr]], obs$X), nrow = length(times))) %*% 
        solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% t(m)
  })) / n
  
  resids1 = cbind(resids, S[,1:(ncol(S)-1)])
  nvars = (ncol(S)-1)
  
  mst <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      t(matrix(c(dmat[[dtr]], obs$X), nrow = length(times))) %*% 
        solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    m %*% as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])
  })) / n
  
  smt <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    m <- Reduce("+", lapply(1:nDTR, function(dtr) {
      t(matrix(c(dmat[[dtr]], obs$X), nrow = length(times))) %*% 
        solve(V[[dtr]]) %*%
        as.matrix(obs[[paste0("dtr", dtr)]], ncol = 1)
    }))
    t(as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])) %*%
      t(m)
  })) / n
  
  sst <- Reduce("+", lapply(split.data.frame(resids1, resids1$id), function(obs) {
    t(as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])) %*% 
      as.matrix(obs[1,(ncol(resids1) - (nvars-1)):ncol(resids1)])
  })) / n
  
  if(sum(sst) == 0){
    mmt
  } else{
    mmt - mst%*%solve(sst)%*%smt
    
  }
  
}

########################### Generic Functions Needed ###########################
### Note: there are more because like scenario 2, we need to calculate the
### estimated weights within this function. 

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

mod.derivs <- function(times, spltime, design) {
  nDTR <- switch(design, 8, 4, 3)
  A <- dtrIndex(design)
  
  dmat <- lapply(1:nDTR, function(dtr) {
    do.call(rbind, lapply(times, function(ti) {
      matrix(c(
        1,
        (ti <= spltime) * matrix(c(ti, ti * A$a1[dtr], rep(0, 6)), nrow = 1) +
          (ti > spltime) * matrix(
            c(
              spltime,
              spltime * A$a1[dtr],
              (ti - spltime),
              (ti - spltime) * A$a1[dtr],
              (ti - spltime) * A$a2R[dtr],
              (ti - spltime) * A$a2NR[dtr],
              (ti - spltime) * A$a1[dtr] * A$a2R[dtr],
              (ti - spltime) * A$a1[dtr] * A$a2NR[dtr]
            ),
            nrow = 1
          )
      ),
      nrow = 1)
    }))
  })
  
  # Since we constructed dmat using the model for design 1, and the models for
  # designs 2 and 3 are nested inside the design 1 model (see paper supplement),
  # we can remove the columns of dmat that are all zeros, as these correspond to
  # parameters that aren't in the nested model.
  
  # Find columns where all entries are zero for all DTRs
  removeCols <- vector("integer")
  if (design == 2)
    removeCols <- c(6, 8)
  else if (design == 3)
    removeCols <- c(6, 8, 9)
  
  # Remove those columns (if necessary)
  if (length(removeCols) != 0) {
    lapply(dmat, function(x) x[, -removeCols])
  } else
    dmat
}

varmat_adj <- function(sigma2, rho, times, design,
                       corstr = c("identity", "exchangeable", 
                                  "ar1", "unstructured")) {
  
  nDTR <- switch(design, 8, 4, 3)
  corstr <- match.arg(corstr)
  
  sigma2 <- reshapeSigma(sigma2, times, design)
  
  if (length(rho) == 1 & corstr == "unstructured") {
    warning(
      paste(
        "rho is length 1 with unstructured corstr.",
        "Setting corstr to exchangeable."
      ))
    corstr <- "exchangeable"
  } else if (corstr == "unstructured") {
    if (!is.matrix(rho) | (is.matrix(rho) &
                           !all.equal(dim(rho), c(length(times), nDTR)))) {
      stop("rho is not of proper dimension: must be T by nDTR")
    } else {
      lapply(1:nDTR, function(dtr) {
        diag(sqrt(sigma2[, dtr])) %*%
          cormat(rho[, dtr], length(times), corstr) %*% 
          diag(sqrt(sigma2[, dtr]))
      })
    }
  } else {
    if(length(rho) == 1){
      lapply(1:nDTR, function(dtr) {
        diag(sqrt(sigma2[, dtr])) %*% 
          cormat(rho, length(times), corstr) %*% diag(sqrt(sigma2[, dtr]))
      })
    } else{
      lapply(1:nDTR, function(dtr) {
        diag(sqrt(sigma2[, dtr])) %*% 
          cormat(rho[dtr], length(times), corstr) %*% diag(sqrt(sigma2[, dtr]))
      })
    }
  }
}
