# This script contains auxilliary functions. 

dtrNames <- function(design) {
  DTRs <- dtrIndex(design)
  
  dtrTriples <- as.data.frame((do.call(cbind, DTRs) + 1) / 2)
  dtrTriples[dtrTriples == 0.5] <- 0
  dtrTriples <- apply(dtrTriples, 1, function(x) paste(x, collapse = ""))
  dtrTriples
}

dtrIndex <- function(design) {
  if (design == 1) {
    a1   <- rep(c(1, -1), each = 4)
    a2R  <- rep(rep(c(1, -1), each = 2), 2)
    a2NR <- rep(c(1, -1), 4)
  } else if (design == 2) {
    a1   <- rep(c(1, -1), each = 2)
    a2R  <- rep(0, 4)
    a2NR <- rep(c(1, -1), 2)
  } else if (design == 3) {
    a1   <- c(1, 1, -1)
    a2R  <- c(0, 0, 0)
    a2NR <- c(1, -1, 0)
  }
  return(list("a1" = a1, "a2R" = a2R, "a2NR" = a2NR))
}

response.beta <- function(d, gammas, r1, r0, respDirection = c("high", "low"),
                          sigma, causal = F) {
  shape1 <- r1 / (1 - r1)
  shape0 <- r0 / (1 - r0)
  
  respDirection <- match.arg(respDirection)
  tail <- switch(respDirection, "high" = F, "low" = T)
  
  if (causal) {
    x1 <- pnorm(d$Y1.1, mean = gammas[1] + gammas[2] + gammas[3],
                sd = sigma, lower.tail = tail)
    x0 <- pnorm(d$Y1.0, mean = gammas[1] + gammas[2] - gammas[3],
                sd = sigma, lower.tail = tail)
    
    d$respProb.1 <- qbeta(x1, shape1 = shape1, shape2 = 1)
    d$respProb.0 <- qbeta(x0, shape1 = shape0, shape2 = 1)
    
    d$R.1 <- sapply(1:nrow(d), function(i) rbinom(1, 1, d$respProb.1[i]))
    d$R.0 <- sapply(1:nrow(d), function(i) rbinom(1, 1, d$respProb.0[i]))
  } else {
    x <- pnorm(d$Y1, mean = gammas[1] + gammas[2] + gammas[3]*d$A1, 
               sd = sigma, lower.tail = tail)
    respProb <- vector("numeric", nrow(d))
    
    respProb[d$A1 ==  1] <- qbeta(x[d$A1 ==  1], shape1 = shape1, shape2 = 1)
    respProb[d$A1 == -1] <- qbeta(x[d$A1 == -1], shape1 = shape0, shape2 = 1)
    d$respProb <- respProb
    d$R <- sapply(1:nrow(d), function(i) rbinom(1, 1, respProb[i]))
  }
  
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.indep <- function(d, gammas, r1, r0, respDirection = NULL,
                           causal = F, ...) {
  if (causal) {
    d$R.1 <- rbinom(nrow(d), 1, r1)
    d$R.0 <- rbinom(nrow(d), 1, r0)
  } else{
    d$R <- NA
    d$R[d$A1 ==  1] <- rbinom(sum(d$A1 ==  1), 1, r1)
    d$R[d$A1 == -1] <- rbinom(sum(d$A1 == -1), 1, r0)
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.oneT <- function(d, gammas, r1, r0, respDirection = c("high", "low"),
                          sigmaY1, causal = F) {
  respDirection <- match.arg(respDirection)
  tail <- switch(respDirection, "high" = F, "low" = T)
  
  if (length(sigmaY1) == 1) {
    sigmaY1.1 <- sigmaY1.0 <- sigmaY1
  } else if (length(sigmaY1) == 2) {
    sigmaY1.0 <- sigmaY1[1]
    sigmaY1.1 <- sigmaY1[2]
  } else {
    warning("sigmaY1 can have length at most 2; ignoring subsequent elements.")
  }
  
  upsilon <- qnorm(r1, as.numeric(sum(gammas[1:3])), 
                   sigmaY1.1, lower.tail = tail)
  if (causal) {
    d$R.0 <- as.numeric(d$Y1.0 >= upsilon)
    d$R.1 <- as.numeric(d$Y1.1 >= upsilon)
  } else {
    d$R <- as.numeric(d$Y1 >= upsilon)
  }
  r0temp <- pnorm(upsilon, sum(gammas[1:2]) - gammas[3],
                  sigmaY1.0, lower.tail = tail)
  if (r0temp != r0) {
    warning(paste("Overwriting the provided value of r0 to accomodate the",
                  "single-threshold response function.",
                  "The provided value is", round(r0, 3), "and the new value is",
                  round(r0temp, 3)))
    r0 <- r0temp
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.sq <- function(d, gammas, r1, r0, respDirection = c("high", "low"), 
                        sigma, causal = F) {
  respDirection <- match.arg(respDirection)
  
  upsilon1.UB <- qnorm(r1/2, sum(gammas[1:3]),sigma,
                       lower.tail = F)
  upsilon0.UB <- qnorm(r0/2, sum(gammas[1:2]) - gammas[3], sigma,
                       lower.tail = F)
  upsilon1.LB <- qnorm(r1/2, sum(gammas[1:3]), sigma)
  upsilon0.LB <- qnorm(r0/2, sum(gammas[1:2]) - gammas[3], sigma)
  
  if (causal) {
    d$R.1 <- d$R.0 <- NA
    if (respDirection == "high") {
      d$R.1 <- as.numeric(d$Y1.1 >= upsilon1.UB | d$Y1.1 <= upsilon1.LB)
      d$R.0 <- as.numeric(d$Y1.0 >= upsilon0.UB | d$Y1.0 <= upsilon0.LB)
    } else {
      d$R.1 <- as.numeric(d$Y1.1 >= -upsilon1.LB & d$Y1.1 <= upsilon1.UB)
      d$R.0 <- as.numeric(d$Y1.0 >= -upsilon0.LB & d$Y1.0 <= upsilon0.UB)
    }
  } else {
    d$R <- NA
    if (respDirection == "high") {
      d$R[d$A1 ==  1] <-
        as.numeric(d$Y1[d$A1 ==  1] >= upsilon1.UB |
                     d$Y1[d$A1 ==  1] <= upsilon1.LB)
      d$R[d$A1 == -1] <-
        as.numeric(d$Y1[d$A1 == -1] >= upsilon0.UB |
                     d$Y1[d$A1 == -1] <= upsilon0.LB)
    } else {
      d$R[d$A1 ==  1] <-
        as.numeric(d$Y1[d$A1 ==  1] >= -upsilon1.LB &
                     d$Y1[d$A1 ==  1] <= upsilon1.UB)
      d$R[d$A1 == -1] <-
        as.numeric(d$Y1[d$A1 == -1] >= -upsilon0.LB &
                     d$Y1[d$A1 == -1] <= upsilon0.UB)
    }
  }
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

response.twoT <- function(d, gammas, r1, r0, respDirection = c("high", "low"),
                          sigmaY1, causal = F) {
  respDirection <- match.arg(respDirection)
  tail <- switch(respDirection, "high" = F, "low" = T)
  
  if (length(sigmaY1) == 1) {
    sigmaY1.1 <- sigmaY1.0 <- sigmaY1
  } else if (length(sigmaY1) == 2) {
    sigmaY1.0 <- sigmaY1[1]
    sigmaY1.1 <- sigmaY1[2]
  } else {
    warning("sigmaY1 can have length at most 2; ignoring subsequent elements.")
  }
  
  upsilon1 <- qnorm(r1, sum(gammas[1:3]), sigmaY1.1,
                    lower.tail = tail)
  upsilon0 <- qnorm(r0, sum(gammas[1:2]) - gammas[3], sigmaY1.0,
                    lower.tail = tail)
  
  if (causal) {
    d$R.1 <- d$R.0 <- NA
    d$R.1 <- as.numeric(d$Y1.1 >= upsilon1)
    d$R.0 <- as.numeric(d$Y1.0 >= upsilon0)
  } else {
    d$R <- NA
    d$R[d$A1 ==  1] <- as.numeric(d$Y1[d$A1 ==  1] >= upsilon1)
    d$R[d$A1 == -1] <- as.numeric(d$Y1[d$A1 == -1] >= upsilon0)
  }
  
  return(list("data" = d, "r1" = r1, "r0" = r0))
}

reshapeSigma <- function(sigma, times, design) {
  nDTR <- switch(design, 8, 4, 3)
  
  sigma <- as.matrix(sigma)
  
  # Restructure sigma to a length(times)-by-nDTR matrix
  if (sum(dim(sigma) == c(1, 1)) == 2) {
    # scalar case
    sigma <- matrix(rep(sigma, length(times) * nDTR), nrow = length(times))
  } else if (nrow(sigma) == length(times) & 
             length(grep("time", rownames(sigma))) == length(times) & 
             ncol(sigma) == 1) {
    # case in which there is one row per time
    # (i.e., pool.time = F, pool.dtr = T)
    sigma <- matrix(rep(sigma, nDTR), ncol = nDTR, byrow = F)
  } else if (nrow(sigma) == nDTR & 
             length(grep("dtr", rownames(sigma))) == nDTR & 
             ncol(sigma) == 1) {
    # case in which there is one row per DTR
    # (i.e., pool.time = T, pool.dtr = F)
    sigma <- t(matrix(rep(sigma, length(times)),
                      ncol = length(times), byrow = F))
  }
  
  sigma
}

createDTRIndicators <- function(d, design) {
  if (sum(c("A1", "A2R", "A2NR") %in% names(d)) != 3)
    stop("Must provide data from a SMART.")
  if (design == 1) {
    d$dtr1 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR ==  1)))
    d$dtr2 <- as.numeric(with(d, A1 ==  1 & (A2R ==  1 | A2NR == -1)))
    d$dtr3 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR ==  1)))
    d$dtr4 <- as.numeric(with(d, A1 ==  1 & (A2R == -1 | A2NR == -1)))
    d$dtr5 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR ==  1)))
    d$dtr6 <- as.numeric(with(d, A1 == -1 & (A2R ==  1 | A2NR == -1)))
    d$dtr7 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR ==  1)))
    d$dtr8 <- as.numeric(with(d, A1 == -1 & (A2R == -1 | A2NR == -1)))
  } else if (design == 2) {
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr4 <- as.numeric(with(d, (A1 == -1) * (R + (1 - R) * (A2NR == -1))))
  } else if (design == 3) {
    d$dtr1 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR ==  1))))
    d$dtr2 <- as.numeric(with(d, (A1 ==  1) * (R + (1 - R) * (A2NR == -1))))
    d$dtr3 <- as.numeric(with(d, (A1 == -1)))
  } else stop("design must be one of 1, 2, 3.")
  
  return(d)
}

validTrial <- function(d, design) {
  l <- split.SMART(d)
  if (design == 1) {
    length(l) == 8
  } else if (design == 2) {
    length(l) == 6
  } else if (design == 3) {
    length(l) == 5
  } else stop("design must be in 1-3.")
}

split.SMART <- function(d, marginal = FALSE) {
  if (marginal) {
    l <- lapply(grep("dtr", names(d), value = TRUE), function(s) {
      subset(d, get(s) == 1)
    })
    names(l) <- grep("dtr", names(d), value = TRUE)
  } else {
    l <- split.data.frame(d, list(d$A1, d$R, d$A2R, d$A2NR))
    x <- apply(unique(subset(d, select = c("A1", "R", "A2R", "A2NR"))), 1,
               function(v) {
                 paste(v, collapse = ".")
               })
    l <- lapply(1:length(x), function(i) l[[x[i]]])
    names(l) <- x
  }
  l
}

### Construct a t-by-t correlation matrix
cormat <-
  function(rho,
           t,
           corstr = c("identity", "exchangeable", "ar1", "unstructured")) {
    corstr <- match.arg(corstr)
    if (corstr == "identity") {
      diag(rep(1, t))
    } else if (corstr == "exchangeable") {
      m <- matrix(rep(rho, t ^ 2), nrow = t)
      diag(m) <- rep(1, t)
      m
    } else if (corstr == "ar1") {
      m <- diag(t)
      rho ^ (abs(row(m) - col(m)))
    } else if (corstr == "unstructured") {
      if (length(rho) != t)
        stop("for unstructured corstr, must have rho of length t")
      m <- diag(3)
      tpairs <- combn(1:t, 2)
      for (j in 1:ncol(tpairs)) {
        m[tpairs[1, j], tpairs[2, j]] <-
          m[tpairs[2, j], tpairs[1, j]] <- rho[j]
      }
      m
    }
  }


