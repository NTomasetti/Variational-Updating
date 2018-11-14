rm(list = ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id <- as.numeric(repenv)

library(plyr, lib.loc = 'packages')
library(ks, lib.loc = 'packages')
library(tibble, lib.loc = 'packages')
library(mvtnorm, lib.loc = 'packages')

eucdist <- function(y,z, w) {
  t(y-z) %*% w %*% (y-z)
}

drawGK <-function(theta, T, e0 = 0){
  
  a <- theta[1]
  b <- theta[2]
  g <- theta[3]
  k <- theta[4]
  
  xy <- rnorm(T)
  eps <- b*(1+.8*((1-exp(-g*xy))/(1+exp(-g*xy))))*((1+xy^2)^k)*xy
  y <- a * e0 + eps[1]
  
  if(T > 1){
    for(t in 2:T){
      y[t] <- a * eps[t-1] + eps[t]
    }
  }
  data.frame(y = y, eps = eps)
}

abcVB <- function(lambda, z, draws, post, dim = 4, mix = ncol(lambda)){
  draws <- as.matrix(draws)
 
  
  Jacobian <- apply(draws, 1, function(x){
    (1  / x[1] + 1 / (1 - x[1]))
  })
  draws[, 1] <-  -log((1 - draws[,1]) / draws[,1])
  
  maxIter <- 1000
  threshold <- 0.0001
  alpha <- 0.01
  beta1 <- 0.9
  beta2 <- 0.99
  e <- 1e-8
  
  LB <- rep(0, maxIter)
  meanLB <- 5
  iter <- 1
  diff <- threshold + 1
  MZ <- VZ <- rep(0, mix)
  M <- V <- matrix(0, nrow(lambda), mix)
  
  while(diff > threshold & iter <= maxIter){
    
    mean <- lambda[1:dim, ]
    Sigma <- array(0, dim = c(dim, dim, mix))
    SigInv <- Sigma
    
    for(m in 1:mix){
      U <- matrix(lambda[dim + 1:dim^2, m], dim)
      Sigma[,,m] <- t(U) %*% U
      if(det(Sigma[,,m]) > 1e-8){
        SigInv[,,m] <- solve(Sigma[,,m])
      } else {
        SigInv[,,m] <- diag(1/diag(Sigma[,,m]))
      }
    }
   
  
    pi <- exp(z) / sum(exp(z))
    qComp <- matrix(0, nrow(draws), mix)
    for(m in 1:mix){
      qComp[,m] <- mvtnorm::dmvnorm(draws, mean[,m], Sigma[,,m])
    }
    qDens <- c(qComp %*% pi)
    qDensTransf <- qDens * Jacobian
    
    scorePi <- qComp / qDens
    
    scoreZ <- matrix(0, nrow(draws), mix)
    denom <- sum(exp(z))^2
  
    for(i in 1:mix){
      for(j in 1:mix){
        if(i == j){
          scoreZ[, i] <- scoreZ[, i] + scorePi[, j] * sum(exp(z[i] + z[-i])) / denom
        } else {
          scoreZ[, i] <- scoreZ[, i] - scorePi[, j] * exp(z[i] + z[j]) / denom
        }
      }
    }
    
    score <- array(0, dim = c(nrow(lambda), mix, nrow(draws)))
    for(m in 1:mix){
      for(n in 1:nrow(draws)){
        
        meandiff <- (draws[n, ] - mean[,m])
        score[1:dim, m, n] <- SigInv[,,m] %*% meandiff
        product <- SigInv[,,m] %*% meandiff %*% t(meandiff) %*% SigInv[,,m] 
        scoreSig <- -SigInv[,,m] + diag(diag(SigInv[,,m]) )/ 2 + product - diag(diag(product))/2
        
        score[5, m, n] <- lambda[c(5, 9, 13, 17), m] %*% scoreSig[1, 1:4] + lambda[5, m] * scoreSig[1, 1] # U_11
        score[9, m, n] <- lambda[c(5, 9, 13, 17), m] %*% scoreSig[2, 1:4] + lambda[9, m] * scoreSig[2, 2] # U_21
        score[10, m, n] <- lambda[c(10, 14, 18), m] %*% scoreSig[2, 2:4] + lambda[10, m] * scoreSig[2, 2] # U_22
        score[13, m, n] <- lambda[c(5, 9, 13, 17), m] %*% scoreSig[3, 1:4] + lambda[13, m] * scoreSig[3, 3] # U_31
        score[14, m, n] <- lambda[c(10, 14, 18), m] %*% scoreSig[3, 2:4] + lambda[14, m] * scoreSig[3, 3] # U_32
        score[15, m, n] <- lambda[c(15, 19), m] %*% scoreSig[3, 3:4] + lambda[15, m] * scoreSig[3, 3] # U_33
        score[17, m, n] <- lambda[c(5, 9, 13, 17), m] %*% scoreSig[4, 1:4] + lambda[17, m] * scoreSig[4, 4] # U_41
        score[18, m, n] <- lambda[c(10, 14, 18), m] %*% scoreSig[4, 2:4] + lambda[18, m] * scoreSig[4, 4] # U_42
        score[19, m, n] <- lambda[c(15, 19), m] %*% scoreSig[4, 3:4] + lambda[19, m] * scoreSig[4, 4] # U_43
        score[20, m, n] <- 2 * lambda[20, m] * scoreSig[4, 4] # U_44
        
        score[, m, n] <- score[, m, n] * pi[m] * scorePi[n, m] # Convert from dlog N_i / dlam to dlog N_sum /dlam
      }
    }


    w <- qDensTransf / post
    LB[iter] <- mean(w * (log(post) - log(qDensTransf)))
    
    gradient <- w * score * (log(post) - log(qDensTransf))
    gradient <- apply(gradient, 1:2, mean)
    gradientSq <- gradient^2
    
    gradZ <- w * scoreZ * (log(post) - log(qDensTransf))
    gradZ <- apply(gradZ, 2, mean)
    gradZSq <- gradZ^2
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mstar <- M / (1 - beta1^iter)
    Vstar <- V / (1 - beta2^iter)
    update <- alpha * Mstar / (sqrt(Vstar) + e)
    if(any(is.na(update))){
      print('Break Lambda')
      break
    }
    lambda <- lambda + update
  
    
    MZ <- beta1 * MZ + (1 - beta1) * gradZ
    VZ <- beta2 * VZ + (1 - beta2) * gradZSq
    Mstar <- MZ / (1 - beta1^iter)
    Vstar <- VZ / (1 - beta2^iter)
    update <- alpha * Mstar / (sqrt(Vstar) + e)
    if(any(is.na(update))){
      print('Break Z')
      break
    }
    z <- z + update
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter- 4)])
      diff <- abs(meanLB - oldMeanLB)
    }
    #if(iter %% 100 == 0){
    #  print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    #}
    iter <- iter + 1
  }
 # print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
  list(lambda = lambda, z = z, iter = iter-1, LB = LB[1:(iter-1)])
}

set.seed(id)
a <- 0.5 #AR parameter
b <- 1 #scale
g <- 2 #skewness
k <- 0.5 # kurtosis

T <- 500
Tseq <- c(0, seq(50, T, 10))
reps <- 50000
keep <- 0.01

data <- drawGK(c(a, b, g, k), T)
y <- data$y
eps <- data$eps

s0 <- function(x, qr) {T <- length(x); 1/(T-1) * sum(x[1:(T-1)] * x[2:T])}
s1 <- function(x, qr) {T <- length(x); (1/T*sum(x))}
s3 <- function(x, qr) {T <- length(x); 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)}
s4 <- function(x, qr) {T <- length(x); (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5}
s1r <- function(x, qr) qr[4]
s2r <- function(x, qr) qr[6]-qr[2]
s3r <- function(x, qr) (qr[6]+qr[2]-2*qr[4])/(qr[6]-qr[2])
s4r <- function(x, qr) (qr[7]-qr[5]+qr[3]-qr[1])/(qr[6]-qr[2])
sumStats <- plyr::each(s0 = s0, s1 = s1r, s2 = s2r, s3 = s3, s4 = s4, s5 = s3r, s6 = s4r)

results <- tibble()
diverge <- FALSE

for(t in 2:8){#(length(Tseq)-1)){
  print(paste(t, Sys.time()))
  ySub <- y[(Tseq[t-1]+1):Tseq[t]]
  yFull <- y[1:Tseq[t]]
  
  # Full ABC Fit
  qrY <- quantile(yFull, c(.08,.25,.36,.5,.6,.75,.875))
  yss <- sumStats(yFull, qrY)
  
  reps <- 50000
  keep <- 0.01
  theta <- cbind(runif(reps, 0, 1), runif(reps, 0, 10), runif(reps, 0, 10), runif(reps, 0, 10))
  zss <- matrix(0, reps, 7)
  #loop through simulation of data
  for(iter in 1:reps){
    z <- drawGK(theta[iter, ], Tseq[t])$y
    qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
    zss[iter, ] <- sumStats(z, qrZ)
  }
  varZ <- apply(zss, 2, var)
  distance <- rep(0, reps)
  for(iter in 1:reps){
    distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
  }
  
  ABC <- as.tibble(cbind(theta, distance))
  names(ABC) = c('a', 'b', 'g', 'k', 'dist')
  
  accept <- ABC[ABC$dist < quantile(ABC$dist, keep), 1:4]
  
  if(t == 2){
    acceptVB <- accept
    acceptSABC <- accept
    acceptSABC2 <- accept
  } else {
    
    # Sequential ABC
    qrY <- quantile(ySub, c(.08,.25,.36,.5,.6,.75,.875))
    yss <- sumStats(ySub, qrY)
    
    # Draw Theta
    reps <- 5000
    keep <- 0.1
    subset <- sample(1:nrow(acceptSABC), reps, replace = T)
    theta <- as.matrix(acceptSABC)[subset, ]
    zss <- matrix(0, reps, 7)
    #loop through simulation of data
    for(iter in 1:reps){
      z <- drawGK(theta[iter, ], Tseq[t] - Tseq[t-1], eps[Tseq[t-1]])$y
      qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
      zss[iter, ] <- sumStats(z, qrZ)
    }
    varZ <- apply(zss, 2, var)
    distance <- rep(0, reps)
    for(iter in 1:reps){
      distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
    }
    SABC <- as.tibble(cbind(theta, distance))
    names(SABC) = c('a', 'b', 'g', 'k', 'dist')
    
    acceptSABC <- SABC[SABC$dist < quantile(SABC$dist, keep), 1:4]
    
    # Sequential ABC with intermittent variatonal updates
    distinct <- length(unique(acceptSABC2$a))
    
    if(distinct >= nrow(acceptSABC) / 10){
      if(!diverge){
        acceptSABC2 <- acceptSABC
      } else {
        # Normal S-ABC Iteration
        reps <- 5000
        keep <- 0.1
        subset <- sample(1:nrow(acceptSABC2), reps, replace = T)
        theta <- as.matrix(acceptSABC2)[subset, ]
        zss <- matrix(0, reps, 7)
        #loop through simulation of data
        for(iter in 1:reps){
          z <- drawGK(theta[iter, ], Tseq[t] - Tseq[t-1], eps[Tseq[t-1]])$y
          qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
          zss[iter, ] <- sumStats(z, qrZ)
        }
        varZ <- apply(zss, 2, var)
        distance <- rep(0, reps)
        for(iter in 1:reps){
          distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
        }
        SABC2 <- as.tibble(cbind(theta, distance))
        names(SABC2) = c('a', 'b', 'g', 'k', 'dist')
        
        acceptSABC2 <- SABC2[SABC2$dist < quantile(SABC2$dist, keep), 1:4]
        
      }
     
    } else {
      diverge <- TRUE
      # VS-ABC Iteration
      # Fit KDE
      grid <- do.call(expand.grid, lapply(acceptSABC2, quantile, prob=seq(0.01, 0.99, length.out = 14)))
      
      postDensity <- ks::kde(as.matrix(acceptSABC2), eval.points = grid)
      grid$post <- c(postDensity$estimate)
      
      # Subset from KDE for VB algo
      normProb <- grid$post
      normProb[normProb <= 0] <- 0
      normProb <- normProb / sum(normProb)
      subset <- sample(1:nrow(grid), 50, prob = normProb)
      draws <- grid[subset, 1:4] 
      postDens <- grid[subset, ]$post 
      
      # Fit VB
      mix <- 5
      z <- rep(0, mix)
      lambda <- matrix(0, 4 * 5, mix)
      lambda[,1] <- c(colMeans(acceptSABC2), chol(cov(acceptSABC2)))
      for(m in 2:mix){
        lambda[,m] <- c(colMeans(acceptSABC2) +  c(0.5, 0.3, 0.5, 0.3) * (-1)^m + rnorm(4, 0, 0.1), chol(cov(acceptSABC2)))
      }
      vb2 <- abcVB(lambda, z, draws, postDens)
      
      # Set up ABC
      qrY <- quantile(ySub, c(.08,.25,.36,.5,.6,.75,.875))
      yss <- sumStats(ySub, qrY)
      
      # Draw Theta
      weights <- vb2$z
      pi <- exp(weights) / sum(exp(weights))
      mean <- vb2$lambda[1:4, ]
      Sig <- array(0, dim = c(4, 4, mix))
      for(m in 1:mix){
        U <- matrix(vb2$lambda[5:20, m], 4)
        Sig[,,m]<- t(U) %*% U
      }
      reps <- 5000
      keep <- 0.1
      theta <- matrix(0, reps, 4)
      
      for(iter in 1:reps){
        okay <- FALSE
        # Truncate distribution to [0, 1], [0, 10]^3
        while(!okay){
          u <- runif(1)
          component <- min(which(cumsum(pi) > u))
          draw <- mvtnorm::rmvnorm(1, mean[,component], Sig[,,component])
          if(all(draw[2:4] > 0) & all(draw[2:4] < 10) & draw[1] > 0 & draw[1] < 1){
            okay <- TRUE
          }
        }
        theta[iter, ] <- draw
      }
      zss <- matrix(0, reps, 7)
      #loop through simulation of data
      for(iter in 1:reps){
        z <- drawGK(theta[iter, ], Tseq[t] - Tseq[t-1], eps[Tseq[t-1]])$y
        qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
        zss[iter, ] <- sumStats(z, qrZ)
      }
      varZ <- apply(zss, 2, var)
      distance <- rep(0, reps)
      for(iter in 1:reps){
        distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
      }
      
      SABC2 <- as.tibble(cbind(theta, distance))
      names(SABC2) = c('a', 'b', 'g', 'k', 'dist')
      
      acceptSABC2 <- SABC2[SABC2$dist < quantile(SABC2$dist, keep), 1:4]
    }
    
    # Fit KDE
    grid <- do.call(expand.grid, lapply(acceptVB, quantile, prob=seq(0.01, 0.99, length.out = 14)))
  
    postDensity <- ks::kde(as.matrix(acceptVB), eval.points = grid)
    grid$post <- c(postDensity$estimate)
    
    # Subset from KDE for VB algo
    normProb <- grid$post
    normProb[normProb <= 0] <- 0
    normProb <- normProb / sum(normProb)
    subset <- sample(1:nrow(grid), 100, prob = normProb)
    draws <- grid[subset, 1:4] 
    postDens <- grid[subset, ]$post 
    
    # Fit VB
    mix <- 5

    drawsTransf <- matrix(0, nrow(acceptVB), 4)
    acceptVB <- as.matrix(acceptVB)
    drawsTransf[,2:4] <- acceptVB[,2:4]
    drawsTransf[,1] <- -log((1 - acceptVB[,1]) / acceptVB[,1])

    z <- rep(0, mix)
    lambda <- matrix(0, 4 * 5, mix)
    lambda[,1] <- c(colMeans(drawsTransf), chol(cov(drawsTransf)))
    for(m in 2:mix){
      lambda[,m] <- c(colMeans(drawsTransf) + 1 * (-1)^m + rnorm(4, 0, 0.1), chol(cov(drawsTransf)))
    }
    vb <- abcVB(lambda, z, draws, postDens)
    
    # Set up ABC
    qrY <- quantile(ySub, c(.08,.25,.36,.5,.6,.75,.875))
    yss <- sumStats(ySub, qrY)
    
    # Draw Theta
    weights <- vb$z
    pi <- exp(weights) / sum(exp(weights))
    mean <- vb$lambda[1:4, ]
    Sig <- array(0, dim = c(4, 4, mix))
    for(m in 1:mix){
      U <- matrix(vb$lambda[5:20, m], 4)
      Sig[,,m]<- t(U) %*% U
    }
    reps <- 5000
    keep <- 0.1
    theta <- matrix(0, reps, 4)
    
    for(iter in 1:reps){
      okay <- FALSE
      u <- runif(1)
      component <- min(which(cumsum(pi) > u))
    
      while(!okay){
        theta[iter, ] <- mvtnorm::rmvnorm(1, mean[,component], Sig[,,component])
        if(all(theta[iter, 2:4] > 0) & all(theta[iter, 2:4] < 0)){
          okay <- TRUE
        }
      }
       
    }
    theta[,1] <- 1 / (1 + exp(-theta[,1]))
    
    zss <- matrix(0, reps, 7)
    #loop through simulation of data
    for(iter in 1:reps){
      z <- drawGK(theta[iter, ], Tseq[t] - Tseq[t-1], eps[Tseq[t-1]])$y
      qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
      zss[iter, ] <- sumStats(z, qrZ)
    }
    varZ <- apply(zss, 2, var)
    distance <- rep(0, reps)
    for(iter in 1:reps){
      distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
    }
    
    
    ABCVB <- as.tibble(cbind(theta, distance))
    names(ABCVB) = c('a', 'b', 'g', 'k', 'dist')
    
    acceptVB <- ABCVB[ABCVB$dist < quantile(ABCVB$dist, keep), 1:4]
  }
  
  # Forecast ABC results
  fcSims <- matrix(0, nrow(accept), 10)
  for(i in 1:nrow(accept)){
    fcSims[i, ] <- drawGK(unlist(accept[i, ]), 10, eps[Tseq[t]])$y
  }
  densFuns <- apply(fcSims, 2, function(x){
    xSub <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
    dens <- density(xSub)
    approxfun(dens)
  })

  ls <- vapply(1:10, function(x) log(densFuns[[x]](y[Tseq[t]+x])), runif(1))
  results <- rbind(results,
                   tibble(ls = ls,
                          t = Tseq[t],
                          h = 1:10,
                          method = 'ABC',
                          unique = 500))
  if(t > 2){
    fcSims <- matrix(0, nrow(acceptSABC), 10)
    for(i in 1:nrow(acceptSABC)){
      fcSims[i, ] <- drawGK(unlist(acceptSABC[i, ]), 10, eps[Tseq[t]])$y
    }
    densFuns <- apply(fcSims, 2, function(x){
      xSub <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
      dens <- density(xSub)
      approxfun(dens)
    })
    
    ls <- vapply(1:10, function(x) log(densFuns[[x]](y[Tseq[t]+x])), runif(1))
    results <- rbind(results,
                     tibble(ls = ls,
                            t = Tseq[t],
                            h = 1:10,
                            method = 'S-ABC',
                            unique = length(unique(acceptSABC$a))))
    
    fcSims <- matrix(0, nrow(acceptSABC2), 10)
    for(i in 1:nrow(acceptSABC2)){
      fcSims[i, ] <- drawGK(unlist(acceptSABC2[i, ]), 10, eps[Tseq[t]])$y
    }
    densFuns <- apply(fcSims, 2, function(x){
      xSub <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
      dens <- density(xSub)
      approxfun(dens)
    })
    
    ls <- vapply(1:10, function(x) log(densFuns[[x]](y[Tseq[t]+x])), runif(1))
    results <- rbind(results,
                     tibble(ls = ls,
                            t = Tseq[t],
                            h = 1:10,
                            method = 'VS-ABC-2',
                            unique = length(unique(acceptSABC2$a))))
    
    
    
    fcSims <- matrix(0, nrow(acceptVB), 10)
    for(i in 1:nrow(acceptVB)){
      fcSims[i, ] <- drawGK(unlist(acceptVB[i, ]), 10, eps[Tseq[t]])$y
    }
    densFuns <- apply(fcSims, 2, function(x){
      xSub <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
      dens <- density(xSub)
      approxfun(dens)
    })
    
    ls <- vapply(1:10, function(x) log(densFuns[[x]](y[Tseq[t]+x])), runif(1))
    results <- rbind(results,
                     tibble(ls = ls,
                            t = Tseq[t],
                            h = 1:10,
                            method = 'VS-ABC',
                            unique = 500))
  }
}
write.csv(results, paste0('VSABC/rep', id, '.csv'), row.names = FALSE)

library(tidyverse)
results <- tibble()
for(i in 1:200){
  file <- paste0('ma1gk/rep', i, '.csv')
  if(file.exists(file)){
    temp <- read_csv(file, col_types = cols())
    temp$id <- i
    results <- rbind(results, temp)
  }
}
length(unique(results$id))

results %>%
  group_by(t, h, method) %>%
  summarise(meanls = mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, meanls, colour = method)) + 
  facet_wrap(~h, scales = 'free', ncol = 5) +
  theme_bw() + theme(legend.position = 'bottom') + labs(colour = 'Method', x = 'T', y = 'Mean Forecast Logscore')

results %>%
  filter(method %in% c('S-ABC', 'VS-ABC-2') & h == 1) %>%
  group_by(t, method) %>%
  summarise(mid = mean(unique),
            lower = min(unique),
            upper = max(unique)) %>%
  ggplot() + geom_line(aes(t, mid, colour = method)) + geom_ribbon(aes(t, ymin =  lower, ymax = upper, colour = method), alpha = 0.2)

results %>%
  filter(method %in% c('S-ABC', 'VS-ABC-2') & h == 1) %>%
  select(-ls) %>%
  spread(method, unique) %>%
  mutate(ratio = `VS-ABC-2` / `S-ABC`) %>%
  select(t, ratio, id) %>%
  inner_join(
    results %>%
      filter(method %in% c('S-ABC', 'VS-ABC-2') & h == 1) %>%
      select(-unique) %>%
      spread(method, ls) %>%
      mutate(diff = `VS-ABC-2` - `S-ABC`) %>%
      select(t, diff, id)
  ) %>%
  group_by(t) %>%
  summarise(diff = mean(diff, na.rm = TRUE),
            ratio = mean(ratio)) %>%
  gather(var, value, -t) %>%
  ggplot() + geom_line(aes(t, value)) + facet_wrap(~var, scales = 'free')
  
  
accept %>% 
  as.data.frame() %>%
  mutate(method = 'ABC') %>% 
  rbind(acceptVB %>% 
          as.data.frame() %>%
          mutate(method = 'VB'),
        acceptSABC %>% 
          as.data.frame() %>%
          mutate(method = 'SABC'),
        acceptSABC2 %>% 
          as.data.frame() %>%
          mutate(method = 'SABC2')) %>%
  gather(var, value, -method) %>%
  ggplot() + geom_density(aes(value, colour =method)) + facet_wrap(~var, scales = 'free')

# Forecast ABC results
fcSims <- matrix(0, nrow(accept), 10)
for(i in 1:nrow(accept)){
  fcSims[i, ] <- drawGK(unlist(accept[i, ]), 10, eps[Tseq[t]])$y
}


fcSimsSABC <- matrix(0, nrow(acceptSABC), 10)
for(i in 1:nrow(acceptSABC)){
  fcSimsSABC[i, ] <- drawGK(unlist(acceptSABC[i, ]), 10, eps[Tseq[t]])$y
}
fcSimsSABC2 <- matrix(0, nrow(acceptSABC2), 10)
for(i in 1:nrow(acceptSABC2)){
  fcSimsSABC2[i, ] <- drawGK(unlist(acceptSABC2[i, ]), 10, eps[Tseq[t]])$y
}

  
fcSimsVB <- matrix(0, nrow(acceptVB), 10)
for(i in 1:nrow(acceptVB)){
  fcSimsVB[i, ] <- drawGK(unlist(acceptVB[i, ]), 10, eps[Tseq[t]])$y
}
acceptT <- data.frame(a = rep(a, 500), b = rep(b, 500), g = rep(g, 500), k = rep(k, 500))

fcSimsT <- matrix(0, nrow(acceptT), 10)
for(i in 1:nrow(accept)){
  fcSimsT[i, ] <- drawGK(unlist(acceptT[i, ]), 10, eps[Tseq[t]])$y
}

acceptVB2 <- acceptVB
acceptVB2$b <- b
acceptVB2$a <- a
fcSimsVB2 <- matrix(0, nrow(acceptVB2), 10)
for(i in 1:nrow(acceptVB)){
  fcSimsVB2[i, ] <- drawGK(unlist(acceptVB2[i, ]), 10, eps[Tseq[t]])$y
}



fc <- data.frame(ABC = fcSims[,10], SABC = fcSimsSABC[,10], SABC2 = fcSimsSABC2[,10], VB = fcSimsVB[,10], True = fcSimsT[,10],
                 VB2 = fcSimsVB2[,10])
fc %>% 
  gather(var, value) %>%
  ggplot() + geom_density(aes(value, colour = var)) +
  geom_vline(aes(xintercept = y[190])) + xlim(-10, 30)

