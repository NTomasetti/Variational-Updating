rm(list = ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id <- as.numeric(repenv)

library(tibble, lib.loc = 'packages')
library(VineCopula, lib.loc = 'packages')
library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
sourceCpp('particleFilter.cpp')

drawSV <- function(theta, T, v0 = NULL){
  
  t1 <- theta[1]
  t2 <- theta[2]
  
  if(is.null(v0)){
    v0 <- rnorm(1, 0, sqrt(t2 / (1 - t1^2)))
  }
  v <- t1 * v0 +  rnorm(1, 0, sqrt(t2))
  y <- rnorm(1, 0, sqrt(exp(v)))
  
  if(T > 1){
    for(t in 2:T){
      v[t] <- t1 * v[t-1] + rnorm(1, 0, sqrt(t2))
      y[t] <- sqrt(exp(v[t])) * rnorm(1, 0, 1)
    }
  }
  
  data.frame(y = y, v = v)
}

forecastSV <- function(y, theta, H = 10, S = 50){
  v <- rnorm(S, 0, sqrt(theta[2] / (1 - theta[1]^2)))
  w <- dnorm(y[1], 0, sqrt(exp(v)), log = TRUE)
  w <- w - max(w)
  w <- exp(w) / sum(exp(w))
  v <- sample(v, S, prob = w, replace = TRUE)
  
  for(t in 2:length(y)){
    v <- theta[1] * v + rnorm(S, 0, sqrt(theta[2]))
    w <- dnorm(y[t], 0, sqrt(exp(v)), log = TRUE)
    w <- w - max(w)
    w <- exp(w) / sum(exp(w))
    v <- sample(v, S, prob = w, replace = TRUE)
  }
  vFC <- vFC <- theta[1] * v[1] + rnorm(1, 0, sqrt(theta[2]))
  for(t in 2:H){
    vFC <- c(vFC, theta[1] * vFC[t-1] + rnorm(1, 0, sqrt(theta[2])))
  }
  y <- rnorm(H, 0, sqrt(exp(vFC)))
  y
}

mcmcVB <- function(lambda, draws, post, maxIter = 1000, threshold = 0.001, alpha = 0.01){
  
  draws <- as.matrix(draws)
  beta1 <- 0.9
  beta2 <- 0.99
  e <- 1e-8
  
  LB <- rep(0, maxIter)
  meanLB <- 5
  iter <- 1
  diff <- threshold + 1
  M <- V <- rep(0, length(lambda))
  
  while(diff > threshold & iter <= maxIter){

    mean <- lambda[1:3]
    U <- matrix(lambda[4:12], 3)
    Sigma <- t(U) %*% U
    SigInv <- solve(Sigma)
    
    #qDensTransf <- mvtnorm::dmvnorm(drawsTransf, mean, Sigma)
    #qDens <- qDensTransf * Jacobian
    qDens <- mvtnorm::dmvnorm(draws, mean, Sigma)
    
    score <- matrix(0, 12, nrow(draws))
    for(n in 1:nrow(draws)){
      #meandiff <- drawsTransf[n, ] - mean
      meandiff <- draws[n, ] - mean
      score[1:3, n] <- SigInv %*% meandiff
      
      product <- SigInv %*% meandiff %*% t(meandiff) %*% SigInv 
      
      scoreSig <- -SigInv + diag(diag(SigInv) )/ 2 + product - diag(diag(product))/2
    
      score[4, n] <- lambda[c(4, 7, 10)] %*% scoreSig[1, 1:3] + lambda[4] * scoreSig[1, 1] # U_11
      score[7, n] <- lambda[c(4, 7, 10)] %*%  scoreSig[2, 1:3] + lambda[7] * scoreSig[2, 2] # U_21
      score[8, n] <- lambda[c(8, 11) ] %*% scoreSig[2, 2:3] + lambda[8] * scoreSig[2, 2] # U_22
      score[10, n] <- - lambda[c(4, 7, 10)] %*%  scoreSig[3, 1:3] + lambda[10] * scoreSig[3, 3] # U_31
      score[11, n] <-  lambda[c(8, 11) ] %*% scoreSig[3, 2:3] + lambda[11] * scoreSig[3, 3] # U_32
      score[12, n] <- 2 * lambda[12] %*% scoreSig[3, 3] #U_33
    }
    
    w <- qDens / post
    LB[iter] <- mean(w * (log(post) - log(qDens)))
    
    gradient <- w * t(score) * (log(post) - log(qDens))
    gradient <- colMeans(gradient)
    gradientSq <- gradient^2
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
  #print(paste0('Iteration: ', iter-1, ' ELBO: ', meanLB))
  list(lambda = lambda, iter = iter-1, LB = LB[1:(iter-1)])
}

PMCMC <- function(y, lambda, iter, start, initial, stepsize = 0.01, 
                  targetAcceptance = 0.234, S = 200, suppressProgress = FALSE){
  dim <- length(start)
  saveDraws <- matrix(0, iter, dim+ 1)
  alpha <- -qnorm(targetAcceptance/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2 * 3.141598) * exp(alpha^2/2)/(2 * alpha) + 1/(dim * 0.234 * (1 - 0.234))
  draw <- start
  
  oldPF <- particleFilter(y, draw, lambda, initial, S)
  oldDens <- oldPF$dens
  state <- oldPF$state
  
  for (i in 1:iter) {
    
    
    if(!suppressProgress){
      if (i == 50) {
        startTime <- Sys.time()
      }
      else if (i == 150) {
        timePerIter <- (Sys.time() - startTime)/100
        class(timePerIter) <- "numeric"
        print(paste0("Estimated Finishing Time: ", Sys.time() + 
                       timePerIter * (iter - 150)))
        if (attr(timePerIter, "units") == "mins") {
          attr(timePerIter, "units") = "secs"
          timePerIter <- timePerIter * 60
        }
      }
      if (i %% 1000 == 0) {
        mins <- (iter - i) * timePerIter[1]/60
        if (mins > 180) {
          print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                       round(mins/60, 2), " hours."))
        }
        else if (mins > 1) {
          print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                       round(mins, 2), " minutes."))
        }
        else {
          print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                       ceiling(mins * 60), " seconds."))
        }
      }
    }
   
    logitDraw <- log(draw / (1 - draw))
    logitCan <- logitDraw + stepsize * rnorm(dim)
    candidate <- 1 / (1 + exp(-logitCan))
    
    canProp <- sum(log(1 / candidate + 1 / (1 - candidate)) + dnorm(logitCan, logitDraw, stepsize, log = TRUE))
    oldProp <- sum(log(1 / draw + 1 / (1 - draw)) + dnorm(logitDraw, logitCan, stepsize, log = TRUE))
    
    canPF <- particleFilter(y, candidate, lambda, initial, S)
    canDens <- canPF$dens
    ratio <- exp(canDens - oldDens + oldProp - canProp)
    if(is.na(ratio)){
      ratio <- 0
    }
    c <- stepsize * stepsizeCons
    if (runif(1) < ratio) {
      draw <- candidate
      state <- canPF$state
      oldDens <- canDens
      stepsize <- stepsize + c * (1 - targetAcceptance)/(28 + i)
    } else {
      stepsize <- stepsize - c * targetAcceptance/(28 + i)
    }
    saveDraws[i, ] <- c(draw, state)
  }
  draws <- data.frame(saveDraws)
  if (is.null(names(start))) {
    colnames(draws)[1:dim] <- paste0("V", 1:dim)
  }  else {
    colnames(draws)[1:dim] <- names(start)
  }
  colnames(draws)[dim+1] <- 'X[T]'
  draws
}

set.seed(id)
t1 <- 0.9
t2 <- 0.2

T <- 1000
Tseq <- c(0, seq(100, T, 100))

data <- drawSV(c(t1, t2), T+10)
y <- data$y
v <- data$v

priorLambda <- c(rep(0, 3), diag(10, 3))

for(t in 2:length(Tseq)){
  print(paste(id, t))
  
  yFull <- y[1:Tseq[t]]
  
  startTime <- Sys.time()
  
  if(t == 2){
    draw <- c('theta[1]' = 0.8, 'theta[2]' = 0.2)
  } else {
    draw <- colMeans(mcmc[,1:2])
  }
  mcmc <- PMCMC(y = yFull,
                lambda = priorLambda,
                iter = 30000, 
                start = draw,
                initial = TRUE,
                suppressProgress = TRUE)[20001:30000, ]
  
  mcmcTime <- Sys.time() - startTime
  if(attr(mcmcTime, 'units') == 'mins'){
    mcmcTime <- as.numeric(mcmcTime)
    mcmcTime <- mcmcTime * 60
  } else {
    mcmcTime <- as.numeric(mcmcTime)
  }
  
  
  if(t == 2){
    vmcmc <- mcmc
  } else {
    ySub <- y[(Tseq[t-1]+1):Tseq[t]]
    
    startTime <- Sys.time()
    
    margins <- apply(vmcmc, 2, function(x){
      dens <- density(x)
      approxfun(dens)
    })
    cdf <- apply(vmcmc, 2, function(x){
      y <- seq(1/length(x), 1, 1/length(x))
      y[rank(x)]
    })
    vine <- VineCopula::RVineStructureSelect(cdf[sample(1:10000, 1000), ])
    
    subset <- sample(1:10000, 100)
    draws <- vmcmc[subset, ]
  
    marginDensity <- apply(vapply(1:3, function(x) margins[[x]](draws[,x]), runif(100)), 1, prod)
    
    copDens <- VineCopula::RVinePDF(cdf[subset, ], vine)
    
    postDens <- marginDensity * copDens

    # Fit VB
    lambda <- c(colMeans(vmcmc), chol(cov(vmcmc)))
    
    vb <- mcmcVB(lambda, draws, postDens, alpha = 0.001, threshold = 0.0005)
    
    
    vmcmc <- PMCMC(y = ySub,
                   lambda = vb$lambda,
                   iter = 30000,
                   start = vb$lambda[1:2],
                   initial = FALSE,
                   suppressProgress = TRUE)[20001:30000, ]
    
    vmcmcTime <- Sys.time() - startTime
    if(attr(vmcmcTime, 'units') == 'mins'){
      vmcmcTime <- as.numeric(vmcmcTime)
      vmcmcTime <- vmcmcTime * 60
    } else {
      vmcmcTime <- as.numeric(vmcmcTime)
    }
  }
  
  # Forecast    
  
  if(t == 2){
  
    N <- 500
    subset <- sample(1:10000, N)
    fcSims <- matrix(0, N, 10)
    for(i in 1:N){
      fcSims[i, ] <- forecastSV(y[1:Tseq[t]], unlist(mcmc[subset[i], ]))
    }
    densFuns <- apply(fcSims, 2, function(x){
      xSub <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
      dens <- density(xSub)
      approxfun(dens)
    })
    
    ls <- vapply(1:10, function(x) log(densFuns[[x]](y[Tseq[t]+x])), runif(1))
    
    results <- tibble(ls = rep(ls, 2),
                      t = Tseq[t],
                      h = rep(1:10, 2),
                      method = rep(c('MCMC', 'V-MCMC'), rep(10, 2)),
                      time = mcmcTime)
    
  } else {
    
    subset <- sample(1:10000, N)
    fcSims <- matrix(0, N, 10)
    for(i in 1:N){
      fcSims[i, ] <- forecastSV(y[1:Tseq[t]], unlist(mcmc[subset[i], ]))
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
                      method = 'MCMC',
                      time = mcmcTime))
    
    
    for(i in 1:N){
      fcSims[i, ] <- forecastSV(y[1:Tseq[t]], unlist(vmcmc[subset[i], ]))
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
                      method = 'V-MCMC',
                      time = vmcmcTime))
  }
  
}

write.csv(results, paste0('VMCMC/rep', id, '.csv'), row.names = FALSE)
