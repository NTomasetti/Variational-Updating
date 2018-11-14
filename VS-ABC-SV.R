library(tidyverse)

eucdist <- function(y,z, w) {
  t(y-z) %*% w %*% (y-z)
}

abcVB <- function(lambda, z, draws, post, maxIter = 1000, threshold = 0.001, alpha = 0.01){

  draws <- as.matrix(draws)
  drawsTransf <- matrix(0, nrow(draws), 2)
  drawsTransf[,1] <- -log((1 - draws[,1]) / (draws[,1] - 0.5))
  drawsTransf[,2] <- -log((1 - 2 * draws[,2]) / (2 * draws[,2]))
  Jacobian <- apply(draws, 1, function(x){
    (1  / (x[1] - 0.5) + 1 / (1 - x[1])) * (1 / x[2] + 2 * x[2] / (1 - 2 * x[2]))
  })
  
  mix <- ncol(lambda)
  beta1 <- 0.9
  beta2 <- 0.99
  e <- 1e-8
  
  LB <- rep(0, maxIter)
  meanLB <- 5
  iter <- 1
  diff <- threshold + 1
  M <- V <- matrix(0, nrow(lambda), ncol(lambda))
  MZ <- VZ <- rep(0, mix)
  
  while(diff > threshold & iter <= maxIter){
    
    mean <- lambda[1:2, ]
    Sigma <- array(0, dim = c(2, 2, mix))
    SigInv <- Sigma
    
    for(m in 1:mix){
      U <- matrix(lambda[3:6, m], 2)
      Sigma[,,m] <- t(U) %*% U
      SigInv[,,m] <- solve(Sigma[,,m])
    }
    
    
    pi <- exp(z) / sum(exp(z))
    qComp <- matrix(0, nrow(draws), mix)
    for(m in 1:mix){
      qComp[,m] <- mvtnorm::dmvnorm(drawsTransf, mean[,m], Sigma[,,m])
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
    
    score <- array(0, dim = c(6, mix, nrow(draws)))
    for(m in 1:mix){
      for(n in 1:nrow(draws)){
        meandiff <- (drawsTransf[n, ] - mean[,m])
        score[1:2, m, n] <- SigInv[,,m] %*% meandiff
        
        product <- SigInv[,,m] %*% meandiff %*% t(meandiff) %*% SigInv[,,m] 
        
        scoreSig <- -SigInv[,,m] + diag(diag(SigInv[,,m]) )/ 2 + product - diag(diag(product))/2
        
        score[3, m, n] <- 2 * lambda[c(3, 5), m] %*% scoreSig[1, 1:2] # U_11
        score[5, m, n] <- 2 * lambda[c(3, 5), m] %*% scoreSig[2, 1:2] # U_21
        score[6, m, n] <- 2 * lambda[6, m] %*% scoreSig[2, 2] # U_22
        
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
    if(iter %% 100 == 0){
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('Iteration: ', iter-1, ' ELBO: ', meanLB))
  list(lambda = lambda, z = z, iter = iter-1, LB = LB[1:(iter-1)])
}


numericalDeriv <- function(lambda, z, draw, post){

  mix <- ncol(lambda)
  mean <- lambda[1:2, ]
  Sigma <- array(0, dim = c(2, 2, mix))
  SigInv <- Sigma
  
  for(m in 1:mix){
    U <- matrix(lambda[3:6, m], 2)
    Sigma[,,m] <- t(U) %*% U
    SigInv[,,m] <- solve(Sigma[,,m])
  }
  
  
  pi <- exp(z) / sum(exp(z))
  qComp <- rep(0, mix)
  for(m in 1:mix){
    qComp[m] <- mvtnorm::dmvnorm(draw, mean[,m], Sigma[,,m])
  }
  qDens <- c(qComp %*% pi)

  w <- qDens / post
  LB <- w * (log(post) - log(qDens))
  

  
  deriv <- matrix(0, nrow(lambda), ncol(lambda))
  dz <- rep(0, mix)
  for(n in 1:mix){
    
    z2 <- z
    z2[n] <- z2[n] + 1e-13
    pi <- exp(z2) / sum(exp(z2))
    qComp <- rep(0, mix)
    for(m in 1:mix){
      qComp[m] <- mvtnorm::dmvnorm(draw, mean[,m], Sigma[,,m])
    }
    qDens <- c(qComp %*% pi)
    
    w <- qDens / post
    LB2 <- w * (log(post) - log(qDens))
    
    dz[n] <- (LB2 - LB) / 1e-13
    
    for(elem in c(1, 2, 3, 5, 6)){
      
      l2 <- lambda
      l2[elem, n] <- l2[elem, n] + 1e-13
      
      mean2 <- l2[1:2, ]
      Sigma2 <- array(0, dim = c(2, 2, mix))
      SigInv2 <- Sigma
      
      for(m in 1:mix){
        U2 <- matrix(l2[3:6, m], 2)
        Sigma2[,,m] <- t(U2) %*% U2
        SigInv2[,,m] <- solve(Sigma2[,,m])
      }
      
      pi <- exp(z) / sum(exp(z))
      qComp <- rep(0, mix)
      for(m in 1:mix){
        qComp[m] <- mvtnorm::dmvnorm(draw, mean2[,m], Sigma2[,,m])
      }
      qDens <- c(qComp %*% pi)
      
      w <- qDens / post
      LB2 <- w * (log(post) - log(qDens))
      deriv[elem, n] <- (LB2 - LB) / 1e-13
    }
  }
  rbind(dz, deriv)
}

analyticalDeriv <- function(lambda, z, draw, post){
  mix <- ncol(lambda)
  mean <- lambda[1:2, ]
  Sigma <- array(0, dim = c(2, 2, mix))
  SigInv <- Sigma
  
  for(m in 1:mix){
    U <- matrix(lambda[3:6, m], 2)
    Sigma[,,m] <- t(U) %*% U
    SigInv[,,m] <- solve(Sigma[,,m])
  }
  
  
  pi <- exp(z) / sum(exp(z))
  qComp <- rep(0, mix)
  for(m in 1:mix){
    qComp[m] <- mvtnorm::dmvnorm(draw, mean[,m], Sigma[,,m])
  }
  qDens <- c(qComp %*% pi)
  
  scorePi <- qComp / qDens
  
  scoreZ <- rep(0, mix)
  denom <- sum(exp(z))^2
  
  for(i in 1:mix){
    for(j in 1:mix){
      if(i == j){
        scoreZ[i] <- scoreZ[i] + scorePi[j] * sum(exp(z[i] + z[-i])) / denom
      } else {
        scoreZ[i] <- scoreZ[i] - scorePi[j] * exp(z[i] + z[j]) / denom
      }
    }
  }
  
  score <- matrix(0, 6, mix)
  for(m in 1:mix){
    score[1:2, m] <- SigInv[,,m] %*% t(draw - mean[,m])
      
    scoreSig <- -SigInv[,,m] + diag(diag(SigInv[,,m]) )/ 2 + SigInv[,,m] %*% t(draw - mean[,m]) %*% t(t(draw - mean[,m])) %*% SigInv[,,m] -
      diag(diag(SigInv[,,m] %*% t(draw - mean[,m]) %*% t(t(draw - mean[,m])) %*% SigInv[,,m]))/2
    
    score[3, m] <- 2 * lambda[c(3, 5), m] %*% scoreSig[1, 1:2] # U_11
    score[5, m] <- 2 * lambda[c(3, 5), m] %*% scoreSig[2, 1:2] # U_21
    score[6, m] <- 2 * lambda[6, m] %*% scoreSig[2, 2] # U_22
    
    score[, m] <- score[, m] * pi[m] * qComp[m] / qDens # Convert from dlog N_i / dlam to dlog N_sum /dlam
  }
  
  w <- qDens / post

  gradient <- w * score * (log(post) - log(qDens))
  
  gradZ <- w * scoreZ * (log(post) - log(qDens))
  rbind(gradZ, gradient)
}

garchScore <- function(y, beta){
  T <- length(y)
  deriv <- rep(0, 3)
  
  # V_1
  v <- beta[1] / (1 - beta[2] - beta[3])
  # d V_1 / d Beta
  dvdb <- c(1 / (1 - beta[2] - beta[3]),  beta[1] / (1 - beta[2] - beta[3])^2,  beta[1] / (1 - beta[2] - beta[3])^2)
  # d Y_1 / dBeta
  deriv <- dvdb * (-1 / (2 * v) + y[1]^2 / (2 * v^2))
  
  for(t in 2:T){
    v <- c(v, beta[1] + beta[2] * y[t-1]^2 + beta[3] * v[t-1])
    dvdb <- c(1, y[t-1]^2, v[t-1]) + beta[3] * dvdb
    deriv <- deriv  + dvdb * (-1 / (2 * v[t]) + y[t]^2 / (2 * v[t]^2))
  }
  deriv
}

qDens <- function(lambda, mix, theta){
  weights <- lambda[seq(1, dim*(mix-1)+1, dim)]
  pi <- exp(weights) / sum(exp(weights))
  density <- 0
  for(m in 1:mix){
    mean <- lambda[dim*(mix-1) + 1 + 1:ncol(accept)]
    U <- matrix(lambda[dim*(mix-1) + (ncol(accept)+2):dim], ncol(accept))
    Sig <- t(U) %*% U
    density <- density + pi[m] * mvtnorm::dmvnorm(theta, mean, Sig)
  }
  density
}

qMarginal <- function(lambda, z, mix, supports, names){
  n <- ncol(supports)
  dim <- n * (n + 1)
  pi <- exp(z) / sum(exp(z))
  mean <- lambda[1:n, ]
  Sig <- array(0, dim = c(n, n, mix))
  for(m in 1:mix){
    U <- matrix(lambda[(n+1):dim, m], n)
    Sig[,,m]<- t(U) %*% U
  }
  density <- NULL
  for(i in 1:n){
    dens <- rep(0, nrow(supports))
    for(m in 1:mix){
      dens <- dens + pi[m] * dnorm(supports[,i], mean[i, m], sqrt(Sig[i, i, m]))
    }
    density <- rbind(density,
                     tibble(var = names[i],
                            support = supports[,i],
                            dens = dens))
  }
  density
}

qBivariate <- function(lambda, z, mix, supports, names){
  n <- ncol(supports)
  dim <- n * (n + 1)
  pi <- exp(z) / sum(exp(z))
  mean <- lambda[1:n, ]
  Sig <- array(0, dim = c(n, n, mix))
  for(m in 1:mix){
    U <- matrix(lambda[(n+1):dim, m], n)
    Sig[,,m]<- t(U) %*% U
  }
  density <- NULL
  index <- combn(1:n, 2)
  
  for(i in 1:ncol(index)){
    grid <- expand.grid(supports[, index[1, i]], supports[, index[2, i]])
    dens <- rep(0, nrow(grid))
    for(m in 1:mix){
      dens <- dens + pi[m] * mvtnorm::dmvnorm(grid,
                                              mean[c(index[1, i], index[2, i]), m],
                                              matrix(c(Sig[index[1, i], index[1, i], m],
                                                       Sig[index[1, i], index[2, i], m],
                                                       Sig[index[1, i], index[2, i], m],
                                                       Sig[index[2, i], index[2, i], m]), 2))
    }
    density <- rbind(density,
                     tibble(var1 = names[index[1, i]],
                            var2 = names[index[2, i]],
                            support1 = grid[, 1],
                            support2 = grid[, 2],
                            dens = dens))
  }
  density
}


set.seed(20)
t1 <- 0.6
t2 <- 0.1

T <- 500

v0 <- rnorm(1, t1, sqrt(t2 / (1 - t1^2)))
v <- t1 * v0 + rnorm(1, 0, sqrt(t2))
y <- rnorm(1, 0, sqrt(exp(v)))

for(t in 2:T){
  v[t] <- t1 * v[t-1] + rnorm(1, 0, sqrt(t2))
  y[t] <- sqrt(exp(v[t])) * rnorm(1, 0, 1)
}

#garchMod <- rugarch::ugarchspec(mean.model = list(armaOrder = c(0, 0), 
#                                                  include.mean = FALSE))

#garch <- rugarch::ugarchfit(spec = garchMod, data = y)
#yss <- colSums(garch@fit$scores)

garch <- fGarch::garchFit(data = y, include.mean = FALSE)

yss <- garchScore(y, garch@fit$coef)

reps <- 50000
keep <- 0.01

theta <- matrix(0, reps, 2)
theta[,1] <- runif(reps, 0.5, .99)
theta[,2] <- runif(reps, 0.05, 0.5)

zss <- matrix(0, reps, 3)
#loop through simulation of data
for(iter in 1:reps){
  
  if(iter %% 500 == 0){
    print(iter)
  }
  
  vz <- rnorm(1, t1, sqrt(theta[iter, 2] / (1 - theta[iter, 1]^2)))
  vz <- t1 * vz + rnorm(1, 0, sqrt(theta[iter, 2]))
  z <- rnorm(1, 0, sqrt(exp(vz)))
  
  for(t in 2:T){
    vz[t] <- theta[iter, 1] * vz[t-1] + rnorm(1, 0, sqrt(theta[iter, 2]))
    z[t] <- sqrt(exp(vz[t])) * rnorm(1, 0, 1)
  }
  zss[iter, ] <- garchScore(z, garch@fit$par)
}

varZ <- apply(zss, 2, var, na.rm = TRUE)
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

ABC <- as.tibble(cbind(theta, distance))
names(ABC) = c('t1', 't2', 'dist')

ABC %>%
  filter(dist < quantile(dist, keep, na.rm = TRUE)) %>% 
  select(-dist) -> accept

GGally::ggpairs(accept) + theme_bw()

grid <-expand.grid(t1 = seq(0.5001, 0.999, length.out = 100), t2 = seq(0.05, 0.4999, length.out = 100))


accept %>%
  as.matrix() %>%
  ks::kde(eval.points = grid) -> postDensity

grid$post <- c(postDensity$estimate)

grid %>%
  gather(var, support, -post) %>%
  group_by(var, support) %>%
  summarise(post = sum(post)) %>%
  ungroup() %>%
  group_by(var) %>%
  mutate(cons = sum((support - lag(support))*post, na.rm = TRUE),
         post = post / cons) %>%
  ggplot() + geom_line(aes(support, post), colour = 'red') + 
  geom_density(data = accept %>% gather(var, draw), aes(draw)) + facet_wrap(~var, scales = 'free')

normProb <- grid$post
normProb[normProb <= 0] <- 0
normProb <- normProb / sum(normProb)

subset <- sample(1:nrow(grid), 70, prob = normProb)
draws <- grid[subset, 1:2] 

postDens <- grid[subset, ]$post 

draws %>%
  gather(var, draw) %>%
  ggplot() + geom_histogram(aes(draw)) + facet_wrap(~var, scales = 'free')


drawsTransf <- matrix(0, nrow(draws), 2)
drawsTransf[,1] <- -log((1 - draws[,1]) / (draws[,1] - 0.5))
drawsTransf[,2] <- -log((1 - 2 * draws[,2]) / (2 * draws[,2]))

mix <- 3
lambda <- matrix(0, 6, mix)
lambda[,1] <- c(colMeans(drawsTransf), chol(cov(drawsTransf)))
for(m in 2:mix){
  lambda[,m] <- c(colMeans(drawsTransf) + c(0.5, 0.5) * (-1)^m + rnorm(2, 0, 0.15), chol(cov(drawsTransf)))
}
z <- rep(0, 3)
dim <- nrow(lambda) / mix

vb <- abcVB(lambda, z, draws, postDens, alpha = 0.002, threshold = 0.0001)

ggplot() + geom_line(aes(1:vb$iter, vb$LB))

pi <- exp(vb$z) / sum(exp(vb$z))
mean <- vb$lambda[1:2, ]
Sig <- array(0, dim = c(2, 2, mix))
for(m in 1:mix){
  U <- matrix(vb$lambda[3:6, m], 2)
  Sig[,,m]<- t(U) %*% U
}

reps <- 2000
keep <- 0.05
theta <- matrix(0, reps, 2)

for(iter in 1:reps){
  u <- runif(1)
  component <- min(which(cumsum(pi) > u))
  theta[iter, ] <- mvtnorm::rmvnorm(1, mean[,component], Sig[,,component])
}

theta[,1] <- 0.5 / (1 + exp(-theta[,1])) + 0.5
theta[,2] <- 0.5 / (1 + exp(-theta[,2]))

theta %>%
  as.data.frame() %>%
  rename(t1 = V1, t2 = V2) %>%
  gather(var, value) %>%
  ggplot() + geom_density(aes(value), colour = 'red') + 
  geom_density(data = accept %>% gather(var, value), aes(value)) + facet_wrap(~var, scales = 'free')




T2 <- 100
v2 <- t1 * v[T] + rnorm(1, 0, sqrt(t2))
y2 <- rnorm(1, 0, sqrt(exp(v2)))

for(t in 2:T2){
  v2[t] <- t1 * v2[t-1] + rnorm(1, 0, sqrt(t2))
  y2[t] <- sqrt(exp(v2[t])) * rnorm(1, 0, 1)
}

yFull <- c(y, y2)
garchMod <- rugarch::ugarchspec(mean.model = list(armaOrder = c(0, 0), 
                                                  include.mean = FALSE))

garch <- rugarch::ugarchfit(spec = garchMod, data = yFull)
yss <- colSums(garch@fit$scores)
omega <- garch@fit$coef[1]
alpha <- garch@fit$coef[2]
beta <- garch@fit$coef[3]

reps <- 10000
keep <- 0.02

theta <- matrix(0, reps, 2)
theta[,1] <- runif(reps, 0.5, .99)
theta[,2] <- runif(reps, 0.05, 0.5)

zss <- matrix(0, reps, 3)
#loop through simulation of data
for(iter in 1522:reps){
  
  if(iter %% 500 == 0){
    print(iter)
  }
  
  vz <- rnorm(1, t1, sqrt(theta[iter, 2] / (1 - theta[iter, 1]^2)))
  vz <- t1 * vz + rnorm(1, 0, sqrt(theta[iter, 2]))
  z <- rnorm(1, 0, sqrt(exp(vz)))
  
  for(t in 2:(T+T2)){
    vz[t] <- theta[iter, 1] * vz[t-1] + rnorm(1, 0, sqrt(theta[iter, 2]))
    z[t] <- sqrt(exp(vz[t])) * rnorm(1, 0, 1)
  }
  garchZ <- rugarch::ugarchfit(spec = garchMod, data = z, fixed.pars = list(omega = omega, alpha1 = alpha, beta1 = beta))
  zss[iter, ] <- colSums(garchZ@fit$scores)
}

zss[which(rowSums(zss) == 0), ] <- NA

varZ <- apply(zss, 2, var, na.rm = TRUE)
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

ABC2 <- as.tibble(cbind(theta, distance))
names(ABC2) = c('theta[1]', 'theta[2]', 'dist')

ABC2 %>%
  filter(dist < quantile(dist, keep, na.rm = TRUE)) %>%
  select(-dist) -> accept2


n <- ncol(accept)
dim <- 1 + n * (n + 1)
weights <- lambda[seq(1, dim*(mix-1)+1, dim)]
pi <- exp(weights) / sum(exp(weights))
mean <- matrix(0, n, mix)
Sig <- array(0, dim = c(n, n, mix))
for(m in 1:mix){
  mean[,m] <- lambda[dim*(m-1) + 1 + 1:ncol(accept)]
  U <- matrix(lambda[dim*(m-1) + (ncol(accept)+2):dim], ncol(accept))
  Sig[,,m]<- t(U) %*% U
}

reps <- 2000
keep <- 0.05
theta <- matrix(0, reps, 2)

for(iter in 1:reps){
  okay <- FALSE
  while(!okay){
    u <- runif(1)
    component <- min(which(cumsum(pi) > u))
    draw <- mvtnorm::rmvnorm(1, mean[,component], Sig[,,component])
    if(draw[1] > 0.5 & draw[1] < 0.99 & draw[2] > 0.05 & draw[2] < 0.5){
      okay <- TRUE
    }
  }
  theta[iter, ] <- draw
}

zss <- matrix(0, reps, 3)
#loop through simulation of data
for(iter in 1:reps){
  if(iter %% 500 == 0){
    print(iter)
  }
  
  vz <- t1 * v[T] + rnorm(1, 0, sqrt(t2))
  z <- rnorm(1, 0, sqrt(exp(vz)))
  
  for(t in 2:T2){
    vz[t] <- theta[iter, 1] * vz[t-1] + rnorm(1, 0, sqrt(theta[iter, 2]))
    z[t] <- sqrt(exp(vz[t])) * rnorm(1, 0, 1)
  }
  garchZ <- rugarch::ugarchfit(spec = garchMod, data = z, fixed.pars = list(omega = omega, alpha1 = alpha, beta1 = beta))
  zss[iter, ] <- colSums(garchZ@fit$scores)
  
  
}
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

ABCVB <- as.tibble(cbind(theta, distance))
names(ABCVB) = c('theta[1]', 'theta[2]', 'dist')

ABCVB %>%
  filter(dist < quantile(dist, keep)) %>%
  select(-dist) -> acceptVB


accept2 %>%
  mutate(method = 'ABC') %>%
  rbind(acceptVB %>%
          mutate(method = 'UV-ABC')) %>%
  gather(var, draw, -method) %>%
  ggplot() + geom_density(aes(draw, colour = method)) + 
  geom_vline(data = tibble(var = c('theta[1]', 'theta[2]'), actual = c(t1, t2)), aes(xintercept = actual)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed) + theme_bw() + 
  labs(x = NULL, y = NULL) + theme(legend.position = 'bottom')


qSeq <- seq(0, 0.99, 0.01)

tibble(q = qSeq,
       ABC = quantile(ABC2$dist, qSeq, na.rm = TRUE),
       `UV-ABC` = quantile(ABCVB$dist, qSeq, na.rm = TRUE)) %>% 
  gather(method, quantile, -q) %>%
  ggplot() + geom_point(aes(quantile, q, colour = method)) + 
  labs(x = 'Distance', y = 'CDF') + 
  theme_bw()


quantile(ABC2$dist, 0.01)
sum(ABCVB$dist <= quantile(ABC2$dist, 0.01)) / nrow(ABCVB)  
