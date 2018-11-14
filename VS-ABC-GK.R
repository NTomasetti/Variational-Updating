library(tidyverse)
Rcpp::sourceCpp('VSABC.cpp')

eucdist <- function(y,z, w) {
  t(y-z) %*% w %*% (y-z)
}

drawGK <-function(theta,T){
  
  a <- theta[1]
  b <- theta[2]
  g <- theta[3]
  k <- theta[4]
  
  xy <- rnorm(T)
  y <- a+b*(1+.8*((1-exp(-g*xy))/(1+exp(-g*xy))))*((1+xy^2)^k)*xy
  y
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

qMarginal <- function(lambda, mix, supports, names){
  n <- ncol(supports)
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

qBivariate <- function(lambda, mix, supports, names){
  n <- ncol(supports)
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
a <- 3 #mean
b <- 1 #scale
g <- 2 #skewness
k <- 0.5 # kurtosis

T <- 200
reps <- 50000
keep <- 0.01

y <- drawGK(c(a, b, g, k), T)

s1=function(x, qr) (1/T*sum(x))
s3=function(x, qr) 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)
s4=function(x, qr) (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5
s1r = function(x, qr) qr[4]
s2r = function(x, qr) qr[6]-qr[2]
s3r = function(x, qr) (qr[6]+qr[2]-2*qr[4])/(qr[6]-qr[2])
s4r = function(x, qr) (qr[7]-qr[5]+qr[3]-qr[1])/(qr[6]-qr[2])


sumStats <- plyr::each(s1r, s2 = s2r, s3 = s3, s4 = s4, s5 = s3r, s6 = s4r)
qrY <- quantile(y, c(.08,.25,.36,.5,.6,.75,.875))
yss <- sumStats(y, qrY)

theta <- matrix(runif(4*reps, 0, 10), nrow=reps)
zss <- matrix(0, reps, 6)
#loop through simulation of data
for(iter in 1:reps){
  z <- drawGK(theta[iter, ], T)
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

ABC %>%
  filter(dist < quantile(dist, keep)) %>%
  dplyr::select(a, b, g, k) -> accept

GGally::ggpairs(accept) + theme_bw()

grid <- do.call(expand.grid, lapply(accept, quantile, prob=seq(0.01, 0.99, length.out = 15)))

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
  ggplot() + geom_line(aes(support, post)) + 
  geom_density(data = accept %>% gather(var, draw), aes(draw), colour = 'red') + facet_wrap(~var, scales = 'free')

normProb <- grid$post
normProb[normProb <= 0] <- 0
normProb <- normProb / sum(normProb)

subset <- sample(1:nrow(grid), 100, prob = normProb)
draws <- grid[subset, 1:4] 

postDens <- grid[subset, ]$post 

draws %>%
  gather(var, draw) %>%
  ggplot() + geom_histogram(aes(draw)) + facet_wrap(~var, scales = 'free')



mix <- 5

maxIter <- 1000
threshold <- 0.00001
alpha <- 0.01
beta1 <- 0.9
beta2 <- 0.99
e <- 1e-8

LB <- rep(0, maxIter)
meanLB <- 5
iter <- 1
diff <- threshold + 1
lambda <- c(0, colMeans(accept), t(chol(cov(accept))))
for(m in 2:mix){
  lambda <- c(lambda, 0, colMeans(accept) + c(0.5, 2.5, 0.5, 0.3) * (-1)^m + rnorm(4, 0, 0.1), t(chol(cov(accept))))
}
dim <- nrow(lambda) / mix


lambda <- matrix(lambda, ncol = 1)
M <- V <- rep(0, nrow(lambda))

while(diff > threshold & iter <= maxIter){
  
  score <- matrix(0, nrow(draws), nrow(lambda))
  logq <- rep(0, nrow(draws))
  for(i in 1:nrow(draws)){
    deriv <- mixNormScore(unlist(draws[i, ]), lambda, mix)
    score[i,] <- deriv$grad
    logq[i] <- deriv$val
  }
  
  w <- exp(logq) / postDens
  LB[iter] <- mean(w * (log(postDens) - logq))
  
  gradient <- w * score * (log(postDens) - logq)
  gradient <- apply(gradient, 2, mean)
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
  if(iter > 1){
    lambda <- lambda + update
  }
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
print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
ggplot() + geom_line(aes(1:(iter-1), LB[1:(iter-1)]))


gridVB <- cbind(a = seq(0, 6, length.out = 100),
                b = seq(0, 10, length.out = 100),
                g = seq(0, 10, length.out = 100),
                k = seq(0, 1.6, length.out = 100))
  

densVB <- qMarginal(lambda, mix, gridVB, colnames(gridVB))

ggplot() + geom_line(data = densVB, aes(support, dens), colour = 'red') + 
  geom_density(data = accept %>% gather(var, draw), aes(draw)) + 
  facet_wrap(~var, scales = 'free') + theme_bw() + labs(x = NULL, y = NULL) -> margins

bivVB <- qBivariate(lambda, mix, gridVB, colnames(gridVB))

pairs <- combn(1:4, 2)

plots <- apply(pairs, 2, function(x){
  samples <- accept[,x]
  vars <- colnames(samples)
  names(samples) <- c('V1', 'V2')
  
  bivVB %>%
    filter(var1 == vars[1] & var2 == vars[2]) %>%
    ggplot()  + geom_tile(aes(support1, support2, fill = dens)) + 
    geom_density2d(data = samples, aes(V1, V2), colour = 'black') + 
    labs(x = vars[1], y = vars[2]) + 
    theme_bw() +
    scale_fill_gradient(low = scales::muted('darkred'), high = 'white') + 
    theme(legend.position = 'none') -> p
})

gridExtra::grid.arrange(margins, plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
                        layout_matrix = matrix(c(1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7), byrow = TRUE, ncol = 3))








y2 <- drawGK(c(a, b, g, k), T)

yFull <- c(y, y2)
qrY <- quantile(yFull, c(.08,.25,.36,.5,.6,.75,.875))
yss <- sumStats(yFull, qrY)

theta <- matrix(runif(4*reps, 0, 10), nrow=reps)
zss <- matrix(0, reps, 6)
#loop through simulation of data
for(iter in 1:reps){
  z <- drawGK(theta[iter, ], 2*T)
  qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
  zss[iter, ] <- sumStats(z, qrZ)
}
varZ <- apply(zss, 2, var)
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

ABC2 <- as.tibble(cbind(theta, distance))
names(ABC2) = c('a', 'b', 'g', 'k', 'dist')

ABC2 %>%
  filter(dist < quantile(dist, keep)) %>%
  dplyr::select(a, b, g, k) -> accept2


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

qrY <- quantile(y2, c(.08,.25,.36,.5,.6,.75,.875))
yss <- sumStats(y2, qrY)

reps <- 5000
theta <- matrix(0, reps, 4)

for(iter in 1:reps){
  okay <- FALSE
  while(!okay){
    u <- runif(1)
    component <- min(which(cumsum(pi) > u))
    draw <- mvtnorm::rmvnorm(1, mean[,component], Sig[,,component])
    if(all(draw > 0) & all(draw < 10)){
      okay <- TRUE
    }
  }
  theta[iter, ] <- draw
}

zss <- matrix(0, reps, 6)
#loop through simulation of data
for(iter in 1:reps){
  z <- drawGK(theta[iter, ], T)
  qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
  zss[iter, ] <- sumStats(z, qrZ)
}
#varZ <- apply(zss, 2, var)
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

ABCVB <- as.tibble(cbind(theta, distance))
names(ABCVB) = c('a', 'b', 'g', 'k', 'dist')

ABCVB %>%
  filter(dist < quantile(dist, keep * 10)) %>%
  dplyr::select(a, b, g, k) -> acceptVB


theta <- as.matrix(accept[sample(1:nrow(accept), reps, replace = TRUE), ])
zss <- matrix(0, reps, 6)
#loop through simulation of data
for(iter in 1:reps){
  z <- drawGK(theta[iter, ], T)
  qrZ <- quantile(z, c(.08,.25,.36,.5,.6,.75,.875))
  zss[iter, ] <- sumStats(z, qrZ)
}
#varZ <- apply(zss, 2, var)
distance <- rep(0, reps)
for(iter in 1:reps){
  distance[iter] <- eucdist(yss, zss[iter, ], diag(1/varZ))
}

SABC <- as.tibble(cbind(theta, distance))
names(SABC) = c('a', 'b', 'g', 'k', 'dist')

SABC %>%
  filter(dist < quantile(dist, keep * 10)) %>%
  dplyr::select(a, b, g, k) -> acceptSABC



accept2 %>%
  mutate(method = 'ABC') %>%
  rbind(acceptVB %>%
          mutate(method = 'VS-ABC'),
        acceptSABC %>%
          mutate(method = 'S-ABC')) %>%
  gather(var, draw, -method) %>%
  ggplot() + geom_density(aes(draw, colour = method)) + 
  geom_vline(data = tibble(var = c('a', 'b', 'g', 'k'), actual = c(a, b, g, k)), aes(xintercept = actual)) + 
  facet_wrap(~var, scales = 'free') + theme_bw() + 
  labs(x = NULL, y = NULL, colour = 'Method') + theme(legend.position = 'bottom')


qSeq <- seq(0, 0.99, 0.01)

tibble(q = qSeq,
       ABC = quantile(ABC2$dist, qSeq),
       `S-ABC` = quantile(SABC$dist, qSeq), 
       `VS-ABC` = quantile(ABCVB$dist, qSeq)) %>%
  gather(method, quantile, -q) %>%
  ggplot() + geom_point(aes(quantile, q, colour = method)) + 
  labs(x = 'Distance', y = 'CDF', colour = 'Method') + 
  theme_bw()


quantile(ABC2$dist, 0.01)
sum(ABCVB$dist <= quantile(ABC2$dist, 0.01)) / nrow(ABCVB)  
sum(SABC$dist <= quantile(ABC2$dist, 0.01)) / nrow(SABC)  


