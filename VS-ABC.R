library(tidyverse)
Rcpp::sourceCpp('VSABC.cpp')
set.seed(5)
mu <- 2
sigma <- 1
T <- 50

y <- rnorm(T, mu, sigma)
ssY <- c(mean(y), sd(y))
keep <- 0.0025

reps <- 100000
ABC <- list()

for(i in 1:reps){
  if(i %% 5000 == 0){
    print(i)
  }
  draw <- c(rnorm(1, 0, 10), rnorm(1, 0, 1))
  z <- rnorm(T, draw[1], exp(draw[2]))
  ssZ <- c(mean(z), sd(z))
  ABC[[i]] <- tibble(mu = draw[1],
                     logsd = draw[2],
                     ssZ = ssZ,
                     iter = i)
}

ABC <- bind_rows(ABC)

ABC %>%
  group_by(iter) %>%
  summarise(dist = sqrt((ssY[1] - ssZ[1])^2 + (ssY[2] - ssZ[2])^2)) %>%
  arrange(dist) %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  filter(order <= reps * keep) %>%
  left_join(ABC) -> accept

accept %>%
  dplyr::select(mu, logsd) %>%
  gather(var, draw) %>%
  ggplot() + geom_density(aes(draw)) + facet_wrap(~var, scales = 'free')


accept %>%
  dplyr::select(mu, logsd) %>%
  as.matrix() %>%
  ks::kde() -> postDensity

gridDensity <- expand.grid(mu = postDensity$eval.points[[1]],
                           logsd = postDensity$eval.points[[2]])
gridDensity$post <- c(postDensity$estimate)
gridDensity$post[gridDensity$post <= 0] <- 0

ggplot(gridDensity) + geom_tile(aes(mu, logsd, fill = post)) + 
  geom_density2d(data = accept, aes(mu, logsd))

normProb <- gridDensity$post
normProb[normProb <= 0] <- 0
normProb <- normProb / sum(normProb)

subset <- sample(1:nrow(gridDensity), 100, prob = normProb)
draws <- gridDensity[subset, c('mu', 'logsd')] 

postDens <- gridDensity[subset, ]$post 

ggplot(gridDensity) + geom_tile(aes(mu, logsd, fill = post)) + 
  geom_density2d(data = draws, aes(mu, logsd))


mix <- 2

maxIter <- 1000
threshold <- 0.0005
diff <- threshold + 1
LB <- rep(0, maxIter)
alpha <- 0.02
beta1 <- 0.9
beta2 <- 0.99
meanLB <- 5
iter <- 1
e <- 1e-8

lambda <- c(0, mean(accept$mu), mean(accept$logsd), t(chol(cov(cbind(accept$mu, accept$logsd)))))
for(m in 2:mix){
  lambda <- c(lambda, 0, mean(accept$mu) + 0.1 * m, mean(accept$logsd) + 0.1 * m, t(chol(cov(cbind(accept$mu, accept$logsd)))))
}

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
    print(paste0('Iteration: ', iter, ' ELBO: ', meanLB, ', Last Difference: ', diff))
  }
  iter <- iter + 1
}
print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))


qDens <- function(lambda, mix, mu, sd){
  weights <- lambda[seq(1, 7*(mix-1)+1, 7)]
  pi <- exp(weights) / sum(exp(weights))
  density <- 0
  for(m in 1:mix){
    mean <- lambda[7*(mix-1)+2:3]
    U <- matrix(lambda[7*(mix-1) + 4:7], 2)
    Sig <- t(U) %*% U
    density <- density + pi[m] * mvtnorm::dmvnorm(c(mu, sd), mean, Sig)
  }
  density
}

gridDensity %>%
  group_by(mu, logsd) %>%
  mutate(VBpost = qDens(lambda, mix, mu, logsd)) -> gridDensity
  
gridDensity %>%
  ggplot() + geom_tile(aes(mu, logsd, fill = VBpost)) + 
  geom_density2d(data = accept, aes(mu, logsd)) 

thetaSABC <- as.matrix(accept[sample(1:nrow(accept), 50000, replace = TRUE), 4:5])


y2 <- rnorm(T, mu, sigma)
yFull <- c(y, y2)
ssY <- c(mean(yFull), sd(yFull))
keep <- 0.0025

reps <- 100000
ABC <- list()

for(i in 1:reps){
  if(i %% 5000 == 0){
    print(i)
  }
  draw <- c(rnorm(1, 0, 10), rnorm(1, 0, 1))
  z <- rnorm(2*T, draw[1], exp(draw[2]))
  ssZ <- c(mean(z), sd(z))
  ABC[[i]] <-  tibble(mu = draw[1],
                      logsd = draw[2],
                      ssZ = ssZ,
                      iter = i)
}

ABC <- bind_rows(ABC)

ABC %>%
  group_by(iter) %>%
  summarise(dist = sqrt((ssY[1] - ssZ[1])^2 + (ssY[2] - ssZ[2])^2)) %>%
  arrange(dist) %>%
  left_join(ABC) -> ABC

ABC %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  filter(order <= reps * keep) -> accept

accept %>%
  dplyr::select(mu, logsd) %>%
  gather(var, draw) %>%
  ggplot() + geom_density(aes(draw)) + facet_wrap(~var, scales = 'free')

mean <- matrix(0, 2, mix)
Sigma <- array(0, dim = c(2, 2, mix))
z <- lambda[seq(1, 7 * (mix-1) + 1, 7)]
pi <- exp(z) / sum(exp(z))

for(m in 1:mix){
  mean[,m] <- lambda[7*(mix-1) + 2:3]
  U <- matrix(lambda[7 * (mix -1) + 4:7], 2)
  Sigma[,,m] <- t(U) %*% U
}

ssY <- c(mean(y2), sd(y2))
reps <- 50000
keep <- 0.05
ABCVB <- list()
for(i in 1:reps){
  if(i %% 5000 == 0){
    print(i)
  }
  u <- runif(1)
  component <- min(which(cumsum(pi) > u))
  draw <- mvtnorm::rmvnorm(1, mean[,component], Sigma[,,component])
  z <- rnorm(T, draw[1], exp(draw[2]))
  ssZ <- c(mean(z), sd(z))
  ABCVB[[i]] <- tibble(mu = draw[1],
                        logsd = draw[2],
                        ssZ = ssZ,
                        iter = i)
  
}

ABCVB <- bind_rows(ABCVB)

ABCVB %>%
  group_by(iter) %>%
  summarise(dist = sqrt((ssY[1] - ssZ[1])^2 + (ssY[2] - ssZ[2])^2)) %>%
  arrange(dist) %>%
  left_join(ABCVB) -> ABCVB
  
ABCVB %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  filter(order <= reps * keep) -> acceptVB

SABC <- list()
for(i in 1:reps){
  if(i %% 5000 == 0){
    print(i)
  }
  u <- runif(1)
  draw <- thetaSABC[i, ]
  z <- rnorm(T, draw[1], exp(draw[2]))
  ssZ <- c(mean(z), sd(z))
  SABC[[i]] <- tibble(mu = draw[1],
                       logsd = draw[2],
                       ssZ = ssZ,
                       iter = i)
  
}

SABC <- bind_rows(SABC)

SABC %>%
  group_by(iter) %>%
  summarise(dist = sqrt((ssY[1] - ssZ[1])^2 + (ssY[2] - ssZ[2])^2)) %>%
  arrange(dist) %>%
  left_join(SABC) -> SABC

SABC %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  filter(order <= reps * keep) -> acceptSABC




accept %>%
  mutate(method = 'ABC') %>%
  rbind(acceptVB %>%
          mutate(method = 'ABC-VB')) %>%
  dplyr::select(method, mu, logsd) %>%
  gather(var, draw, -method) %>%
  ggplot() + geom_density(aes(draw, colour = method)) + facet_wrap(~var, scales = 'free')


qSeq <- seq(0, 0.995, 0.005)

tibble(q = qSeq,
       ABC = quantile(ABC$dist, qSeq),
       `VS-ABC` = quantile(ABCVB$dist, qSeq),
       `S-ABC` = quantile(SABC$dist, qSeq)) %>% 
  gather(method, quantile, -q) %>%
  ggplot() + geom_point(aes(quantile, q, colour = method)) + 
  labs(x = 'Distance', y = 'CDF', colour = 'Method') + 
  theme_bw() + theme(legend.position = 'bottom')
sum(ABCVB$dist <= quantile(ABC$dist, 0.01)) / nrow(ABCVB)  
sum(SABC$dist <= quantile(ABC$dist, 0.01)) / nrow(SABC)  


muSeq <- seq(1.5, 2.5, length.out = 250)
sigSeq <- seq(0.6, 1.5, length.out = 250)
grid <- expand.grid(mu = muSeq, sig = sigSeq)
grid %>%
  group_by(mu, sig) %>%
  mutate(dens = prod(dnorm(c(y, y2), mu, sig)) * dnorm(mu, 0, 10) * dlnorm(sig, 0, 1)) -> grid

grid %>% 
  ungroup() %>%
  group_by(mu) %>%
  summarise(dens = sum(dens)) %>%
  ungroup() %>%
  mutate(dens = dens / sum(dens) / (muSeq[2] - muSeq[1])) %>%
  rename(support = mu) -> muMarginal

grid %>% 
  ungroup() %>%
  group_by(sig) %>%
  summarise(dens = sum(dens)) %>%
  ungroup() %>%
  mutate(dens = dens / sum(dens) / (sigSeq[2] - sigSeq[1]) ) %>%
  rename(support = sig) %>%
  rbind(muMarginal) %>%
  mutate(var = rep(c('sigma', 'mu'), rep(250, 2)),
         method = 'Exact') -> marginals

accept %>%
  mutate(method = 'ABC') %>%
  rbind(acceptVB %>%
          mutate(method = 'VS-ABC'),
        acceptSABC %>%
          mutate(method = 'S-ABC')) %>%
  mutate(sd = exp(logsd)) %>%
  select(method, mu, sd) %>%
  rename(sigma = sd) %>%
  gather(var, draw, -method) %>%
  ggplot() + geom_density(aes(draw, colour = method)) + 
  geom_line(data = marginals, aes(support, dens, colour = method)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed) +
  labs(x = NULL, y = NULL, colour = 'Method') + 
  theme_bw() + theme(legend.position = 'bottom')




