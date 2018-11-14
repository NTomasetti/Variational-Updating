rm(list = ls())
library(tibble)

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


results <- list()
for(id in 1:200){
  set.seed(id)
  t1 <- 0.9
  t2 <- 0.2
  
  T <- 2000
  Tseq <- c(0, seq(100, T, 50))
  reps <- 50000
  keep <- 0.01
  
  data <- drawSV(c(t1, t2), T+10)
  y <- data$y

  res <- tibble()
  
  for(t in 2:length(Tseq)){
    yFull <- y[1:Tseq[t]]
    mean <- mean(yFull)
    sd <- sd(yFull)
    ls <- dnorm(y[Tseq[t]+1:10], mean, sd, log = TRUE)
    res <- rbind(res,
                 tibble(ls = ls,
                        t = Tseq[t],
                        h = 1:10,
                        method = 'Baseline',
                        id = id))
  }
  results[[id]] <- res
}
results <- dplyr::bind_rows(results)

# sv from vs resultsR.

sv %>%
  group_by(t, method) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  ungroup() %>%
  rbind(results %>%
          group_by(t, method) %>%
          summarise(ls = mean(ls, na.rm = TRUE)) %>%
          ungroup()) %>%
  spread(method, ls) %>%
  gather(method, ls, -t, -Baseline) %>%
  mutate(ls = ls - Baseline) %>%
  ggplot() + geom_line(aes(t, ls, colour = method))


