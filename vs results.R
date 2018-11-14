library(tidyverse)

# GK Model

gk <- tibble()
for(i in 1:200){
  file <- paste0('ma1gk/rep', i, '.csv')
  if(file.exists(file)){
    tmp <- read_csv(file, col_types = cols())
    tmp$id <- i
    gk <- rbind(gk, tmp)
  }
}

length(unique(gk$id))

gk %>%
  select(-unique) %>%
  spread(method, ls) %>%
  gather(method, ls, `S-ABC`, `VS-ABC`, `VS-ABC-2`) %>%
  mutate(ls = ls - ABC) %>%
  group_by(t, method) %>%
  summarise(ls = 10*mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, ls, colour = method))

gk %>% 
  group_by(t, method) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, ls, colour = method)) + 
  theme_bw() + theme(legend.position = 'bottom') + 
  labs(colour = 'Method', x = 'T', y = 'Mean Forecast Logscore') 

gk %>%
  filter(method %in% c('S-ABC', 'VS-ABC-2')) %>%
  group_by(t, method) %>%
  summarise(min = min(unique), 
            max = max(unique),
            mean = mean(unique)) %>%
  ggplot() + geom_line(aes(t, mean, colour = method)) +
  geom_ribbon(aes(t, ymin = min, ymax = max, fill = method), alpha = 0.3) + 
  theme_bw() + theme(legend.position = 'bottom') + 
  labs(colour = 'Method', x = 'T', y = 'Number of Distinct Particles') + guides(fill=FALSE)





# SV with Normal Approx
svN <- tibble()
for(i in 1:200){
  file <- paste0('SV/rep', i, '.csv')
  if(file.exists(file)){
    tmp <- read_csv(file, col_types = cols())
    tmp$id <- i
    svN <- rbind(svN, tmp)
  }
}

# SV Beta
svB <- tibble()
for(i in 1:200){
  file <- paste0('SV_Beta/rep', i, '.csv')
  if(file.exists(file)){
    tmp <- read_csv(file, col_types = cols())
    tmp$id <- i
    svB <- rbind(svB, tmp)
  }
}

svN <- rbind(svN,
             svB %>% filter(method == 'ABC'))


svN$approx <- 'Normal'
svB$approx <- 'Beta'

sv <- rbind(svN, svB)

sv %>%
  select(-unique) %>%
  mutate(method = paste(method, approx)) %>%
  select(-approx) %>%
  #filter(!method %in% c('ABC Normal', 'ABC Beta', 'S-ABC Normal')) %>%
 #select(-unique, -approx) %>%
  spread(method, ls) %>%
  gather(method, ls, -t, -h, -id, -`ABC Beta`) %>%
  mutate(ls = ls - `ABC Beta`) %>%
  group_by(method, t) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, ls, colour = method)) #+ facet_wrap(~approx)  + 
  theme_bw() + theme(legend.position = 'bottom') + 
  labs(colour = 'Method', x = 'T', y = 'Mean Forecast Logscore') 

sv %>%
  select(-unique) %>% 
  spread(method, ls) %>%
  gather(method, ls, `S-ABC`, `VS-ABC`, `VS-ABC-2`) %>%
  mutate(ls = ls - ABC) %>%
  group_by(method, approx, t) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  spread(approx, ls) %>%
  mutate(ls = Normal - Beta) %>%
  ggplot() + geom_line(aes(t, ls, colour = method)) #+ facet_wrap(~approx, ncol = 1)

# vmcmc

vmcmc <- tibble()
for(i in 1:200){
  file <- paste0('VMCMC/rep', i, '.csv')
  if(file.exists(file)){
    tmp <- read_csv(file, col_types = cols())
    tmp$id <- i
    vmcmc <- rbind(vmcmc, tmp)
  }
}

vmcmc %>%
  group_by(method) %>%
  summarise(ls = mean(ls, na.rm = TRUE))

vmcmc %>%
  group_by(method, t, h) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, ls, colour = method)) + facet_wrap(~h)

vmcmc %>%
  filter(h == 1) %>%
  group_by(id, method) %>%
  mutate(time = ifelse(method == 'MCMC', time, cumsum(time))) %>%
  ungroup() %>%
  group_by(method, t) %>%
  summarise(meantime = mean(time), 
            upper = quantile(time, 0.975),
            lower = quantile(time, 0.025)) %>%
  ggplot() + geom_line(aes(t, meantime / 60, colour = method)) + 
  geom_ribbon(aes(t, ymin = lower / 60, ymax = upper / 60, fill = method), alpha = 0.2) + 
  labs(x = 'T', y = 'Mean (Cumulative) Run Time, Minutes', colour = 'Method') + 
  theme_bw() + theme(legend.position = 'bottom') +  guides(fill=FALSE) + 
  scale_x_continuous(breaks = seq(100, 1000, 100)) + scale_y_continuous(breaks = seq(5, 30, 5))

vmcmc %>%
  select(-time) %>%
  spread(method, ls) %>%
  mutate(ls = MCMC - `V-MCMC`) %>%
  group_by(t) %>%
  summarise(meanls = mean(ls, na.rm = TRUE), 
            upper = quantile(ls, 0.975,  na.rm = TRUE),
            lower = quantile(ls, 0.025,  na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, meanls)) + 
  geom_ribbon(aes(t, ymin = lower, ymax = upper), alpha = 0.2)+ 
  labs(x = 'T', y = 'Mean Difference in Logscore (MCMC - V-MCMC)', colour = 'Method') + 
  theme_bw() + theme(legend.position = 'bottom') + 
  scale_x_continuous(breaks = seq(100, 1000, 100))

  