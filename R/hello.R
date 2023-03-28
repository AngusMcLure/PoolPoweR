




specificity <- 0.99
for(N in c(100, 90, 50, 20,10, 5, 2, 1)){
  plot.function(function(p)ifelse(fi_ratio_imperfect(N,p,sensitivity,specificity,F) <N,
                                  fi_ratio_imperfect(N,p,sensitivity,specificity,F),
                                  N),
                to = 0.1,add = T, log = 'y')
}

for(N in c(2:10)){
  plot.function(function(p)fi_ratio_imperfect(N,p,sensitivity,specificity,F),to = 0.1,add = T)
}

for(p in 10^((-3:0)/3)){
  plot.function(function(s)fi_ratio_imperfect(s,p,sensitivity,specificity),to = 100,add = T,log = 'y')
}


### Testing out Fisher information

N <- 4; corr <- 0.4; s = 4; theta <- 0.01; sens <- 1; spec <- 1;
rho. <- corr^-1 - 1
Alpha <- theta * rho.; Alpha
Beta <- (1 - theta) * rho.; Beta


fi_pool_imperfect(s,theta,sens,spec) * N
fi_pool_imperfect_cluster(s,N,theta,sens,spec,corr, form = 'beta')
fi_pool_imperfect_cluster(s,1,theta,sens,spec,corr,form = 'beta') * N

solve(fi_pool_imperfect_cluster(    2,        6,    theta,sens,spec,corr,form = 'logitnorm'))
solve(fi_pool_imperfect_cluster(    2,        6,    theta,sens,spec,corr,form = 'cloglognorm'))



solve(fi_pool_imperfect_cluster(    2,        6,    theta,sens,spec,corr, form = 'beta'))
solve(fi_pool_imperfect_cluster(c(1,2,3), c(1,4,1), theta,sens,spec,corr, form = 'beta'))
solve(fi_pool_imperfect_cluster(c(1,2,3), c(2,2,2), theta,sens,spec,corr, form = 'beta'))
solve(fi_pool_imperfect_cluster(c(1,  3), c(3,  3), theta,sens,spec,corr, form = 'beta'))

design_effect_cluster_fisher(s,N,theta,sens,spec,corr, form = 'beta')
design_effect_cluster_fisher(s,N,theta,sens,spec,corr, form = 'logitnorm')
design_effect_cluster_fisher(s,N,theta,sens,spec,corr, form = 'cloglognorm')
design_effect_cluster_fisher(s,N,theta,sens,spec,corr, form = 'discrete')
(1 + (N*s - 1) * corr)
(1 + (N*s - 1) * corr) * fi_ratio(s,theta)


plot_data <- expand.grid(#prevalence = seq(0.01,0.11,by = 0.02),
  prevalence = 0.02,
  sensitivity = c(0.9,0.99,1),
  specificity = c(0.9,0.99,1),
  correlation = seq(0.1,0.3,by = 0.1),
  s = 2^(0:4),
  N = 2^(0:4),
  form = c('beta'),
  reference = c(T,F)) %>%
  mutate(de = Vectorize(design_effect_fisher)(s,prevalence,sensitivity,specificity, correlation, N, form, reference))

FigDE <- plot_data %>%
  subset(correlation == 0.2 & sensitivity == 1) %>%
  #mutate(sensitivity = forcats::fct_rev(factor(sensitivity))) %>%
  ggplot(aes(y = s, x = de, color = factor(N))) +
  geom_line() + facet_grid(#sensitivity ~ specificity,
    specificity ~ reference,
    labeller = label_both) +
  scale_y_log10()
#scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +

FigDE




#############################################
# Optimising design for estimating prevalence
#############################################

cu <- 1.4; ct <- 4; cl <- 40;
cu <- 1.4; ct <- 4; cl <- 400;

P <- 2; catch.mean <- 20; cc <- catch.mean * cu * 0.5 / P
theta <- 0.01; rho <- 0.1
optimise_s_prevalence (theta,cu,ct,cl,rho, N = 4,max.s = 250,interval = 0.2)
optimise_N_prevalence(theta,cu-cc/catch.mean,ct,cl,cc,catch.mean,rho,P,max.s = 50)
optimise_sN_prevalence(theta,cu,ct,cl,rho,       max.s = 50)
optimise_NP_prevalence(theta,cu-cc/catch.mean,ct,cl,cc,catch.mean, rho,max.s = 50)


opt <- optimise_NP_prevalence(theta,cu-cc/catch.mean,ct,cl,cc,catch.mean, rho,max.s = 50)
sqrt(solve(fi_pool_imperfect_cluster(opt$s,opt$N,theta,1,1,0.1))[1,1]) * 1.96
sqrt(1/(opt$N * fi_pool_imperfect(opt$s,theta,1,1))) * 1.96
sqrt(1/(opt$catch * fi_pool_imperfect(1,theta,1,1))) * 1.96

design_effect_cluster_fisher(opt$s,opt$N,theta,1,1,rho)

unit_fi_cost_clustered(4,4,theta,rho,1,1,cu,ct,cl)

Vectorize(optimise_s_prevalence)(theta,cu,ct,cl,rho, N = 2:10,max.s = 50)


Vectorize(optimise_s_prevalence)(prevalence = c(0.0001,0.001,0.002,0.005,0.01),
                                 sensitivity = 1, specificity = 1,
                                 cost.unit = 3.3 + 0.68 + 0.45,
                                 cost.pool = 4.66 + 2.73, max.s = 200)

Vectorize(optimise_s_prevalence)(prevalence = c(0.0001,0.001,0.002,0.005,0.01),
                                 sensitivity = 1, specificity = 1,
                                 cost.unit = 3.3 + 0.68 + 0.45,
                                 cost.pool = 4.66 + 2.73,
                                 cost.location = 8.93, correlation = 0.01, N = 5,
                                 max.s = 200)

Vectorize(optimise_sN_prevalence)(prevalence = c(0.0001,0.001,0.002,0.005,0.01),
                                  sensitivity = 1, specificity = 1,
                                  cost.unit = 3.3 + 0.68 + 0.45,
                                  cost.pool = 4.66 + 2.73,
                                  cost.location = 8.93, correlation = 0.01,
                                  max.s = 200)

Vectorize(optimise_sN_prevalence)(prevalence = 0.01,
                                  sensitivity = 1, specificity = 1,
                                  cost.unit = 3.3 + 0.68 + 0.45,
                                  cost.pool = 4.66 + 2.73,
                                  cost.location = 8.93,
                                  correlation = c(0,0.01,0.1,0.3),
                                  max.s = 20)


Vectorize(unit_fi_cost_clustered)(s = 1:20, N = 1,prevalence = 0.01, correlation = 0.001,
                                  sensitivity = 1, specificity = 1,
                                  cost.unit =  3.3 + 0.68 + 0.45,
                                  cost.pool = 4.66 + 2.73,
                                  cost.location = 8.93, form = 'beta')

N <- 10; Vectorize(optimise_s_prevalence)(0.01, 0.9, 0.9, 1, 2, 4, 0.1,N,max.s = 30)

theta <- 0.01; 1-(1-theta)^optimise_s_prevalence(theta, 0.9, 0.99, 1, 0, max.s = 30)








plot_data_alt <- expand.grid(prevalence = c(0.001,0.003,0.01,0.03),
                             sensitivity = c(0.9,0.99,0.999,1),
                             specificity = seq(0.9,1,by = 0.001),
                             cost.unit = 1,
                             cost.pool = c(0,0.1,1,10,Inf)) %>%
  mutate(s = Vectorize(optimise_s_prevalence)(prevalence,sensitivity, specificity, cost.unit,cost.pool, max.s = 1000))

FigOptimals_alt <- plot_data_alt %>%
  subset(sensitivity == 1 ) %>%
  mutate(prevalence = factor(paste0(prevalence*100,'%'))) %>%
  ggplot(aes(x = specificity, y = s, color = factor(cost.pool))) +
  geom_line() + facet_grid(
    . ~ prevalence,
    labeller = label_both) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, 50))
FigOptimals_alt
ggsave("FigOptimals_alt.png",FigOptimals_alt)


plot_data_PPV <- expand.grid(prevalence = seq(0.001,0.1,by = 0.001),
                             specificity = seq(0.9,1,by = 0.001),
                             sensitivity = 1,
                             cost.unit = 1,
                             cost.pool = c(0,0.1,1,10,Inf)) %>%
  mutate(PPV = (1 + (1 - specificity)/(prevalence * sensitivity))^-1,
         s = Vectorize(optimise_s_prevalence)(prevalence,sensitivity, specificity, cost.unit,cost.pool, max.s = 300))

plot_data_PPV %>%
  subset(sensitivity == 1) %>%
  ggplot(aes(x= PPV, y = s,color = prevalence)) +
  geom_point() +
  facet_grid(.~cost.pool)


