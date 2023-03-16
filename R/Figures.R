# Figures for papers
library(tidyverse)
library(latex2exp)

#Design effect as a function of N, s, and theta and specificity
#This fixes the total number of units per site, but varies the pool size
de.data.spec <- expand.grid(`Pool Size` = c(1,2,4,8),
                            catch = c(16,24,32),
                            theta = seq(0.005,0.15,0.005),
                            specificity= c(0.98,1)) %>%
  mutate(N = catch/`Pool Size`,
         de = Vectorize(design_effect_cluster_fisher)(`Pool Size`,N,theta,1,specificity,0.1))
de.data.spec %>%
  mutate(across(c(`Pool Size`,catch,specificity), as.factor)) %>%
  mutate(specificity = forcats::fct_rev(specificity)) %>%
  ggplot(aes(x = theta, y = de, color = `Pool Size`)) +
  geom_line() +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0,NA)) +
  #scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cols = vars(catch), rows = vars(specificity), labeller = label_both) +
  xlab('Prevalence') +
  ylab('Design Effect')

## Similar to above but varying corr instead of specificity
de.data.corr <- expand.grid(`Pool Size` = c(1,2,4,8),
                            catch = c(16,32),
                            theta = seq(0.005,0.05,0.005),
                            correlation= c(0.01,0.1)) %>%
  mutate(N = catch/`Pool Size`,
         de = Vectorize(design_effect_cluster_fisher)(`Pool Size`,N,theta,1,1,correlation,'logitnorm'),
         cost = 1.2 * `Pool Size` * N + 4.33 * N + 40,
         var = de /(N * `Pool Size` * fi_pool_imperfect(1,theta,1,1)),
         unitcost = var * cost)

de.fig.corr <- de.data.corr %>%
  mutate(across(c(`Pool Size`,catch,correlation), as.factor)) %>%
  mutate(correlation = forcats::fct_rev(correlation)) %>%
  ggplot(aes(x = theta, y = de, color = `Pool Size`)) +
  geom_line() +
  geom_hline(yintercept = 1,linetype = 2) +
  coord_cartesian(ylim = c(0.8,NA)) +
  #scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cols = vars(catch), rows = vars(correlation), labeller = label_both) +
  xlab('Prevalence') +
  ylab('Design Effect') +
  theme(text = element_text(size = 12))
de.fig.corr

ggsave('./Figures/FigDesignEffectCorr.png',de.fig.corr, width = 6, height = 5, unit = 'in',dpi = 'retina')

unitcost.fig.corr <- de.data.corr %>%
  mutate(across(c(`Pool Size`,catch,correlation), as.factor)) %>%
  mutate(correlation = forcats::fct_rev(correlation)) %>%
  ggplot(aes(x = theta, y = unitcost, color = `Pool Size`)) +
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cols = vars(catch), rows = vars(correlation), labeller = label_both) +
  xlab('Prevalence') +
  ylab('Unit Cost of Information')
unitcost.fig.corr

optimalN.data.corr <-  de.data.corr %>%
  select(-c(N,de,`Pool Size`)) %>% unique %>%
  mutate(out = Vectorize(optimise_N_prevalence)(theta,
                                                cost.unit = 1.2,
                                                cost.test = 4.33,
                                                cost.location = 40,
                                                cost.collect = 0,
                                                catch.collect = catch,
                                                correlation = correlation, P = 1,
                                                form = 'logitnorm') %>% t %>% as.data.frame())
optimalN.data.corr %>%
  mutate(N = unlist(out$N),
         across(c(catch,correlation), as.factor),
         correlation = forcats::fct_rev(correlation)) %>%
  ggplot(aes(x = theta, y = N)) +
  geom_line() +
  geom_hline(yintercept = 1,linetype = 2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cols = vars(catch), rows = vars(correlation), labeller = label_both) +
  xlab('Prevalence') +
  ylab('Optimal number of pools')







#### Figures for technical paper

#plot ratio: how many time more total units do we need to test if are using group testing

plot_data <- expand.grid(prevalence = seq(0,0.1,by = 0.001),
                         sensitivity = c(0.9,0.99,0.999,1),
                         specificity = c(0.9,0.99,0.999,1),
                         s = c(50, 20,10, 5, 2, 1)) %>%
  mutate(ratio = fi_ratio_imperfect(s, prevalence, sensitivity, specificity, FALSE))

FigDesignEffect <- plot_data %>%
  subset(sensitivity == 1 ) %>%
  mutate(sensitivity = paste0("'Sensitivity: ' *varphi*' = ",sensitivity,"'"),
         specificity = paste0("'Specificity: ' *psi*' = "   ,specificity,"'")) %>%
  mutate(s = factor(s)) %>%
  ggplot(aes(x = prevalence, y = ratio, color = s)) +
  geom_line() +
  facet_grid(. ~ specificity ,
             labeller = label_parsed) +
  scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),breaks = seq(0,0.1,0.02)) +
  ylab(latex2exp::TeX(r"(Design effect: $D(\theta, N, s, 1, \psi)$)")) +
  xlab(latex2exp::TeX(r"(Prevalence: $\theta$)")) +
  guides(color = guide_legend(title = 'Group/pool size: s',nrow = 1)) +
  theme(legend.position = 'bottom')
FigDesignEffect
ggsave("./Figures/FigDesignEffect.png",FigDesignEffect)

FigDesignEffectSupp <- plot_data %>%
  mutate(sensitivity = paste0("'Sensitivity: ' *varphi*' = ",sensitivity,"'"),
         specificity = paste0("'Specificity: ' *psi*' = "   ,specificity,"'")) %>%
  mutate(s = factor(s),
         sensitivity = forcats::fct_rev(sensitivity)) %>%
  ggplot(aes(x = prevalence, y = ratio, color = s)) +
  geom_line() +
  facet_grid(sensitivity ~ specificity,scales = 'free_y',
             labeller = label_parsed) +
  scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),breaks = seq(0,0.1,0.02)) +
  ylab(latex2exp::TeX(r"(Design effect: $D(\theta, N, s, \varphi, \psi)$)")) +
  xlab(latex2exp::TeX(r"(Prevalence: $\theta$)")) +
  guides(color = guide_legend(title = 'Group/pool size: s', nrow = 1)) +
  theme(legend.position = 'bottom')
FigDesignEffectSupp
ggsave("./Figures/FigDesignEffectSupp.png",FigDesignEffectSupp)

#Optimal s for simple random sampling

optimal_s_data <- expand.grid(prevalence = seq(0.001,0.1,by = 0.001),
                              sensitivity = c(0.9,0.99,0.999,1),
                              specificity = c(0.9,0.99,0.999,1),
                              cost.unit = 1,
                              cost.test = c(0,0.1,1,10,Inf)) %>%
  mutate(out = Vectorize(optimise_s_prevalence)(prevalence, cost.unit, cost.test,
                                                sensitivity = sensitivity,specificity = specificity,
                                                max.s = 350, interval = 0) %>% t %>% as.data.frame()) %>%  mutate(s = unlist(out$s),
                                                                                                                  cost = unlist(out$cost))

fig_optimal_s <- optimal_s_data %>%
  #subset(sensitivity == 1 ) %>%
  mutate(specificity = paste0("'Specificity: ' *psi*' = "   ,specificity,"'")) %>%
  mutate(sensitivity = paste0("'Sensitivity: ' *varphi*' = ",sensitivity,"'")) %>%
  mutate(sensitivity = forcats::fct_rev(factor(sensitivity))) %>%
  ggplot(aes(x = prevalence, y = s, color = factor(cost.test))) +
  geom_line() + facet_grid(sensitivity ~ specificity,
                           #. ~ specificity,
                           labeller = label_parsed) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0,0.1,0.02)) +
  # scale_x_log10(labels = scales::percent_format(accuracy = 0.1),
  #               breaks = c(0.001,0.003,0.01,0.03,0.1)) +
  coord_cartesian(ylim = c(1, 250))+
  scale_y_log10() +
  ylab('Optimal pool/group size') +
  xlab(latex2exp::TeX(r"(Prevalence: $\theta$)")) +
  guides(color = guide_legend(title = latex2exp::TeX(r"(Cost ratio: $c_p/c_u$)"))) +
  theme(legend.position = 'bottom')
fig_optimal_s
ggsave("./Figures/FigOptimals.png",fig_optimal_s,height = 24, width= 24, unit = 'cm')


plot_data_interval <- expand.grid(prevalence = seq(0.001,0.1,by = 0.001),
                                  sensitivity = c(1),
                                  specificity = c(0.9,0.99,0.999,1),
                                  cost.unit = 1,
                                  cost.test = c(0.1,1,10),
                                  interval = c(0.1)) %>%
  mutate(out = Vectorize(optimise_s_prevalence)(prevalence, cost.unit, cost.test,
                                                sensitivity = sensitivity,specificity = specificity,
                                                max.s = 350, interval = interval) %>% t %>% as.data.frame()) %>%
  mutate(s = unlist(out$s),
         cost = unlist(out$cost),
         min_s = unlist(map(out$s_interval, ~.x[1])),
         max_s = unlist(map(out$s_interval, ~.x[2])))



FigOptimalsInterval <- plot_data_interval %>%
  subset(sensitivity == 1 &
           interval == 0.1) %>%
  mutate(cost.test = paste0("'Cost ratio: '*c[p] / c[u] * {phantom() == phantom()} * ",cost.test),
         specificity = paste0("'Specificity: ' *psi*' = "   ,specificity,"'")) %>%
  ggplot(aes(x = prevalence, y = s, ymin = min_s, ymax = max_s#,
             #color = cost.test, fill = cost.test
  )) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  facet_grid(#sensitivity ~ specificity,
    cost.test ~ specificity,
    labeller = label_parsed) +
  # scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),
  #                    breaks = seq(0,0.02,0.005)) +
  # scale_x_log10(labels = scales::percent_format(accuracy = 0.1),
  #                    breaks = c(0.001,0.003,0.01,0.03,0.1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0,0.1,0.02)) +
  scale_y_log10() +
  ylab('Pool/group size') +
  xlab(latex2exp::TeX(r"(Prevalence: $\theta$)")) +
  coord_cartesian(ylim = c(1, 320))
FigOptimalsInterval
ggsave("./Figures/FigOptimalsInterval.png",FigOptimalsInterval, height = 11.66, width= 24, unit = 'cm')

## Demonstrate that Fisher information for cluster sample sometimes but not always converges towards simple random for small rho

#Full fisher information matrix and top left element of inverse:
data_rho_limit <- expand.grid(rho = c(-1,10^(seq(-0.01,-2.2,-0.01))),
                              theta = c(0.01,0.1,0.3,0.4,0.5),
                              #theta = seq(0.01,0.5,0.01),
                              s = c(10),
                              N = c(5),
                              form = c('beta', 'discrete','logitnorm'),
                              stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(fi = ifelse(rho == -1,
                     list(N * fi_pool_imperfect(s,theta,1,1)),
                     list(fi_pool_imperfect_cluster(s,N,theta,1,1,rho,form = form))),
         var = solve(fi)[1,1],
         fi_theta = ifelse(rho == -1, fi,fi[1,1]),
         fi_rho = ifelse(rho == -1, NA,fi[2,2]),
         fi_cross = ifelse(rho == -1, NA,-fi[1,2])) %>%
  ungroup %>%
  pivot_longer(cols = var:fi_cross,names_to = 'measure')

data_rho_limit %>%
  group_by(form) %>%
  group_map(~{
    fisherplot <- .x %>%
      mutate(measure = factor(measure, levels = c('fi_theta','var','fi_cross','fi_rho')),
             measure = recode(measure,
                              fi_theta = paste0(as.character(latex2exp::TeX(r"($I(\theta,rho)_{1,1}$)"))),
                              var = as.character(latex2exp::TeX(r"($\Sigma_{1,1} \equiv I(\theta,\rho)^{-1}_{1,1}$)")),
                              fi_cross = as.character(latex2exp::TeX(r"($-I(\theta,rho)_{2,1} = -I(\theta,rho)_{1,2}$)")),
                              fi_rho = as.character(latex2exp::TeX(r"($I(\theta,rho)_{2,2}$)")))) %>%
      ggplot(aes(x = rho)) +
      geom_line(aes(y = value, color = factor(theta)),data = ~subset(.x,rho!=-1)) +
      geom_hline(aes(yintercept = value,
                     color = factor(theta)),
                 data = ~subset(.x, rho ==-1),
                 linetype = 2) +
      facet_wrap(facets = vars(measure),scale = 'free_y',labeller = label_parsed,strip.position = 'left') +
      scale_y_log10() +
      ylab('') +
      xlab(latex2exp::TeX(r"($\rho$)")) +
      guides(color = guide_legend(title = latex2exp::TeX(r"($\theta$)"))) +
      theme(strip.placement = "outside",
            strip.background = element_rect(fill='white'),
            text = element_text(size = 12))
      ggsave(paste0('FisherMatrix',.y,'.png'),fisherplot,width = 6, height = 5)
})

#Just the top left element of inverse Fisher information for a range of s and N

data_rho_limit_alt <- expand.grid(rho = c(-1,10^(seq(-0.01,-2.2,-0.01))),
                              theta = c(0.01,0.02,0.05,0.1),
                              #theta = seq(0.01,0.5,0.01),
                              s = c(50,100,200,400,1000,2000),
                              N = 2,
                              form = c('beta'),
                              stringsAsFactors = FALSE) %>%
  subset(!(theta > 0.01 & s >= 1000)) %>%
  rowwise() %>%
  mutate(fi = ifelse(rho == -1,
                     list(N * fi_pool_imperfect(s,theta,1,1)),
                     tryCatch(list(fi_pool_imperfect_cluster(s,N,theta,1,1,rho,form = form)), error = function(e){NA})),
         var = tryCatch(solve(fi)[1,1], error = function(e){NA})) %>%
  ungroup

data_rho_limit_alt %>%
  group_by(form) %>%
  group_map(~{
    varplot <- .x %>%
      mutate(s = fct_rev(factor(s)),
             theta = paste(as.character(latex2exp::TeX(r"($\theta=$)")),theta, sep = ' * ')) %>%
      ggplot(aes(x = rho)) +
      geom_line(aes(y = var, color = s),data = ~subset(.x,rho!=-1)) +
      geom_hline(aes(yintercept = var,
                     color = s),
                 data = ~subset(.x, rho ==-1),
                 linetype = 2) +
      facet_wrap(facets = vars(theta),labeller = label_parsed, scales = 'free_y') +
      scale_y_log10() +
      xlab(latex2exp::TeX(r"($\rho$)")) +
      guides(color = guide_legend(title = latex2exp::TeX(r"($s$)"))) +
      theme(strip.placement = "outside",
            strip.background = element_rect(fill='white'),
            text = element_text(size = 12)) +
      ylab(latex2exp::TeX(r"($\Sigma_{1,1} \equiv I(\theta,\rho)^{-1}_{1,1}$)"))
    ggsave(paste0('Variance',.y,'.png'),varplot,width = 6, height = 5)
  })



rho <- 10^(seq(-0.01,-2.8,-0.1))
theta <- 0.5
s <- 10
N <- 5
fi <- c(map(rho,~fi_pool_imperfect_cluster(s,N,theta,1,1,.x, form = 'beta')),
        (N*fi_pool_imperfect(s,theta,1,1)))
fi
plot(c(rho,0),unlist(map(fi,~solve(.x)[1,1])),log = 'y')
