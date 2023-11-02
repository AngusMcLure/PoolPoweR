#plot ratio: how many time more total units do we need to test if are using group testing

plot_data <- expand.grid(prevalence = seq(0,0.1,by = 0.001),
                         sensitivity = c(0.9,0.99,0.999,1),
                         specificity = c(0.9,0.99,0.999,1),
                         s = c(50, 20,10, 5, 2, 1)) %>%
  mutate(ratio = fi_ratio(s, prevalence, sensitivity, specificity, FALSE))

FigDesignEffectSupp <- plot_data %>%
  mutate(s = factor(s),
         sensitivity = fct_rev(factor(sensitivity))) %>%
  ggplot(aes(x = prevalence, y = ratio, color = s)) +
  geom_line() +
  facet_grid(sensitivity ~ specificity, scales = 'free_y',
             labeller = label_bquote(varphi == .(as.character(sensitivity)),
                                     psi == .(specificity))) +
  scale_y_log10() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),breaks = seq(0,0.1,0.02)) +
  ylab(latex2exp::TeX(r"(Design effect: $D(\theta, N, s, \varphi, \psi)$)")) +
  xlab(bquote(theta)) +
  guides(color = guide_legend(title = 's', nrow = 1)) +
  theme(legend.position = 'bottom',
        text = element_text(size = 12))
FigDesignEffectSupp
ggsave("./Figures/FigDesignEffectSupp.png",FigDesignEffectSupp,
       height = 7, width= 7, unit = 'in')

FigDesignEffectSimple <- FigDesignEffectSupp
FigDesignEffectSimple$data <- FigDesignEffectSimple$data %>%
  subset(sensitivity == 1 & specificity == 1 & s !=50)

(FigDesignEffectSimple +
  labs(x = 'Prevalence',
       y = 'Design Effect') +
    facet_null() +
    scale_y_continuous()+
    guides(color = guide_legend(title = 'Pool Size', nrow = 1))) %>%
  ggsave("./Figures/FigDesignEffectSimple.png",.,
         height = 4, width= 4, unit = 'in')


#Optimal s for simple random sampling
optimal_s_data <- expand.grid(prevalence = seq(0.001,0.1,by = 0.001),
                              sensitivity = c(0.9,0.99,0.999,1),
                              specificity = c(0.9,0.99,0.999,1),
                              cost.unit = 1,
                              cost.pool = c(0,0.1,1,10,Inf)) %>%
  rowwise() %>%
  mutate(out = optimise_s_prevalence(prevalence, cost.unit, cost.pool,
                                     sensitivity = sensitivity,specificity = specificity,
                                     max.s = 350, interval = 0) %>% t %>% as.data.frame()) %>%
  mutate(s = unlist(out$s),
         cost = unlist(out$cost))

fig_optimal_s <- optimal_s_data %>%
  mutate(sensitivity = forcats::fct_rev(factor(sensitivity))) %>%
  ggplot(aes(x = prevalence, y = s, color = factor(cost.pool))) +
  geom_line() +
  facet_grid(sensitivity ~ specificity,
             labeller = label_bquote(varphi == .(as.character(sensitivity)),
                                     psi == .(specificity))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0,0.1,0.02)) +
  coord_cartesian(ylim = c(1, 250))+
  scale_y_log10() +
  ylab('Optimal s') +
  xlab(bquote(theta)) +
  guides(color = guide_legend(title = latex2exp::TeX(r"(Cost ratio: $c_p/c_u$)"))) +
  theme(legend.position = 'bottom', text = element_text(size = 12))
fig_optimal_s
ggsave("./Figures/FigOptimals.png",fig_optimal_s,
       height = 7, width= 7, unit = 'in')

fig_optimal_s_simple <- fig_optimal_s
fig_optimal_s_simple$data <- fig_optimal_s_simple$data %>% subset(sensitivity == 1 & specificity == 1)
fig_optimal_s_simple + facet_null()
ggsave("./Figures/FigOptimalsSimple.png",fig_optimal_s_simple,
       height = 7, width= 7, unit = 'in')

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
      #scale_y_log10() +
      ylab('') +
      xlab(latex2exp::TeX(r"($\rho$)")) +
      guides(color = guide_legend(title = latex2exp::TeX(r"($\theta$)"))) +
      theme(strip.placement = "outside",
            strip.background = element_rect(fill='white'),
            text = element_text(size = 12))
    ggsave(paste0('./Figures/FisherMatrix',.y,'.png'),fisherplot,
           width = 6, height = 5)
  })

#Just the top left element of inverse Fisher information for a range of s and N

data_rho_limit_alt <- expand.grid(rho = c(-1,10^(seq(-0.01,-2.2,-0.01))),
                                  theta = c(0.01,0.02,0.05,0.1),
                                  #theta = seq(0.01,0.5,0.01),
                                  s = c(50,100,200,400,1000,2000),
                                  N = 2,
                                  form = c('beta','discrete'),
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
    ggsave(paste0('./Figures/Variance',.y,'.png'),varplot,width = 6, height = 5)
  })


# Unit information plots for fixed cost structure, range of prevalence, correlation (amount and form)

plot_data_clustered_alt <- expand.grid(prevalence = c(0.01, 0.05, 0.1),
                                       sensitivity = 1,specificity = 1,
                                       cost.unit = 1,cost.pool = 4, cost.location = 40,
                                       correlation = c(NA,0,0.01,0.1),
                                       N = c(2,4,8), s = 1:50,
                                       form = c('discrete', 'beta', 'logitnorm'),
                                       stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(unitcost = ifelse(is.na(correlation),
                           unit_fi_cost(s, prevalence, sensitivity,specificity,cost.unit, cost.pool),
                           unit_fi_cost_clustered(s = s, N = N, prevalence = prevalence,
                                                  correlation = correlation,
                                                  sensitivity = sensitivity, specificity = specificity,
                                                  cost.unit =  cost.unit,
                                                  cost.pool = cost.pool,
                                                  cost.location = cost.location, form = form)))

plot_data_clustered_alt %>%
  ungroup %>% group_by(form) %>%
  group_map(~{
    .x <- .x %>%
      mutate(correlation = ifelse(is.na(correlation),
                                  'NA (simple random survey)',
                                  ifelse(correlation == 0,
                                         "0 (known)" ,
                                         as.character(correlation))))
    unit_info_cluster <- .x %>%
      mutate(N = as.factor(N)) %>%
      ggplot(aes(x = s, y = unitcost, color = N)) +
      geom_line() +
      facet_grid(prevalence~correlation,
                 labeller = label_bquote(theta == .(prevalence), rho == .(correlation)),
                 scales = 'free') +
      scale_y_log10() +
      xlab('s') + ylab('unit cost of information') +
      theme(text = element_text(size = 12))
    unit_info_cluster %>%
      ggsave(paste('./Figures/unit info cluster',.y,'.png')
             ,., width = 6.5, height = 6.5, unit = 'in')
    unit_info_cluster_alt <- .x %>%
      mutate(correlation = factor(correlation,
                                  levels = c('0.1','0.01',"0 (known)",
                                             'NA (simple random survey)'))) %>%
      ggplot(aes(x = s, y = unitcost, color = correlation)) +
      geom_line() +
      facet_grid(prevalence~N,
                 labeller = label_bquote(theta == .(prevalence * 100) * '%', N == .(N)),
                 scales = 'free') +
      scale_y_log10() +
      labs(x = 's', y = 'Unit cost of information',
           color = latex2exp::TeX(r"($\rho$)")) +
      theme(text = element_text(size = 12),
            legend.position = 'bottom')
    unit_info_cluster_alt %>%
      ggsave(paste0('./Figures/unit info cluster alt ',.y,'.png')
             ,., width = 6.5, height = 6.5, unit = 'in')
  })

## Optimal s for given number of pools per cluster

plot_data_clustered <- expand.grid(prevalence = seq(0.005, 0.1, 0.001),
                                   sensitivity = 1,specificity = 1,
                                   cost.unit = 1,cost.pool = 4, cost.location = 40,
                                   correlation = c(NA,0,0.01,0.1),
                                   N = c(2,4,8),
                                   form = c('beta'),
                                   stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(out = optimise_s_prevalence(prevalence, cost.unit, cost.pool, cost.location, correlation, N = N, form = 'beta', max.s = 200, interval = 0.1) %>% t %>% as.data.frame()) %>%
  mutate(s = unlist(out$s),
         min_s = unlist(map(out$s_interval, ~.x[1])),
         max_s = unlist(map(out$s_interval, ~.x[2])))


optimal_cluster <- plot_data_clustered %>%
  mutate(across(c(min_s,max_s),~{ifelse(correlation %in% c(0,0.1),.x,NA)}),
         correlation = ifelse(is.na(correlation),
                              'NA (simple random survey)',
                              ifelse(correlation == 0,
                                     "0 (known)" ,
                                     as.character(correlation))),
         correlation = factor(correlation, levels = c('0.1','0.01',"0 (known)",'NA (simple random survey)'))) %>%
  ggplot(aes(x = prevalence, y = s, ymax = max_s,ymin = min_s,
             fill = correlation, color = correlation)) +
  geom_line() +
  geom_ribbon(alpha = 0.3,show.legend = FALSE, linetype = 2) +
  facet_wrap(N~., labeller = label_bquote(N == .(N))) +
  labs(y = 'Optimal s',
       x = bquote(theta),
       color = bquote(rho)) +
  theme(legend.position = 'bottom', text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, NA)) +
  scale_x_continuous(labels = scales::percent_format(1))
optimal_cluster
ggsave('./Figures/optimal cluster.png',optimal_cluster, width = 6.5, height = 4.5, unit = 'in')

fig_opmital_s_simple <-
  plot_data_clustered %>% subset(is.na(correlation)) %>%
  ggplot(aes(x = prevalence, y = s, ymin = min_s, ymax = max_s))+
  geom_line()+
  geom_ribbon(alpha = 0.3,show.legend = FALSE, linetype = 2)+
  labs(x = 'Prevalence', y = 'Pool Size') +
  theme(text = element_text(size = 12)) +
  scale_x_continuous(labels = scales::percent_format(1))


ggsave('./Figures/optimal s simple.png',fig_opmital_s_simple, width = 4, height = 4, unit = 'in')



## Optimal s, N and catch (sN)

plot_data_clustered_sN <- expand.grid(prevalence = seq(0.005,0.1,by = 0.001),
                                      cost.unit = 1,
                                      cost.pool = c(4,40),
                                      cost.location = c(4,40),
                                      correlation = c(NA,0,0.01,0.1,0.3),
                                      N = c(2,4,8)) %>%
  rowwise() %>%
  mutate(out = optimise_sN_prevalence(prevalence,
                                      cost.unit, cost.pool, cost.location,
                                      correlation, form = 'beta',
                                      max.s = 200, max.N = 20) %>% t %>% as.data.frame) %>%
  mutate(s = unlist(out$s),
         N = unlist(out$N),
         catch = unlist(out$catch))


optimal_cluster_sN <- plot_data_clustered_sN %>%
  mutate(simplerandom = is.na(correlation)) %>%
  mutate(correlation = ifelse(is.na(correlation),
                              'NA (simple random survey)',
                              ifelse(correlation == 0,
                                     "0 (known)" ,
                                     as.character(correlation)))) %>%
  mutate(correlation = factor(correlation, levels = c('0.3','0.1','0.01',"0 (known)",'NA (simple random survey)'))) %>%
  pivot_longer(c(s, N,catch),names_to = 'outcome') %>%
  mutate(outcome = paste('Optimal', recode(outcome,
                                           catch = 'sN')),
         outcome = fct_relevel(outcome,'Optimal s', 'Optimal N')) %>%
  subset(value != Inf & !is.na(correlation)) %>%
  ggplot(aes(x = prevalence, y = value,
             color = correlation,
             linetype = simplerandom)) +
  geom_line() +
  facet_grid(outcome ~ cost.pool + cost.location, scales = 'free',
             labeller = label_bquote(rows = .(outcome),
                                     cols = c[u] * {phantom() == phantom()} * 1 * ';' *
                                       c[p] * {phantom() == phantom()} * .(cost.pool) * ';' *
                                       c[s] * {phantom() == phantom()} * .(cost.location)),
             switch = 'y'
  ) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0,NA)) +
  scale_x_continuous(breaks = seq(0,0.1, 0.02), labels = scales::percent_format(1)) +
  labs(x = bquote(theta),
       y = element_blank(),
       color = bquote(rho)) +
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        strip.placement = 'outside',
        strip.background = element_rect(fill = 'white')) +
  guides(linetype = 'none')
optimal_cluster_sN
ggsave('./Figures/optimal cluster sN.png',optimal_cluster_sN, width = 6.5, height = 4.5, unit = 'in')

optimal_cluster_sN_simple <- optimal_cluster_sN

optimal_cluster_sN_simple$data <- optimal_cluster_sN$data  %>%
  subset(cost.pool == 4 &
           correlation !='NA (simple random survey)'&
           correlation != '0.01') %>%
  mutate(outcome = recode(outcome,
                          `Optimal s` = 'Optimal\nPool Size',
                          `Optimal N` = 'Optimal\nPool Number',
                          `Optimal sN` = 'Optimal\nCatch'))

optimal_cluster_sN_simple <- optimal_cluster_sN_simple +
  labs(x = 'Prevalence',
       color = 'ICC')+
  facet_grid(outcome ~ cost.pool + cost.location, scales = 'free',
             labeller = label_bquote(rows = .(outcome),
                                     cols = c[u] * {phantom() == phantom()} * 1 * ';' *
                                       c[t] * {phantom() == phantom()} * .(cost.pool) * ';' *
                                       c[s] * {phantom() == phantom()} * .(cost.location)),
             switch = 'y'
  )
optimal_cluster_sN_simple
ggsave('./Figures/optimal cluster sN simple.png',
       optimal_cluster_sN_simple, width = 6.5, height = 4.5, unit = 'in')

### Design effect

de.data.corr <- expand.grid(s = c(1,4,8),
                            catch = c(16,32,64),
                            theta = seq(0.005,0.1,0.0025),
                            correlation= c(0,0.01,0.1),
                            form = c('beta','discrete'),
                            stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(N = catch/s,
         de = design_effect_cluster_fisher(s,N,theta,1,1,correlation,form),
         cost = 1 * s * N + 4 * N + 40,
         var = de /(N * s * fi_pool_imperfect(1,theta,1,1)),
         unitcost = var * cost)

de.data.corr %>%
  group_by(form) %>%
  group_map(~{
    de.fig.corr <- .x %>%
      mutate(s =  as.factor(s),
             correlation = ifelse(correlation == 0,
                                  "0 (known)" ,
                                  as.character(correlation))) %>%
      ggplot(aes(x = theta, y = de, color = s)) +
      geom_line() +
      geom_hline(yintercept = 1,linetype = 2) +
      coord_cartesian(ylim = c(0,NA),xlim = c(0,NA)) +
      #scale_y_log10() +
      scale_x_continuous(breaks = seq(0,0.1, 0.02), labels = scales::percent_format(1)) +
      facet_grid(cols = vars(catch), rows = vars(correlation),
                 labeller = label_bquote(rows = rho == .(correlation), cols = N %*% s == .(catch)), scales = 'free') +
      labs(x = bquote(theta),
           y = 'Design Effect') +
      theme(text = element_text(size = 12),
            legend.position = 'bottom')

    de.fig.corr %>%
      ggsave(paste0('./Figures/FigDesignEffectCorr',.y,'.png'),.,
             width = 6.5, height = 5, unit = 'in',dpi = 'retina')

    de.fig.corr$data <- de.fig.corr$data %>% subset(catch %in% c(16,64) & correlation != '0.01')
    (de.fig.corr +
        coord_cartesian(ylim = c(1,NA)) +
        scale_y_log10() +
        labs(x = 'Prevalence',
             color = 'Pool Size') +
        facet_grid(cols = vars(catch), rows = vars(correlation),
                   labeller = label_bquote(rows = ICC == .(correlation), cols = Catch == .(catch)))) %>%
      ggsave(paste0('./Figures/FigDesignEffectCorr',.y,'Simple.png'),.,
             width = 4, height = 4, unit = 'in',dpi = 'retina')

  })



### Sample size calculations

sample.size.data <- de.data.corr %>%
  mutate(theta.null = 0.01, power = 0.8, sig.level = 0.05) %>%
  rowwise() %>%
  mutate(Arcsine = ifelse(abs(theta - theta.null) < .Machine$double.eps,
                                 NA, pwr::pwr.p.test(h = pwr::ES.h(theta,theta.null),
                                                     power = power, sig.level = sig.level,
                                                     alternative = ifelse(theta < theta.null, 'less', 'greater'))$n),
         `Fisher matrix` = ifelse(abs(theta - theta.null) < .Machine$double.eps,
                                 NA, sample_size_prevalence(s,N,theta.null,theta,
                                                            power,sig.level = sig.level,
                                                            alternative = ifelse(theta < theta.null, 'less', 'greater'),
                                                            correlation = correlation, form = form)),
         `Design effect` = Arcsine * de)

sample.size.data %>%
  pivot_longer(cols = c(Arcsine,`Fisher matrix`,`Design effect`), names_to = 'method',values_to = "n") %>%
  mutate(s =  as.factor(s),
         correlation = ifelse(correlation == 0,
                              "0 (known)" ,
                              as.character(correlation)),
         method = factor(method, levels=c('Design effect','Fisher matrix','Arcsine'))) %>%
  group_by(form) %>%
  group_map(~{
    plot <- .x %>%
      ggplot(aes(x = theta, y = n, color = method, linetype = s)) +
      geom_line() +
      facet_grid(vars(correlation),vars(catch), labeller = label_bquote(rho == .(correlation),N %*% s == .(catch))) +
      scale_y_log10() +
      coord_cartesian(xlim = c(0,NA)) +
      scale_x_continuous(labels = scales::percent_format(1)) +
      labs(x = bquote(theta),
           y = bquote("minimum sample size " * (n == J%*%N%*%s)))
    plot %>%
      ggsave(paste0('./Figures/SampleSize',.y,'.png'),.,
             width = 6.5, height = 5, unit = 'in',dpi = 'retina')
  })


