library(tidyverse)

#Install then load specific GitHub commit. Use commented out version if loading local version
devtools::install_github('AngusMcLure/PoolPoweR@74c0031a1d579e2d62ef357bc6909f85e17d474d')
library(PoolPoweR)
#devtools::load_all("C:/Users/u4859599/Documents/GitHub/PoolPoweR/")



# functions for optimising pool size for identification (Dorfman)

#expected number of tests per unit screened
exp_test_dorf <- function(pool_size, prevalence){
  if(pool_size == 1){
    return(1)
  }else if(pool_size >1){
    (1 - (1-prevalence)^pool_size) + 1/pool_size
  }else{
    stop('pool size must be greater than or equal to 1')
  }
}

#optimal pool_size (minimise expected number of tests per unit screened)

opts_dorf <- function(prevalence){
  
  if(exp_test_dorf(2, prevalence)> 1)return(1) #if prevalence is so high that even a pool of two is less efficient that individual testing, optimal is individual testing!
  
  opts <- optimise(\(s)exp_test_dorf(s, prevalence),
                        c(1, 1/sqrt(prevalence) * 5))$minimum
  exp_test_dorf_upper <- exp_test_dorf(ceiling(opts), prevalence)
  exp_test_dorf_lower <- exp_test_dorf(floor(opts), prevalence)
  
  ifelse(exp_test_dorf_upper < exp_test_dorf_lower, ceiling(opts),floor(opts))
}




#optimal pool sizes (and pool-level positivity) for different prevalence and cost regimes

scenarios <- expand.grid(p = c(10^seq(-2.7, -0.5, 0.01), seq(0.002, 0.4, 0.001)), #set of prevalences that give you a good spread on both log and linear scales
                         c_u = c(0, 0.1, 1, 10)) #ratio of costs unit:test

d_prev <- scenarios %>%
  rowwise() %>%
  mutate(regime = 'Prevalence - no retesting',
         opts = optimise_prevalence(fixed_design(),prevalence = p,
                                    cost_unit = c_u,cost_pool = 1,
                                    cost_cluster = NA,correlation = NA,
                                    max_s = 1000)$sample_design$pool_size, #optimal pool size for estimating prevalence
         optpp = 1 - (1-p)^opts, #optimal pool-level prevalence for estimating prevalence
         cost_unit = factor(c_u))

d_dorf <- scenarios %>%
  subset(c_u == 0) %>% #--- cost regime is irrelevant for identification
  rowwise() %>%
  mutate(regime = 'Identification - Dorfman',
         opts = opts_dorf(p), #optimal pool size for identification (Dorfman) 
         optpp = 1 - (1-p)^opts, #optimal pool-level prevalence for identification (Dorfman)
         cost_unit = NA)

 

# Ignoring cost but comparing identification vs prevalence estimation ----

#optimal pool size
bind_rows(d_prev %>% subset(cost_unit == 0), d_dorf) %>%
  ggplot(aes(x = p, y = opts, colour = regime)) +
  geom_step() + 
  labs(x = 'Unit level prevalence',
       y = 'Optimal pool size') +
  scale_y_log10() +
  scale_x_log10(labels = scales::percent)

#optimal pool-level prevalence
bind_rows(d_prev %>% subset(cost_unit == 0), d_dorf) %>%
  ggplot(aes(x = p, y = optpp, colour = regime)) +
  geom_step() + 
  labs(x = 'Unit level prevalence',
       y = 'Optimal pool level prevalence') +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(labels = scales::percent)


# For different costs but for prevalence estimation only ----

# optimal pool size for given prevalence - prevalence estimation
d_prev %>%
  ggplot(aes(x = p, y = opts, color = cost_unit)) +
  geom_step() +
  labs(x = 'unit level prevalence',
       y = 'Optimal pool size') +
  scale_y_log10()+
  scale_x_continuous(labels = scales::percent)

# optimal pool-level prevalence
dmeans <- d_prev %>% group_by(cost_unit) %>% summarise(m = mean(optpp))
d_prev %>%
  ggplot(aes(x = p, y = optpp, color = cost_unit)) +
  geom_step() +
  geom_hline(data = dmeans %>% mutate(cost_unit = cost_unit),
             aes(yintercept = m, color = cost_unit), linetype = 2) +
  labs(x = 'unit level prevalence',
       y = 'Optimal Pool-level prevalence') +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)
