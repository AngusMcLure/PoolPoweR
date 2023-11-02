#When estimating prevalence, if catches at each site are random, should one use fixed pool sizes or fixed number of pools, or something else?


#simple example:
# Design prevalence was 1% and design correlation was 0.1.
# cost.location:cost.pool:cost.unit = 40:4:1.4
# Therefore plan was to do 4 pools of 4 mozzies each per location
# But 10 locations caught 8 mosquitoes and 10 locations caught 24 mosquitoes each (total catch is expected, but unevenly distributed)
# compare fisher information for:
#     4 pools at each location (pool sizes 2 and 6)
#     4 mozzies per pool (pool numbers of 6 and 2)

theta <- 0.01; rho <- 0.1

var.equal.N <- var.equal.s <- var.unpooled <- c()
for(k in 1:4){
  var.equal.N[k]  <- solve(10*(fi_pool_imperfect_cluster(k,4,theta,1,1,rho) + fi_pool_imperfect_cluster(8-k,4,theta,1,1,rho)))[1,1]
  var.equal.s[k]  <- solve(10*(fi_pool_imperfect_cluster(4,k,theta,1,1,rho) + fi_pool_imperfect_cluster(4,8-k,theta,1,1,rho)))[1,1]
  var.unpooled[k] <- solve(10*(fi_pool_imperfect_cluster(1,k*4,theta,1,1,rho) + fi_pool_imperfect_cluster(1,(8-k)*4,theta,1,1,rho)))[1,1]
}
var.equal.N
var.equal.s
var.unpooled

plot(1:k, var.equal.N)
plot(1:k, var.equal.s)
plot(1:k, var.unpooled)
plot(1:k, var.equal.s/var.equal.N)
var.equal.s/var.equal.N

cu <- 1.4; ct <- 4; cl <- 40;
optimise_s_prevalence (theta,cu,ct,cl,rho, N = 4,max.s = 250,interval = 0.2)
opt <- optimise_sN_prevalence(theta,cu,ct,cl,rho,       max.s = 50)

catch.mean <- opt$catch
catch.cv <- 1
catch.var <- (catch.mean * catch.cv)^2
catch.disp <- catch.mean^2/(catch.var - catch.mean); catch.disp;
qnbinom(c(0.025, 0.5,0.975), size = catch.disp, mu = catch.mean)


FI.equal.Ns <- fi_pool_imperfect_cluster(opt$s,opt$N,
                                        theta, 1, 1,
                                        rho,form = 'beta')
FI.equal.N <- fi_pool_imperfect_cluster_unequal(nb_catch(catch.mean, catch.var),
                                  pool_fixed_N(opt$N), theta, 1, 1,
                                  rho,form = 'beta',max.iter = 300)
FI.equal.s <-fi_pool_imperfect_cluster_unequal(nb_catch(catch.mean, catch.var),
                                  pool_max_s(opt$s), theta, 1, 1,
                                  rho,form = 'beta',max.iter = 300)
solve(FI.equal.Ns)
solve(FI.equal.N)
solve(FI.equal.s)



