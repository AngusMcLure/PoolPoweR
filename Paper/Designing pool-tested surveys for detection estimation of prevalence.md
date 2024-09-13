# Designing pool-tested surveys for detection estimation of prevalence

[TOC]

## Designing surveys for the estimation of prevalence

### Simple random surveys

#### Likelihood and Fisher information

We consider a survey where units are collected from a population using simple random sampling and tested in pools/groups with a binary test in order to estimate the prevalence ($\theta$) of positive units in the population. We assume that the pool/group sizes are fixed by design and the sensitivity and specificity of the test are known. Consider $n$ independent Bernoulli random variables with mean $\theta$ that are randomly assigned to pools/groups of $K$ different sizes. Let the $\underline{s} = [s_1,\dots,s_K]$ be the vector of pool/group sizes and $\underline{N} = [N_1, \dots, N_K]$ be the number of pools/groups of corresponding sizes such that $\underline{s}\cdot\underline{N} = n$. Consider the result of an imperfect test with known sensitivity, $0<\varphi\leq 1$, and specificity, $0<\psi\leq 1$ applied to each pool that (imperfectly) indicates the presence of any positive units in the pool/group. We assume that over the range of pool/groups sizes considered, test sensitivity and specificity are independent of pool/group size and the sensitivity is the same whether one or multiple units in a pool are positive.  In other words we model the probability that a pool of size $s$ tests positive as
$$
\begin{align*}
\phi_s(\theta)
:&= (1-\psi)(1-\theta)^s + \varphi (1 - (1-\theta)^s)\\
&=\varphi + (1-\varphi-\psi)(1-\theta)^s
\end{align*}
$$

Let $\underline{Y} = [Y_1,\dots, Y_K]$ be the number of positive results for pools/groups of each size. The likelihood for the $\theta$ is then
$$
L(\theta|\underline{y}) = P(\underline{Y} = \underline{y}|\theta,\underline{N},\underline{s},\varphi,\psi) = \prod_{k = 1} ^K {N_k \choose y_k} {\phi_{s_k}(\theta)}^{y_k} (1-\phi_{s_k}(\theta))^{N_k-y_k}
$$
Maximum likelihood estimation of the prevalence $\theta$ must be performed numerically for all cases except where pool/group sizes are equal.

Considering the Fisher information for $\theta$ for the simple case where $\underline{N} = 1$ and $\underline{s} = s$ it is simple to show that this is
$$
\begin{align*}
I(\theta|1, s,\varphi,\psi)
= \frac{s^2(1-\theta)^{2s-2}(1-\varphi - \psi)^2}{\phi_{s}(\theta)(1-\phi_{s}(\theta))}\\
\end{align*}
$$
Since the results of pools/groups are independent given $\theta$, the Fisher information for $\theta$ from such a survey is the sum over Fisher information from each pool:
$$
\begin{align*}
I(\theta|\underline{N}, \underline{s},\varphi,\psi)
&= \sum_{k = 1}^K N_k I(\theta|1,s_k,\varphi,\psi)\\
&= (1-\varphi-\psi)^2\sum_{k = 1}^K N_k \frac{s_k^2(1-\theta)^{2s_k-2}}{\phi_{s_k}(\theta)(1-\phi_{s_k}(\theta))}\\
\end{align*}
$$
For large $n$ the maximum likelihood estimator of $\theta$ is normally distributed with variance $I(\theta|\underline{N},\underline{s},\varphi, \psi)^{-1}$.

The Fisher information is well defined for all $\theta$ if the test is imperfectly specific, but diverges as $\theta \to 0$ for a perfectly specific test:
$$
\lim_{\theta \to 0} I(\theta |1,s,\varphi, \psi) =
\begin{cases}
 \frac{s^2 (1-\varphi - \psi)^2}{(1-\psi)\psi}  &  0<\psi<1\\
\infty & \psi =1
\end{cases}
$$

The behaviour of the Fisher information as $\theta \to 1$ depends on test sensitivity, $\varphi$, with different forms for $s =1$, $s = 2$, and $s>2$:
$$
\lim_{\theta \to 1} I(\theta|1,s,\varphi,\psi) =
\begin{cases}
\frac{ (1-\varphi - \psi)^2}{(1-\varphi)\varphi} &  s=1, \ \ \ 0<\varphi<1\\
0 & s>1, \ \ \  0<\varphi<1\\
\infty & s = 1, \ \ \ \varphi = 1 \\
4\psi & s= 2, \ \ \ \varphi = 1 \\
0 & s>2, \ \ \ \varphi = 1 \\
\end{cases}
$$
The point to note here is that a group test yields more than an individual test when prevalence is low but far less information than individual tests when prevalence is high.

#### Design effect

To make comparisons between designs with the same number of vectors but different pooling/grouping strategies it can be helpful to consider the *unit Fisher information*, or Fisher information per unit in the sample $I(\theta|\underline{N}, \underline{s},\varphi,\psi)/n$, where we recall that $n = \underline{N}\cdot \underline{s}$ If $I(\theta|\underline{N},\underline{s}, \varphi, \psi)/n < I(\theta|1,1, \varphi, \psi)$ then pooling/grouping reduces the unit Fisher information and therefore the total Fisher information to be obtained from a given number of units. Therefore to achieve a target variance for the estimator, pool/group testing will require a larger sample size. The *design effect* for a pool/group testing design is then
$$
D(\theta,\underline{N},\underline{s},\varphi,\psi) = \frac{I(\theta|1,1,\varphi,\psi) \underline{N}\cdot\underline{s},}{I(\theta|\underline{N},\underline{s},\varphi,\psi)}.
$$
Note that scaling $\underline{N}$ by a multiplicative factor doesn't change the design effect. Consequently, for surveys with equal sized pools, the design effect is independent of $N$:
$$
D(\theta,N,s,\varphi,\psi)
= \frac{I(\theta|1,1,\varphi,\psi) Ns,}{I(\theta|N,s,\varphi,\psi)}
= \frac{I(\theta|1,1,\varphi,\psi) Ns,}{I(\theta|1,s,\varphi,\psi) \ N}
= \frac{\phi_s(\theta)(1-\phi_s(\theta))}{s(1-\theta)^{2s-2} \phi_1(\theta)(1-\phi_1(\theta))}.
$$
 For a perfectly sensitive and specific test this simplifies to
$$
D(\theta,N,s,1,1) = \left[ \theta^{-1} - 1\right]\left[(1-\theta)^{-s} - 1\right].
$$
For perfectly specific tests, $D\geq 1$ with equality only when $s = 1$ or as $\theta$ approaches 0, i.e. pooling increases required sample size to achieve a given variance. When prevalence is low the effect is smaller and may be offset by the reduced testing costs associated with pool/group testing. However, if the test is not perfectly specific and prevalence ($\theta$) is small group testing may *increase* the Fisher information (i.e. $D<1$). In fact for equal sized pools
$$
\lim_{\theta \to 0} D(\theta, N, s, \varphi,\psi) =
\begin{cases}
1/s & \psi <1\\
1 & \psi = 1
\end{cases}
$$
In other words, for an imperfect test and prevalence sufficiently close to 0, testing single pool of $s$ units renders the nearly same amount of information as $s^2$ individually tested units. ***Figure []*** illustrates this counter-intuitive result by considering the design effect with four different values of specificity, five different pool/groups sizes, and a range of prevalences.

One can make intuitive sense of this result, by considering the positive predictive value of a pooled test:
$$
PPV(\theta, s, \varphi, \psi) = \frac{\varphi(1-(1-\theta)^s)}{\varphi(1-(1-\theta)^s) + (1-\psi)(1-\theta)^s}.
$$
For a perfectly specific test, positive predictive value is always 1. For an imperfectly specific test, $\theta \approx 0$ and $s \ll \frac{1}{\theta}$, $PPV \approx \frac{\varphi s \theta}{1-\psi}$. In other words, positive predictive value approaches 0 for low prevalence, but is increased approximately s-fold when using modestly sized pools. Furthermore, unless $\theta = 0$, PPV increases monotonically with $s$ and $\lim_{s\to \infty} PPV = 1$.  At low prevalence the high proportion of false positives from an imperfectly specific test makes estimation of prevalence very difficult, but the proportion of false positives can be reduced by testing in pools.

Design effect increases with $\theta$ but most rapidly when sensitivity is low and pools/groups are large.

<img src="./Figures/FigDesignEffectSupp.png" style="zoom:50%;" /> 

***Figure []*** Design effect comparing pool-testing to individual testing for different sized pools, different prevalence, and a range of test sensitivity and specificity. Note that in comparing pooled and individual testing, both the pooled and individual tests are assumed to have the same sensitivity and specificity.

#### Optimising simple random surveys with pooled testing

However, when designing a pooling scheme, the relative cost/effort of collecting and testing vectors needs to be considered alongside the design effect. While the total cost of the survey may also include fixed overhead costs, we consider only the costs that depend on sample size and number of pools/groups as these are the ones that can be minimised by appropriate survey design. Let $c_u$ be the marginal cost of sampling one unit and $c_p$ be the marginal cost associated with handling and testing one pool/group (independent of the number of vectors in the pool/group). The total marginal cost of handling, testing, and sampling units for a pool of size $s$ is then $s c_u + c_p$. For a pool/group of size $s$, the marginal cost of a 'unit' of Fisher information is therefore
$$
m(s|\theta, \varphi, \psi, c_u, c_p) = \frac{c_u s + c_p}{I(\theta|1,s, \varphi, \psi)}.
$$

Choosing pool size ($s$) that minimizes $m$ gives the most cost-effective way to reduce the variance of prevalence estimate. This optimum depends on the test characteristic ($\varphi, \psi$), costs ($c_u, c_p$), and prevalence ($\theta$). Note that multiplying costs by a positive constant does not change the optimum $s$, so we need only consider the ratio of the costs $c_p/c_u$ and the special case where $c_u =0$. Extensive numerical experiments suggests that any reasonable prevalence, test characteristics, and costs, $m$ has unique minima $s^*$ for real $s \geq 1$. If generally true, standard optimisation methods can be used to identify $s^*$, and the integer value of $s$ that minimises $m$ is either $\lceil s^* \rceil$ or $\lfloor s^* \rfloor$. Though we do not have a proof for the most general case, we show that for a perfectly sensitive test $m'(s) = 0$ has at most one solution in $[1,\infty)$ and that therefore local minima in $[1,\infty)$ are global minima.
$$
m'(s) = \frac{c_u I - (c_u s + c_p)\frac{\part I}{\part s}}{I^2}
$$
 Since $\frac{c_u s + c_p}{I} > 0$, solutions for $m'(s) = 0$ must also satisfy
$$
\begin{align*}
0
&= \frac{c_u}{c_u s + c_p} - \frac{\part\log(I)}{\part s} \\
&= \frac{c_u}{c_u s + c_p} - \frac{2}{s} - 2\log(1-\theta) + \frac{\part \phi_s}{\part s} \left[ \frac{1}{\phi_s} - \frac{1}{1-\phi_s} \right] \\
&= \frac{c_u}{c_u s + c_p} - \frac{2}{s} - 2\log(1-\theta) + \log(1-\theta)(\phi_s-\varphi) \left[ \frac{1}{\phi_s} - \frac{1}{1-\phi_s} \right] \\
&= \frac{c_u}{c_u s + c_p} - \frac{2}{s} - \log(1-\theta)\left[\frac{\varphi}{\phi_s} + \frac{1-\varphi}{1-\phi_s} \right] \\
\end{align*}
$$
For a perfectly sensitive test ($\varphi = 1$)  this reduces to
$$
\begin{align*}
0

&=  \frac{c_u}{c_u s + c_p} - \frac{2}{s} - \frac{\log(1-\theta)}{1 - \psi (1-\theta)^s}\\
\end{align*}
$$
We prove that the RHS of the above is monotonic increasing in $s$ for $s\geq 1$ and therefore has at most one root on the same domain.
$$
\begin{align*}
\frac{\part}{\part s}\left[\frac{c_u}{c_u s + c_p} - \frac{2}{s} - \frac{\log(1-\theta)}{1 - \psi (1-\theta)^s}\right]
&= -\frac{c_u^2}{(c_u s + c_p)^2} + \frac{2}{s^2} - \psi\frac{(\log(1-\theta))^2(1-\theta)^s}{(1 - \psi (1-\theta)^s)^2}\\
&\geq \frac{1}{s^2} - \psi\frac{(\log(1-\theta))^2(1-\theta)^s}{(1 - \psi (1-\theta)^s)^2}\\
&\geq \frac{1}{s^2} - \frac{(\log(1-\theta))^2(1-\theta)^s}{(1 - (1-\theta)^s)^2}\\
\end{align*}
$$
where the first inequality follows from $s,c_u,c_p \geq 0$ (with equality if $c_p = 0$), the second inequality follows from $\psi \leq 1$. It remains to be shown that the final term is $\geq 0$ for $s\geq 1$ or equivalently that $h(s) :=\frac{(1 - (1-\theta)^s)^2}{(1-\theta)^s} - s^2(\log(1-\theta))^2 \geq 0$ for $s\geq 1$. We proceed by proving that $h(1) \geq 0$, $h'(1) \geq 0$, and $h''(s) \geq 0$ for $s \geq 1$.
$$
\begin{align*}
h(1)
&= \frac{(1 - (1-\theta))^2}{1-\theta} - (\log(1-\theta))^2 \\
&= \frac{\theta^2}{1-\theta} - (\log(1-\theta))^2 \\
&\geq \frac{\theta^2}{1-\theta} - \frac{\theta^2}{1-\theta} = 0
\end{align*}
$$
with the inequality following from $-\log(1-\theta) \leq \frac{\theta}{\sqrt{1-\theta}}$ for $0 \leq \theta < 1$ with equality only if $\theta = 0$.
$$
\begin{align*}
h'(s) &= \left[(1-\theta)^s - (1-\theta)^{-s}\right] \log(1-\theta) - 2s(\log(1-\theta))^2 \\
h'(1)
&= \left[(1-\theta) - \frac{1}{1-\theta}\right]\log(1-\theta) -2(\log(1-\theta))^2\\
&=  - \theta \frac{2-\theta}{1-\theta} \log(1-\theta) -2(\log(1-\theta))^2\\
&\geq  \frac{2-\theta}{1-\theta} (\log(1-\theta))^2 -2(\log(1-\theta))^2\\
&= (\log(1-\theta))^2 \frac{\theta}{1-\theta} \geq  0\\
\end{align*}
$$
where the first inequality follows from $\log(1-\theta) \leq -\theta$ for $\theta < 1$. Finally
$$
\begin{align*}
h''(s)
&= \left[(1-\theta)^s + (1-\theta)^{-s}\right](\log(1-\theta))^2 - 2(\log(1-\theta))^2 \\
&= \left[(1-\theta)^s + (1-\theta)^{-s} - 2\right](\log(1-\theta))^2\\
&\geq \left[1-s\theta + 1+s\theta - 2\right](\log(1-\theta))^2\\
&= 0
\end{align*}
$$
where the inequality follows from two uses of the Bernoulli inequality, i.e. $(1-\theta)^r \geq 1-r\theta$ for $r \notin [0,1]$ and $ \theta <1$.

***Figure []***  shows the result of numerical optimisation of $s$ for a range of test characteristics (sensitivity and specificity), prevalences, and cost ratios. Increasing prevalence, decreasing sensitivity, increasing specificity, and decreasing $c_p/c_u$ generally favours smaller pools. In the case of perfect specificity and $c_p/c_u \ll 1$, then $m(s) \approx c_u\frac{s}{I(s)} > c_u \frac{1}{I(1)}$ so individual tests are always most efficient regardless of sensitivity or the prevalence of the marker in the vector population. In practice no test is perfectly specific; however, even if a test is 99.9% specific and the cost of testing and handling a pool is much less than collecting units ($c_p = 0$ or $c_p/c_u \ll 1$), the optimal pool size will be greater than one if prevalence is sufficiently low (red line **Figure []**), e.g. $s=2$ if prevalence is 2%,  $s=4$ if prevalence is 1%, or  $s=8$ if prevalence is 0.5%. The optimal pool/group size increases further with increasing $c_p/c_u$. For instance if specificity is 99.9% and a $c_p/c_u = 5$  then the optimal pool sizes are $s = 6$ if  $\theta = 10\%$,  $s=17$ if  $\theta = 2\%$,  $s=26$ if $\theta =1\%$, or  $s=38$ if $\theta= 0.5\%$.

However the optimal pool size does not increase without bound for increasing $c_p/c_u$. In the limiting case where $c_u = 0$ (denoted as $c_p/c_u = \infty$ and purple line in **Figure []**), choosing $s$ that minimises the cost of information becomes equivalent to choosing $s$ that maximises the total Fisher information from a pool. The optimal $s$ in this case depends very little on test specificity and only modestly on test sensitivity (increasing with increasing sensitivity). Whenever the cost of sampling units is effectively zero, such as when collection has already occurred and therefore can be considered a 'sunk cost', this provides an upper bound to the size of pools that should be considered.

While specificity has a dramatic effect on the optimal $s$ especially for low $\theta$, there is no matching dramatic effect for sensitivity. If sensitivity were to decline substantially with increasing pool size it may have a greater effect of choice of pool/group size but we do not consider this case. However, reducing sensitivity does decrease Fisher information 

![](./Figures/FigOptimals.png) 

***Figure []*** Pool size that optimises unit cost of Fisher information as a function of prevalence (x-axis), sensitivity (rows), specificity (columns) and the relative cost of testing a pool ($c_p$) and collecting a unit ($c_u$). $c_p/c_u = 0$ and $c_p/c_u =$ Inf indicate the two extreme scenarios where testing pools or collecting units are free (or a sunk cost).   

### Cluster surveys - single level of clustering
If prevalence survey that will used pooled testing adopts a clustered or hierarchical sampling frame, this needs to be taken into account when determining sampling size and target pooling strategy. Correlation of the presence/absence of the marker in units at a location induces correlation of the outcomes of pools, changing the effective sample size and information to obtained from a given number of units. With positively correlated outcomes the total number of units required to achieve a certain level of power or Fisher information increases. For individually tested prevalence surveys a common approach is to calculate a design effect, which gives a multiplicative factor for the sample size that is required to achieve a given level of power or information with a clustered survey vs a simple random survey. Consider the survey where a total of $n$ units are sampled from $J$ locations, with $n/J$ units tested per location. Let $V_{ij}$ be correlated Bernoulli random variables indicating the presence of the marker in unit $i$ from location $j$, where units are correlated only with other units from the same location i.e.
$$
Corr(V_{ij},V_{i'j'}) = \begin{cases}
1, & i = i', j = j' \\
\rho, & i \neq i', j = j' \\
0, &  j \neq j' \\
\end{cases}
$$
Then the standard population prevalence estimator
$$
\hat{\theta} = \frac{1}{n} \sum_i\sum_j V_{ij}
$$
has standard error
$$
se(\hat{\theta}) = \sqrt{\dfrac{\hat{\theta}(1-\hat{\theta})D}{n}}  = \sqrt{\dfrac{\hat{\theta}(1-\hat{\theta})}{n_{eff}}}
$$
with design effect $D = 1 + \rho (n/J - 1)$ and effective sample size $n_{eff} = n/D$. Note that the design effect calculated with this method depends only on the number of units per location ($n/J$) and the pairwise correlation between units ($\rho$) and doesn't require specific knowledge of the joint distribution of units.

Though pool/group tests can be modelled as correlated Bernoulli trials we cannot use the above approach to yield a design effect that depends only on the pairwise correlation between units. The essential difficulty is that individual units within pools/groups will also be correlated. To illustrate why the above approach cannot be used, consider the case of common pool/group size, $s$, and assume that pools/groups contain only units from a single location. Then let $X^s_{ij}$ be the correlated Bernoulli random variables indicating the presence of marker in pool/group  $i$ from location $j$.  Let $\theta_s = E\left[X^s_{ij}\right]$ be the prevalence of positive pools. As with individual units, pools are then only correlated with other pools from the same location:
$$
Corr(X_{ij},X_{i'j'}) = \begin{cases}
1, & i = i', j = j' \\
\rho_s, & i \neq i', j = j' \\
0, &  j \neq j' \\
\end{cases}
$$
The pool-level prevalence, $\theta_s$, can be estimated with
$$
\hat{\theta}_s = \frac{s}{n}\sum_i\sum_j X_{ij}
$$
and has standard error
$$
se(\hat{\theta}_s) = \sqrt{\dfrac{\hat{\theta}_s(1-\hat{\theta}_s) \left[1 + \rho_s (\frac{n}{Ms} - 1)\right]}{n/s}}.
$$
However, we are trying to estimate the prevalence of positive units ($\theta$), not positive pools/groups $\theta_s$. If $h(\theta_s) = \theta$, then $h(\hat{\theta_s})$ is an estimator of $\theta$ with standard error
$$
se(h(\hat{\theta}_s)) = h'(\hat{\theta}_s)\sqrt{\dfrac{\hat{\theta}_s(1-\hat{\theta}_s) \left[1 + \rho_s (\frac{n}{Ms} - 1)\right]}{n/s}}.
$$

However, $h$ and $\rho_s$ both depend not only pairwise correlations ($\rho$), but on the form of the joint distribution of $s$ and $2s$ units:
$$
\begin{align*}
h^{-1}(\theta)
&= \theta_s\\
&= P(X_{ij} = 1|\theta,\rho)\\
&= 1 - P(V_{1j}=\dots=V_{sj} = 0|\theta,\rho)\\
h'(\theta_s) &= -\left[\frac{\part}{\part \theta} \left. P(V_{1j}=\dots = V_{sj} = 0|\theta, \rho)\right\rvert_{\theta = h(\theta_s)}\right]^{-1}
\end{align*}
$$


$$
\begin{align*}
\rho_s
&= \frac{E[X_{ij} X_{i'j}] - E[X_{ij}] E[X_{i'j}]}{\sqrt{Var[X_{ij}] Var[X_{i'j}]}} \\
&= \frac{P(X_{ij} =  X_{i'j} = 1) - E[X_{ij}]^2}{Var[X_{ij}]} \\
&= \frac{1 - P(V_{1,j} = \dots = V_{s,j} = 0) - P(V_{1+s,j} =\dots = V_{2s,j} = 0) + P(V_{1,j} = \dots = V_{2s,j} = 0) - E[X_{ij}]^2}{Var[X_{ij}]} \\
&= \frac{1 - 2P(V_{1,j} = \dots = V_{s,j} = 0) + P(V_{1,j} = \dots = V_{2s,j} = 0) - E[X_{ij}]^2}{Var[X_{ij}]} \\
\end{align*}
$$

Calculation of the above probabilities requires that we specify the correlation structure for $s$ and $2s$ units from the same location. If we want to compare design effects for different pool/group sizes we must specify the correlation structure for all $s$ we wish to consider. Since we must specify the higher degree correlation structure we use parametric means to specify the full joint distribution of $\{V_{ij}\}$, allowing us to calculate Fisher information for very general group/pool testing regimes.

A general approach to specify the joint distribution of $\{V_{ij}\}$ is to model the prevalence as variable across locations, with units being independent conditioned on the location.  Let $\Theta_j$ represent the prevalence at site $j$. We model $\{\Theta_j\}$ as i.i.d random variables with mean $\theta$, support on $[0,1]$, and density $f_\Theta(p)$ and model units such that
$$
V_{ij}|\Theta_j \sim Bern(\Theta_j),
$$
where $\{V_{ij}|\Theta_j\}$ are independent. Then, if $n_j$ is the number of units from site $j$ then
$$
\sum_{i=1}^{n_j} V_{ij} |\Theta_j \sim Bin(\Theta_j,n_j).
$$
This set up preserves $E[V_{ij}] = \theta$. Units from different sites are independent and therefore uncorrelated, while correlation between two distinct units from the same site, $\rho$, is
$$
\begin{align*}
\rho := Corr(V_{ij},V_{i'j}) &= \frac{E[V_{ij}V_{i'j}] - E[V_{ij}][V_{i'j}]}{\sqrt{Var[V_{ij}]Var[V_{i'j}]}} \\
&=  \frac{P(V_{ij} = V_{i'j} = 1) - E[V_{ij}]^2}{Var[V_{ij}]} \\
&=  \frac{E[\Theta_j^2] - E[\Theta_j]^2}{Var[V_{ij}]} \\
&= \frac{Var[\Theta_j]}{\theta (1-\theta)}.
\end{align*}
$$

We can use any two-parameter family with support on $[0,1]$ to model the distribution of prevalences across sites. We consider three examples: the beta distribution, the logit-normal distribution, and a simple discrete distribution. However we first consider some general results.

Note that $\rho$ follows the general definition of the correlation of two random variables. This should not be confused with the related concept of the intra-cluster (intraclass) correlation coefficient (ICC) which has been calculated in many many ways (see for instance Chakraborty and Hossain (2018) or Goldstein et al (2002)). Notably, the *estat icc* command for estimation of ICCs from mixed and random effects generalised linear models in STATA takes the approach of partitioning variance on the real/propensity scale rather than the probability/prevalence scale. When prevalence is close to 50% and correlation is low, then the ICC calculated in STATA is similar to the definition of $\rho$. However, if prevalence is low (as is usually the case when using pooled testing) the estimate of the ICC will be much larger (possibly by more than an order of magnitude) that $\rho$. 

#### General properties of the likelihood function and Fisher information

##### Likelihood

Recall that $X_{ij}$ are correlated Bernoulli random variables indicating the presence of marker in pool/group $i$ from location $j$. As infection status of vectors sampled from different locations are independent, the outcome of test on pools from different locations are also independent. Therefore we can focus our attention on defining the joint distribution of pools from the same location, dropping references to location $j$ for notational simplicity. Consider first the case where there are $N$ pools at a location each with same number of units $s$. Let $Y$ be the number of positive pools i.e. $Y = \sum_i X_{i}$, and recall that the pool positivity function  is $\phi_s(p) = (1-\varphi - \psi) (1-p)^s + \varphi$. Then
$$
L(\theta, \rho|y,N,s,\varphi,\psi) = {N \choose y} \int_0^1 f_\Theta(p) {\phi_s(p)}^y (1-\phi_s(p))^{N-y} dp
$$

In the case where pools/groups from a single location are different sizes, the expressions for the likelihood are similar. Consider a sample with $K$ unique pool/group sizes $\underline{s} = [s_1,\dots,s_K]$. Let $\underline{N} = [N_1,\dots, N_K]$ be the number of pools/groups of each size and  $\underline{Y}=[Y_1,\dots,Y_K]$ be the number of positive pools of each size. Then the likelihood is
$$
L(\theta, \rho|\underline{y},\underline{N},\underline{s},\varphi,\psi) = \int_0^1 f_\Theta(p) \prod_{k=1}^K {N_k \choose y_k} {\phi_{s_k}(p)}^{y_k} (1-\phi_{s_k}(p))^{N_k-y_k} dp.
$$

Consider the case where $\underline{s} \cdot \underline{N} = 2$, i.e. where a total of two units are sampled at each site and they are either tested individually ($s = 1$, $N = 2$) or in a single pool ($s = 2$, $N = 1$). In either case the likelihood (and therefore the Fisher information) does not depend on the distributional form of $\Theta_j$, but only on $\theta$ and $\rho$. The intuition here is that is $\theta$ and $\rho$ uniquely define the joint distribution of two correlated Bernoulli random variables (though not of three or more), so whenever there are only two units per sampling site the likelihood function should only depend on these parameters. To see this note that for $0\leq y \leq N$ and general $s$ and $N$, the term $\phi_s(p)^y (1-\phi_s(p))^{N-y}$ is a polynomial in $p$ of order $sN$, with coefficients that are functions of $s,N,\varphi,$ and $\psi$ and do not depend on the distributional form of $\Theta_j$. Therefore when $sN \leq 2$, we can write the likelihood as
$$
\begin{align*}
L(\theta, \rho|y)
&= \int_0^1 f_\Theta(p) (c_0 + c_1 p + c_2 p^2) dp\\
&= c_0 + c_1 E[\Theta_j] + c_2 E[\Theta_j^2] \\
&= c_0 + c_1 \theta + c_2 (\rho \theta(1-\theta) + \theta^2),
\end{align*}
$$

where $c_0, c_1, c_2$  do not depend on the distributional form of $\Theta_j$. For $\underline{s} \cdot \underline{N} > 2$ the likelihood will depend on higher moments of $\Theta_j$ and therefore specifying $\theta$ and $\rho$ will not be sufficient.

##### Fisher information matrix

With an expression for the likelihood we can now write out the Fisher information matrix. The number of pools ($\underline{N}$) of each size ($\underline{s}$) are fixed by design, and we assume that the sensitivity ($\varphi$) and specificity ($\psi$) are known, therefore the unknown parameters are only the population mean ($\theta$) and the pairwise correlation of units ($\rho$). Let $\Xi_{\underline{N}}$ be the set of all possible outcomes $\underline{y}$ given pool numbers $\underline{N}$, i.e. $\{0,\dots,N_1\} \times \dots \times \{0,\dots,N_K\}$. The Fisher information matrix is then
$$
I(\theta,\rho|\underline{N}, \underline{s}, \varphi, \psi) = \sum_{\underline{y} \in \Xi_{\underline{N}}} L(\theta, \rho|\underline{y})^{-1}
\begin{bmatrix}
 L_\theta(\theta, \rho|\underline{y}) ^2 &
 L_\theta(\theta, \rho|\underline{y}) L_\rho(\theta, \rho|\underline{y})\\
 L_\theta(\theta, \rho|\underline{y}) L_\rho(\theta, \rho|\underline{y}) &
 L_\rho(\theta, \rho|\underline{y}) ^2\\
\end{bmatrix}
$$

where we have suppressed the dependance of the right hand side on $\underline{N}, \underline{s}, \varphi$, and $\psi$ for notational compactness. As the pool probability function $\phi_s(p)$ does not depend on $\theta$ or $\rho$, if $f_\Theta$ is differentiable with respect to $\theta$ and $\rho$, the derivatives of the likelihood are simply

$$
L_\theta(\theta, \rho|\underline{y}) = \int_0^1 \frac{\part f_\Theta(p|\theta, \rho)}{\part \theta} \prod_{k=1}^K {N_k \choose y_k} {\phi_{s_k}(p)}^{y_k} (1-\phi_{s_k}(p))^{N_k-y_k} dp
$$
and
$$
L_\rho(\theta, \rho|\underline{y}) = \int_0^1 \frac{\part f_\Theta(p|\theta, \rho)}{\part \rho} \prod_{k=1}^K {N_k \choose y_k} {\phi_{s_k}(p)}^{y_k} (1-\phi_{s_k}(p))^{N_k-y_k} dp
$$

##### Large sample variance of MLE

As units and results of pool/group tests from different sampling locations are assumed independent, the total Fisher information from more than one sampling location is the sum of the Fisher information from each location. That is for $J$ locations with sampling and pooling designs ($\underline{N}_j$ and $\underline{s}_j$) the total Fisher information is
$$
I(\theta, \rho|\{\underline{N}_j\},\{\underline{s}_j\},\varphi,\psi) = \sum_{j=1}^J I(\theta, \rho|\underline{N}_j,\underline{s}_j,\varphi,\psi)
$$
For identical sampling and pooling designs at each location this is simply $J \cdot I(\theta, \rho|\underline{N},\underline{s},\varphi,\psi)$. For large $J$ the maximum likelihood estimator of $[\theta, \rho]^T$, is approximately multivariate normal distributed with variance-covariance matrix $\Sigma = I(\theta, \rho|\underline{N}, \underline{s}, \varphi, \psi)^{-1}/J$. In many cases we are primarily interested in population prevalence, and within-site correlations are considered a nuisance parameter. In such cases we are seeking to choose survey designs that minimise $\Sigma_{1,1}$, or achieve a target value of $\Sigma_{1,1}$ while minimising cost or effort.

With this in mind we can define the design effect for a cluster survey with samples collected in pools as
$$
D(\theta, \rho, \underline{N},\underline{s},\psi, \varphi) \equiv I(\theta|1,1,\varphi,\psi) \ \underline{N} \cdot \underline{s} \ \Sigma(\theta, \rho|\underline{N}, \underline{s}, \varphi, \psi)_{1,1}.
$$


##### Non-identifiability with $N=1$

Consider the case where $N = 1$. Intuitively, if there is only a single/group per sampling location, then we should expect not to be able to simultaneously estimate prevalence ($\theta$) and correlation between units within sites ($\rho$), i.e. that $I(\theta, \rho|1,s,\varphi,\psi)$ is non-invertible. This is indeed the case as can be seen by noting that $\Xi_1 = \{0,1\}$ and $L(\theta, \rho|1) = 1- L(\theta, \rho|0)$ and therefore $L_\theta(\theta, \rho|1) = -L_\theta(\theta, \rho|0)$ and $L_\rho(\theta, \rho|1) = -L_\rho(\theta, \rho|0)$. With these identities we can expand out the above sum as
$$
\begin{align*}
I(\theta,\rho|1,s,\varphi,\psi) &= 
L(\theta, \rho|0)^{-1}
\begin{bmatrix}
 L_\theta(\theta, \rho|0) ^2 &
 L_\theta(\theta, \rho|0) L_\rho(\theta, \rho|0)\\
 L_\theta(\theta, \rho|0) L_\rho(\theta, \rho|0) &
 L_\rho(\theta, \rho|0) ^2\\
\end{bmatrix} + \\
& \ \ \ \ \ \
L(\theta, \rho|1)^{-1}
\begin{bmatrix}
 L_\theta(\theta, \rho|1) ^2 &
 L_\theta(\theta, \rho|1) L_\rho(\theta, \rho|1)\\
 L_\theta(\theta, \rho|1) L_\rho(\theta, \rho|1) &
 L_\rho(\theta, \rho|1) ^2\\
\end{bmatrix} \\
&= (L(\theta, \rho|0)^{-1} + L(\theta, \rho|1)^{-1})
\begin{bmatrix}
 L_\theta(\theta, \rho|0)^2  &
 L_\theta(\theta, \rho|0) L_\rho(\theta, \rho|0)\\
 L_\theta(\theta, \rho|0) L_\rho(\theta, \rho|0) &
 L_\rho(\theta, \rho|0) ^2\\
\end{bmatrix}
\end{align*}
$$

which clearly has a determinant of zero.

##### Limits for small $\rho$

[I don't have any proofs in this section! I have a sketch of a series of conjectures which I believe to be true based on numerical experiments and intuition, but I haven't spent much time on any of them. Note of course that since this section is about limits, the numerical examples don't prove anything!]

As intracluster correlation $\rho$ approaches 0, the distribution of $\Theta_j$ approaches a point mass at $\theta$, i.e. no variation in prevalence between sampling locations. Intuitively we may expect that as $\rho$ approaches 0, the variance of the MLE of $\theta$ would also decrease. However in general the variance from a cluster survey with pooled testing may be higher, lower, or equal to the variance from a simple random survey with pools of the same size and total sample size.

**Conjecture 1**: For all $s$ the likelihood from a cluster sampling converges to the likelihood from simple random sampling, i.e. $\lim_{\rho \to 0}L(\theta, \rho|\underline{y}, N, s) = L(\theta|\underline{y},N,s)$. Moreover the Fisher information for $\theta$ from a cluster survey converges to the Fisher information from a simple random survey, i.e.  $\lim_{\rho \to 0}I(\theta, \rho|N,s)_{1,1} = I(\theta|N,s)$. Convergence may be from above or below. In other words for a given sample size and pooling setup ($N$ and $s$) and prevalence ($\theta$), a cluster survey may result in *more or less* Fisher information for $\theta$ than a simple random survey [I have numerical examples suggesting both cases].

**Corollary to Conjecture 1:** We have the following inequality
$$
\begin{align*}
\lim_{\rho \to 0} \Sigma(\theta, \rho | N,s)_{1,1}
&\equiv \lim_{\rho \to 0}\frac{1}{I(\theta,\rho|N,s)_{1,1} - I(\theta,\rho|N,s)_{1,2}^2/I(\theta,\rho|N,s)_{2,2}}\\
&\geq \lim_{\rho \to 0} \frac{1}{I(\theta,\rho|N,s)_{1,1}}\\
&= \frac{1}{I(\theta|N,s)},
\end{align*}
$$
with equality if $\lim_{\rho \to 0}I(\theta, \rho|N,s)_{1,2}^2/I(\theta,\rho|N,s)_{2,2} = 0$. However, $\lim_{\rho \to 0}\Sigma_{1,1}$ need not converge. In a subsequent section, we provide a numerical example where $\Sigma_{1,1}$ appears to grow without bound as $\rho$ approaches 0.

**Conjecture 2**: The above is an equality if:

* $s = 1$ for all $N$ and $\theta$,
* $s > 1$ for all $N$ in limit as $\theta \to 0$.


#### Likelihood and Fisher information for specific distributions of $\Theta_j$

Now we consider three two-parameter families of distributions that could be used to model $\Theta_j$. We parameterise each according to $\theta$ and $\rho$, showing correspondences to standard parametrisations where applicable. All distribution reduce to a point mass a $\theta$ as $\rho \to 0$. For each distribution we will begin by writing out the resultant likelihoods for the $\theta$ and $\rho$ given the number of positive pools. We then discuss the how to calculate the Fisher information, though closed form results are only available for the simplest cases.

#####  Simple discrete distribution for $\Theta_j$

Perhaps the simplest two parameter model for $\Theta_j$ is the following discrete distribution
$$
f_\Theta(p) = \delta_0(p) \rho(1-\theta) + \delta_\theta(p) (1-\rho) + \delta_1(p) \rho \theta.
$$
Under this distribution, the prevalence at a site is either equal to the population mean prevalence ($\theta$), exactly 0 (absent), or exactly 1 (ubiquitous). The correlation structure this imposes on units from the same location is equivalent to the correlated Binomial model proposed by Luceño (1995) and discussed by Diniz (2010). The resultant likelihood function is
$$
\begin{align*}
L(\theta, \rho|\underline{y})
&= \prod_{k = 1}^K{N_k \choose y_k} \left[\rho(1-\theta)\prod_{k = 1}^K(1-\psi)^{y_k} \psi^{N_k-y_k} + \rho\theta\prod_{k = 1}^K\varphi^{y_k}(1-\varphi)^{N_k-y_k} + (1-\rho)\prod_{k = 1}^K\phi_{s_k}(\theta)^{y_k}(1-\phi_{s_k}(\theta))^{N_k-y_k}\right] \\
& = \prod_{k = 1}^K{N_k \choose y_k} \left[\rho(1-\theta)(1-\psi)^{y^*} \psi^{N^*-y^*} + \rho\theta\varphi^{y^*}(1-\varphi)^{N^*-y^*} + (1-\rho)\prod_{k = 1}^K\phi_{s_k}(\theta)^{y_k}(1-\phi_{s_k}(\theta))^{N_k-y_k}\right]. \\
\end{align*}
$$
where $y^* = \sum_{k=1}^K y_k$ and $N^* = \sum_{k=1}^K N_k$. The partial derivatives we need for the Fisher information matrix are
$$
L_\rho(\theta,\rho|\underline{y}) = \prod_{k = 1}^K{N_k \choose y_k} \left[(1-\theta)(1-\psi)^{y^*} \psi^{N^*-y^*} + \theta\varphi^{y^*}(1-\varphi)^{N^*-y^*} -\prod_{k = 1}^K\phi_{s_k}(\theta)^{y_k}(1-\phi_{s_k}(\theta))^{N_k-y_k}\right]
$$

$$
L_\theta(\theta, \rho|\underline{y})
= \prod_{k = 1}^K{N_k \choose y_k} \left[-\rho(1-\psi)^{y^*} \psi^{N^*-y^*} + \rho\varphi^{y^*}(1-\varphi)^{N^*-y^*} + (1-\rho)\sum_{k=1}^K \frac{(1-\varphi-\psi)s_k(1-\theta)^{s_k-1}(N_k\phi_{s_k}(\theta) - y_k)}{\phi_{s_k}(\theta)(1-\phi_{s_k}(\theta))} \prod_{k = 1}^K\phi_{s_k}(\theta)^{y_k}(1-\phi_{s_k}(\theta))^{N_k-y_k}\right]. \\
$$
##### Beta distributed $\Theta_j$

We can model $\Theta_j \sim Beta(\alpha, \beta)$ where $\alpha = \theta(\rho^{-1} -1)$ and $\beta = (1-\theta)(\rho^{-1}-1)$. The density and derivatives with respect to $\theta$ and $\rho$ are
$$
\begin{align*}
f_\Theta(p)
&= \frac{p^{\alpha -1}(1-p)^{\beta -1}}{B(\alpha, \beta)} \\
&= \frac{p^{\theta(\rho^{-1} -1) -1}(1-p)^{(1-\theta)(\rho^{-1} -1) -1}}{B(\theta(\rho^{-1} -1), (1-\theta)(\rho^{-1} -1))} \\
\end{align*}
$$

$$
\begin{align*}
\frac{\part f_\Theta}{\part \theta}(p)
&= f_\Theta(p) (\rho^{-1} -1) \left[ log(p) - log(1-p) - \Digamma(\theta(\rho^{-1} -1)) + \Digamma((1-\theta)(\rho^{-1} -1)) \right] \\
&= f_\Theta(p) (\alpha + \beta) \left[ log(p) - log(1-p) - \Digamma(\alpha) + \Digamma(\beta) \right] \\
\end{align*}
$$

$$
\begin{align*}
\frac{\part f_\Theta}{\part \rho}(p)
&= f_\Theta(p) \rho^{-2} \left[\theta\left\{\Digamma(\theta(\rho^{-1} -1)) - log(p)\right\} + (1-\theta)\left\{\Digamma((1-\theta)(\rho^{-1} -1)) - log(1-p)\right\} - \Digamma((\rho^{-1} -1))\right] \\
&= f_\Theta(p) (\alpha + \beta +1)^2 \left[\frac{\alpha}{\alpha + \beta}\left\{\Digamma(\alpha) - log(p)\right\} + \frac{\beta}{\alpha + \beta}\left\{\Digamma(\beta) - log(1-p)\right\} - \Digamma(\alpha + \beta)\right] \\
\end{align*}
$$

where $\Digamma$ is the digamma function.

The resulting likelihood for $\theta$ and $\rho$ given pooled testing is
$$
L(\theta, \rho|y,N, s, \psi, \varphi) = \frac{{N \choose y}}{B(\alpha, \beta)} \int_0^1 p^{\alpha - 1} (1-p)^{\beta - 1} {\phi_s(p)}^y (1-\phi_s(p))^{N-y} dp
$$
where $B(.,.)$ is the beta function. In the case with $s=1$, perfect sensitivity and specificity $\phi_s(p)$ reduces to $p$ and the likelihood reduces to that of the beta binomial distribution:
$$
\begin{align*}
L(\theta, \rho|y, N, 1,1,1)
&= \frac{{N \choose y}}{B(\alpha, \beta)} \int_0^1 p^{\alpha - 1 + y} (1-p)^{\beta - 1 + N - y} dp \\
&= {N \choose y} \frac{B(\alpha + y, \beta + N - y)}{B(\alpha, \beta)}.
\end{align*}
$$

In the case with $s>1$, but perfect sensitivity and specificity $\phi_s(p)$ reduces to $1- (1-p)^s$ we can use a similar method to compute the integral:

$$
\begin{align*}
L(\theta, \rho|y,N,s,1,1)
&= \frac{{N \choose y}}{B(\alpha, \beta)} \int_0^1 p^{\alpha - 1} (1-p)^{\beta - 1 + s(N - y)} (1-(1-p)^s)^y dp \\
&= \frac{{N \choose y}}{B(\alpha, \beta)} \int_0^1 p^{\alpha - 1} (1-p)^{\beta - 1 + s(N - y)} \sum_{x = 0}^y {y \choose x}(-1)^x(1-p)^{sx} dp \\
&= \frac{{N \choose y}}{B(\alpha, \beta)} \sum_{x = 0}^y {y \choose x} \int_0^1 p^{\alpha - 1} (1-p)^{\beta - 1 + s(N - y)}  (-1)^x(1-p)^{sx} dp \\
&= \frac{{N \choose y}}{B(\alpha, \beta)} \sum_{x = 0}^y {y \choose x} (-1)^x B(\alpha, \beta + s(N - y + x)). \\
\end{align*}
$$

Computing derivatives we have
$$
L_\theta(\theta, \rho|y) = \frac{(\alpha + \beta){N \choose y}}{B(\alpha, \beta)}
\sum_{x = 0}^y (-1)^x {y \choose x} B(\alpha, s(N-y +x) + \beta)(\Digamma(\beta) - \Digamma(\beta + s(N - y + x)))
$$

$$
\begin{align*}
L_\rho(\theta, \rho|y) &=
\frac{(\alpha + \beta + 1)^2{N \choose y}}{B(\alpha, \beta)}
\sum_{x = 0}^y (-1)^x {y \choose x} B(\alpha, s(N-y +x) + \beta)\left[\frac{\beta}{\alpha + \beta}\left\{\Digamma(\beta) - \Digamma(\beta + s(N - y + x))\right\} + \Digamma(\alpha + \beta + s(N-y+x)) - \Digamma(\alpha + \beta) \right].
\end{align*}
$$

Care must be taken however in calculating the sums as the alternating sign can cause numerical stability issues for large $N$ and $s$. A similar method could be used to derive an expression for the general case where $s>1$ and sensitivity and specificity are imperfect, but it is too cumbersome for practical use and is likely to exacerbate the numerical stability issues. However, the integral can be solved numerically. Note however that for $\alpha < 1$ or $\beta < 1$, while the integrals still exist, the integrands diverge at $p=0$ and $p=1$ respectively. While standard quadrature methods can sometimes handle these divergences, in practice they can lead to computational challenges. These can be averted by decomposing the integrand into the sum of two components: one component without poles, and one that has poles at $p=0$ and/or $p=1$ but is a multiple of the density of a beta distribution and therefore has an integral that can be expressed in terms of beta functions. [I've coded this up for the package, and can include details in an appendix].

***Figure []*** shows computed of the Fisher information matrix $I(\theta, \rho)$ and the variance of the MLE for $\theta$, $\Sigma(\theta, \rho|N = 5, s = 10)_{1,1}$. This illustrates the convergence $\lim_{\rho \to 0}I(\theta, \rho)_{1,1} = I(\theta)$ for a range of $\theta$ . For $\theta \leq 0.1$  the variance of the MLE appears to decrease towards $\frac{1}{I(\theta)}$ as intra-cluster correlation $\rho$ approaches $0$. However, for $\theta \geq 0.3$ the variance of the MLE appears to increase as intra-cluster correlation approaches $0$, possibly without bound for $\theta \geq 0.4$. Moreover, for $\theta \leq 0.4$, the variance from a cluster survey is greater than a simple random survey for all calculated values of $\rho$. However for $\theta = 0.5$ the variance from the cluster survey is smaller than for a simple random survey for intra-cluster correlation greater than approximately 0.65.

One way to have an intuitive understanding of these surprising results comes from noting that for any given pool size $s$ prevalence can be sufficiently high such there is a very high probability that all pools will return a positive result. All positive pools provide very little information about either prevalence or intra-cluster correlation. However, if there is substantial intra-cluster correlation, this increases the probability that *some* clusters in a survey will return all positive pools, but reduces the probability that all clusters in a survey return all positive pools.

While this behaviour is counter-intuitive and suggest major problems with using cluster surveys with pooled testing to estimate prevalence, it should be noted that the most pathological behaviour is restricted to the case where prevalence is reasonably high — i.e. when one would generally not be considering using pooled testing. While it is possible to have this behaviour at low prevalence it only occurs with very large pools, far larger than what would be a sensible even for a non-cluster survey (see Figure []). Reducing the pool size avoids the pathological behaviour.

[Note also that you don't get this weird behaviour for the simple discrete distribution. Instead it behaves very nicely with $\lim_{\rho \to 1}\Sigma_{1,1}$ appearing to converge to $1/I(\theta)$, i.e. reduces to the simple random sampling case. But the discrete distribution is strange in its own way as will be seen later]

![](./Figures/FisherMatrixbeta.png)

**Figure []**. Fisher information matrix for $\theta$ and $\rho$ in a cluster sample design with pooled testing with ten units per pool ($s = 10$) and five pools per cluster ($N = 5$). Cluster-level prevalence ($\Theta_j$) is assumed to be beta distributed with mean $\theta = \frac{\alpha}{\alpha + \beta}$ and resultant intra-cluster correlation $\rho = \frac{1}{\alpha + \beta + 1}$. Dotted lines in upper left panel is Fisher information for $\theta$ for a simple random survey with the same pool size and total pools, $I(\theta)$. Dotted lines in the upper right panel is $\frac{1}{I(\theta)}$.  

![](./Figures/Variancebeta.png)

**Figure []** Upper left element of inverse Fisher information matrix for $\theta$ and $\rho$ in a cluster sample design with pooled testing with large pools ($50 \leq s\leq 2000$) and two pools per cluster ($N = 2$). Cluster-level prevalence ($\Theta_j$) is assumed to be beta distributed with mean $\theta = \frac{\alpha}{\alpha + \beta}$ and resultant intra-cluster correlation $\rho = \frac{1}{\alpha + \beta + 1}$. Dotted lines is inverse of Fisher information for $\theta$ for a simple random survey with the same pool size and total pools, $\frac{1}{I(\theta)}$.  

##### Logit-normal distributed $\Theta_j$

If we model $\Theta_j$ as logit-normal, i.e. $h(\Theta_j) \sim Normal(\mu,\sigma)$ where $h$ is the logit function, the likelihood is
$$
\begin{align*}
L(\theta, \rho, \psi, \varphi|y,N, s)
&= {N \choose y} \int_0^1 Normal(h(p)|\mu,\sigma)h'(p) {\phi_s(p)}^y (1-\phi_s(p))^{N-y} dp \\
&= {N \choose y} \int_0^1 \frac{Normal(h(p)|\mu,\sigma)}{p(1-p)} {\phi_s(p)}^y (1-\phi_s(p))^{N-y} dp. \\
\end{align*}
$$
There is no analytic expression for the mean and variance of a logit-normal distribution, therefore calculation of $\mu$ and $\sigma$ from $\theta$ and $\rho$ requires numerical inversion of the integral equations defining the moments.

It is difficult to compute the derivatives of $f_\Theta$ with respect to $\rho$ and $\theta$. Instead we compute $I(\mu,\sigma)$, the Fisher information with respect to canonical parameters for the logit-normal distribution, and use a change of variable i.e. $I(\theta, \rho) = A^TI(\mu,\sigma)A$, where $A$ is the Jacobian matrix for the transformation from $\{\theta, \rho\}$ to $\{\mu,\sigma\}$. $A$ is most conveniently calculated in terms of its inverse
$$
A =
\begin{bmatrix}
\frac{\part\theta}{\part \mu} & \frac{\part\theta}{\part \sigma}\\
\frac{\part\rho}{\part \mu} & \frac{\part\rho}{\part \sigma}\\
\end{bmatrix}^{-1}
$$
where
$$
\begin{align*}
\frac{\part\theta}{\part \mu} &= \int_0^1 Normal(h(p)|\mu,\sigma) \frac{h(p) - \mu}{\sigma^2} p dp\\

\frac{\part\theta}{\part \sigma} &= \int_0^1 Normal(h(p)|\mu,\sigma) \frac{(h(p) - \mu)^2 - \sigma^2}{\sigma^3} p dp\\

\frac{\part\rho}{\part \mu} &= \left[\int_0^1 Normal(h(p)|\mu,\sigma) \frac{h(p) - \mu}{\sigma^2} p^2 dp - \frac{\part\theta}{\part \mu}(2\theta + \rho(1-2\theta))\right]\left[\theta(1-\theta)\right]^{-1}\\

\frac{\part\rho}{\part \sigma} &= \left[\int_0^1 Normal(h(p)|\mu,\sigma) \frac{(h(p) - \mu)^2 - \sigma^2}{\sigma^3}  p^2 dp - \frac{\part\theta}{\part \sigma}(2\theta + \rho(1-2\theta))\right]\left[\theta(1-\theta)\right]^{-1}.
\end{align*}
$$

The logit-normal distribution arises naturally in random effect and mixed effect regression models as has been applied to pool-tested data in [references]. While the logit link function is the most common for analysis of binary response data, the above setup generalises readily to any appropriate link function, e.g. the complementary log log function [references] or the blended log-logit link functions proposed by Clark and Barr (2018).

The Fisher information matrix and its inverse display similar quantitative and qualitative results for logit-normal distributed $\Theta_j$ as for beta distributed $\Theta_j$ (compare below with figure from previous section).

![](./Figures/FisherMatrixlogitnorm.png)

**Figure []**. Fisher information matrix for $\theta$ and $\rho$ in a cluster sample design with pooled testing with ten units per pool ($s = 10$) and five pools per cluster ($N = 5$). Cluster-level prevalence ($\Theta_j$) is assumed to be logit-normally distributed with mean $\theta$ and resultant intra-cluster correlation $\rho$. Dotted lines in upper left panel is Fisher information for $\theta$ for a simple random survey with the same pool size and total pools, $I(\theta)$. Dotted lines in the upper right panel is $\frac{1}{I(\theta)}$.  

#### Optimising cluster surveys with pooled testing

For a cluster survey using a test with known sensitivity and specificity, there are three design parameters to be chosen: number of locations ($J$), number of pools/groups per location ($\underline{N}$), and units per pool/group ($\underline{s}$). As before, let $c_u$ be the marginal cost of collecting one unit and $c_p$ the marginal cost of one pool/group. Additionally consider the marginal cost of collecting vectors at an additional location $c_l$. The total costs associated with the collection of vectors, handling and testing of groups from one sampling location is then $J(c_l + N_{tot} c_p + \underline{N}\cdot \underline{s} \ c_u)$, where $N_{tot} = \sum_k N_k$. For large $J$, the variance of the maximum likelihood estimator of the prevalence is approximately $(J \ \Sigma_{1,1})^{-1}$.  In analogy to the the pool-tested simple random survey where we try to minimise the marginal cost of Fisher information, in a pool-tested cluster survey we wish to minimise 

$$
\frac{c_l + N_{tot} c_p + \underline{N}\cdot \underline{s} \ c_v}{\Sigma(\theta, \rho, \underline{N},\underline{s},\psi,\varphi)_{1,1}}.
$$

Though this is no longer strictly the unit cost of Fisher information we call it the *cost of information* to give both measures a unified name. Note however, that it is not always straightforward to make comparisons between the unit information from a simple random survey and the cluster. First, cluster sampling introduces a new cost $c_l$, which did not have to be considered when optimising a simple random survey. Second, the cost of sampling a unit, $c_u$, will usually be different between the the two survey methods. In  fact, cluster surveys are often motivated by reducing $c_u$ or because simple random sampling is so difficult that $c_u$ is very large and unquantified. However, to get an understanding for how unknown non-zero correlation impacts the optimal design of experiments, we consider as a comparator the case where units are collected and test in a pooled-tested cluster framework (with associated cost structure) but that prevalence is estimated from the results assuming that there is no correlation between units at location. In this case the unit cost of information is
$$
\frac{c_l + N_{tot} c_p + \underline{N}\cdot \underline{s} \ c_v}{I(\theta| \underline{N},\underline{s},\psi, \varphi)}.
$$
The figures below illustrates how correlation, prevalence, and relative cost of survey components determine the cost of information and optimal pooling configurations. To reduce the search space for pooling designs we consider designs with $N$ pools of a size $s$ where the above reduces to
$$
\frac{c_l + N c_p + N \cdot s \ c_v}{I(\theta| N,s,\psi, \varphi)} = \frac{c_l/N + c_p + s \ c_v}{I(\theta| 1,s,\psi, \varphi)}.
$$
[I haven't tried to show that any of these have unique $N$ or $s$ that optimise them. It was hard enough in the simple random sampling case! My intuition and numerical experiments suggests its all should be very well behaved, but I could be wrong… My algorithm for identifying optimal s for a given N is to use stats::optimise to identify real value of s that minimises cost of information and then identify which of floor(s) and ceiling(s) is better. To identify optimal s and N simultaneously I start with N = 1, identify optimal s for N=1, then ++N until the unit cost of information for the optimal s stops improving. It's pretty basic but is *seems* to work fine (and is very fast). This assumes unimodality over the $s \times N$ search space which I'm guessing is true].

The first figure compares the cost of information from a cluster survey where there is assumed to be no correlation, to those where correlation is unknown but very small ($\rho = 0.01$) or moderate ($\rho = 0.1)$ and to a simple random survey with pooled testing. In all scenarios costs are $c_l = 40, c_p = 4, c_u = 1$. Individual testing ($s = 1$) is not optimal for the range of scenarios shown below (but would be for large $\theta$), however the optimal value of $s$ is highly dependent on the other factors. When we can assume no correlation, the unit information curve with respect to $s$ has a very broad trough, i.e. a wide range of pool sizes around the optimal value results in similar cost of information. However, accounting for clustering increases the cost of information. While this increase is observed for all $s$, it is marginal for small pools and very marked for large pools. Consequently the unit information curve with respect to $s$ develops a narrower trough and the optimal pool size decreases. 

However, the shape of the cost of unit information function depends on the form of the distribution of prevalence at locations. There is very little difference between the logit-normal and beta distributed cases, however the simple discrete distribution considered in previous sections leads to markedly different results (second figure in this section). There is little to no difference for small pools ($s \approx 1$), however there is a substantial difference for larger pools. Most notably, the unit information function no longer has a single local minima in $s$ for a given $N$, and the minimum of the range of pool sizes considered ($s = 1-50$) is larger than for either simple random sampling or cluster sampling with known 0 correlation. Consequently, the optimal pool size is larger for a cluster survey if assuming the discrete distribution, than for simple random sampling. 

The third figure illustrates that increasing correlation, prevalence, and number of pools per location all favours smaller pools for beta-distributed site prevalence. For example, with the costs assumed above,  2.5% prevalence, and four pools per location, the optimal pool size is 25 for a cluster survey assuming no correlation ($\rho = 0$), 15 for a simple random survey but only 4 if there is even only moderate (unknown) correlation ($\rho = 0.1$). Moreover, the range of pool sizes that are approximately cost-optimal is also greatly reduced when there is correlation. For the same scenario the range of pool sizes that cost no more than 10% more than the optimal pool is 15-42 for a cluster survey assuming no correlation, 8-29 for a simple random survey, and only 3-7 for moderate unknown correlation.

However in general may we wish to simultaneously optimise both the size ($s$) of number of pools per site ($N$). For our comparison case of cluster sampling without correlation at sites, cost of information is maximised as $N$ goes to infinity at which point the cost information and optimum pool size converges to that of obtained for simple random sampling. The results of such an optimisation for a range of cost structures are given in the fourth figure below, again assuming beta-distributed site prevalences. As before, increasing prevalence or correlation decreases the optimal pool size. However, increasing prevalence or decreasing correlation increases optimal number of pools per site. However, the optimal number of  total units to be sampled at each site ($sN$) is less sensitive to prevalence or correlation for the range of scenarios considered. The cost per pool is the primary determining factor for the optimal pool size, with higher costs favouring larger pools. On the other hand, cost per location is the primary determining factor for the optimal number of pools per location. However increasing cost per pool does also have small effect to decrease the optimal number of pools per site also.

The final figure in this section shows the design effect as a function of correlation and number of units per site, pool size, and prevalence. Increasing prevalence increases the design effect unless individual testing used in which case there is no effect (this is for beta case — much more strange for discrete case ). Increasing the number of units caught at a site for a given pool size increases the design effect for moderately correlation ($\rho = 0.1$), has no effect when there is no correlation (and this is known), but has a more complicated relationship for very low correlation ($\rho = 0.01$).


![](./Figures/unit info cluster alt beta.png)

***Figure []*** Unit cost of information for prevalence estimation surveys with pooled testing and cluster or simple random surveys. Survey designs considered are simple random sampling ($\rho =$ NA) and cluster surveys with known zero clustering at site level ($\rho = 0$) or unknown degree of clustering ($\rho = 0.01, 0.1$). Tests are assumed to be perfectly specific and sensitive. The component survey costs are costs per cluster/location ($c_l = 40$; cluster sample designs only), costs per pool/test $c_p = 4$, and costs per unit collected $c_u = 1$. Fisher information calculations for cluster surveys with non-zero clustering assume prevalence is beta-distributed across sampling locations.

![](./Figures/unit info cluster alt discrete.png)

***Figure []*** Same as previous figure but Fisher information calculations for cluster surveys with non-zero clustering assume prevalence follows a simple discrete distribution across sampling locations.

![](./Figures/optimal cluster.png)

***Figure []*** Pool sizes that minimises cost of information in prevalence estimation surveys (solid line) and the range of pool sizes that are no more than 10% more costly than the optimal (shaded area bounded by dotted lines). Survey designs considered are simple random sampling ($\rho =$ NA) and cluster surveys with known zero clustering at site level ($\rho = 0$) or unknown degree of clustering ($\rho = 0.01, 0.1$). Tests are assumed to be perfectly specific and sensitive. The component survey costs are the same as previous figure.  Fisher information calculations for cluster surveys with non-zero clustering assume prevalence is beta-distributed across sampling locations. Note: the range of near-optimal pool sizes are only given for $\rho = 0$ and $\rho = 0.1$ to avoid over-plotting.

![](./Figures/optimal cluster sN.png)

***Figure []*** The pool sizes ($s$), number of pools per sampling site ($N$), and units per location ($N \cdot s$) that minimises cost of information in prevalence estimation surveys. Survey designs considered are simple random sampling ($\rho =$ NA) and cluster surveys with known zero clustering at site level ($\rho = 0$) or unknown degree of clustering ($\rho = 0.01, 0.1$). Tests are assumed to be perfectly specific and sensitive. The cost per pool ($c_p$) and cost per location ($c_l$) vary across columns. Fisher information calculations for cluster surveys with non-zero clustering assume prevalence is beta-distributed across sampling locations. Note: optimal $N$ is infinite for for $\rho = 0$ at which point the cost information and optimal $s$ reduce to that of simple random survey.

 ![](./Figures/FigDesignEffectCorrbeta.png)

***Figure []*** Design effect as a function of prevalence, within-site correlation, pool size and number of units per site. The dotted horizontal line indicates the reference design effect: simple random sampling with individual testing. All tests are assumed to be perfectly sensitive and specific. In each column the total units per site ($sN$) is fixed so $N$ is inversely proportional to $s$.  Variation of prevalence across sampling locations is modelled with the beta distribution. Note the different y-axis scales in each row. 

### Power size calculations based on Fisher information

There's a couple ways to do this

* Use asymptomatic normality of MLE and the Fisher information to build confidence intervals and procedure for accepting rejecting the null hypothesis. 
* Calculate a design effect from Fisher information and multiply sample size from some existing tool (e.g. binomial sample size calculator).

There's also a couple of different set ups where you want to calculate sample sizes:

* Compare prevalence estimated from one sample (cluster or otherwise) to a threshold
* Compare prevalence estimates for different samples (cluster or otherwise) where the samples are independent but might have different pooling designs

Other complications:

* One and two sided tests
* Random catch sizes (see next section)

#### Comparing prevalence to threshold

Consider the case where there is a threshold value of prevalence, $\theta_0$, and we wish to use a cluster pool survey to demonstrate that the prevalence is above or below this threshold. Using asymptotic normality of the MLE we can calculate minimum sample sizes required for target power, $p$, confidence level, $c$, and a design prevalence, $\theta_a$. We assume that the pooling design will consist of some number, $J$, of replicates of a pooling design. In a non-cluster survey this we could summarise this with pooling design $J\underline{N}, \underline{s}$. For a cluster survey this might require $J$ locations each with pooling design $\underline{N},\underline{s}$. In either case the total number of units sampled is $n = J \underline{N} \cdot \underline{s}$​​. The standard error at the design prevalence, $\sigma_a$, and standard error at the threshold prevalence, $\sigma_0$, are
$$
\sigma_x = \frac{s_x}{\sqrt{J}} = 
\begin{cases}
\left(J \ I(\theta_x | \underline{N}, \underline{s} \ \right)^{-\frac{1}{2}} & \text{for simple random surveys,}\\
\left(J^{-1} \ {I(\theta_x,\rho|\underline{N},\underline{s})^{-1}}_{1,1}\right)^{\frac{1}{2}} & \text{for cluster surveys.}
\end{cases}
$$
The power to reject one-sided and two-sided null hypotheses can be approximated as
$$
p \approx
\begin{cases}
z^{-1}\left(\frac{\left(\theta_0-\theta_a\right) - z(c) \sigma_0}{\sigma_a} \right) & H_0: \theta_a > \theta_0 \\
z^{-1}\left(\frac{\left(\theta_a-\theta_0\right) - z(c) \sigma_0}{\sigma_a} \right) & H_0: \theta_a < \theta_0 \\
z^{-1}\left(\frac{\left(\theta_a-\theta_0\right) - z(\frac{c+1}{2}) \sigma_0}{\sigma_a} \right) + z^{-1}\left(\frac{\left(\theta_0-\theta_a\right) - z(\frac{c+1}{2}) \sigma_0}{\sigma_a} \right) & H_0: \theta_a = \theta_0 \\
\end{cases}
$$


where $z$ is the quantile function of the standard normal distribution. For one-sided hypothesis tests, if the power is to exceed $p$, $J$ must be at least
$$
J \geq \left[\frac{z(c) s_0 + z(p) s_a}{\theta_0 - \theta_a}\right]^2.
$$
#### Comparing prevalence in two samples

Now consider the case where we conduct pooled cluster surveys in two populations; for instance, surveys in two separate regions, or surveys at two time points before and after an intervention in a single region. Here we assume that the samples are statistically independent. There are now two design prevalences, $\theta_a$ and $\theta_b$, and two pooling designs consisting of $J_a$ and $J_b$ replicates of represented by $\underline{N}_a,\underline{s}_a$ and $\underline{N}_b,\underline{s}_b$. In the cluster sampling case there may also be two correlations $\rho_a$ and $\rho_b$. In both the simple random and cluster sampling cases, the total number of units across both samples is $n = n_a + n_b = J_a \underline{N}_a \cdot \underline{s}_a + J_b \underline{N}_b \cdot \underline{s}_b$. Using the independence of the samples and the asymptotic normality of each MLE of prevalence, the difference in the two MLEs is approximately normal with standard error, $\sigma_\delta$
$$
\sigma_\delta  = \sqrt{\sigma_a^2 + \sigma_b^2}  = \sqrt{\frac{s_a^2}{n_a} + \frac{s_b^2}{n_b}},
$$
where $s_x$ and $\sigma_x$​ are defined similarly to the previous section, with the generalisation to different pooling designs:
$$
\sigma_x = \frac{s_x}{\sqrt{J_x}} = 
\begin{cases}
\left(J_x \ I(\theta_x | \underline{N}_x, \underline{s}_x \ \right)^{-\frac{1}{2}} & \text{for simple random surveys,}\\
\left(J_x^{-1} \ {I(\theta_x,\rho_x|\underline{N}_x,\underline{s}_x)^{-1}}_{1,1}\right)^{\frac{1}{2}} & \text{for cluster surveys.}
\end{cases}
$$
Under the null hypothesis $\theta_a = \theta_b = \theta_c$, the standard error of the difference is $s_c \sqrt{(\frac{1}{J_a} + \frac{1}{J_b})}$. We reject the null hypothesis with confidence $c$, if $\left| \widehat{\theta_a} - \widehat{\theta_b} \right| > z(\frac{1+c}{2}) \widehat{s_c} \sqrt{(\frac{1}{J_a} + \frac{1}{J_b})}$. 


















and the standard error of the estimate of prevalence yielded from combining the two samples is 
$$
\sigma_x = \frac{s_x}{\sqrt{J}} = 
\begin{cases}
\left(J \ I(\theta_x | \underline{N}, \underline{s} \ \right)^{-\frac{1}{2}} & \text{for simple random surveys,}\\
\left(J^{-1} \ {I(\theta_x,\rho|\underline{N},\underline{s})^{-1}}_{1,1}\right)^{\frac{1}{2}} & \text{for cluster surveys.}
\end{cases}
$$
The power to reject one-sided and two-sided null hypotheses can be approximated as
$$
p \approx
\begin{cases}
z^{-1}\left(\frac{\left(\theta_b-\theta_a\right)}{\sqrt{\sigma_a^2 + \sigma_b^2}} - z(c) \right) & H_0: \theta_a > \theta_b \\
z^{-1}\left(\frac{\left(\theta_a-\theta_b\right)}{\sqrt{\sigma_a^2 + \sigma_b^2}} - z(c) \right) & H_0: \theta_a < \theta_b \\
z^{-1}\left(\frac{\left(\theta_a-\theta_b\right)}{\sqrt{\sigma_a^2 + \sigma_b^2}} - z(\frac{c+1}{2})\right) + z^{-1}\left(\frac{\left(\theta_b-\theta_a\right)}{\sqrt{\sigma_a^2 + \sigma_b^2}} - z(\frac{c+1}{2}) \right) & H_0: \theta_a = \theta_b \\
\end{cases}
$$

The assumption of independent samples could be violated if, for instance, the two samples were collected from the same set of clusters and the cluster-level prevalences for the two populations of interested are correlated. This situation arises naturally in molecular xenomonitoring surveys, where mosquitos are collected in traps from a number of (randomly) selected sites (clusters) but traps catch more than species of mosquito. If these mosquitos are sorted and then pooled and tested for the presence of a pathogen separately, this yields a pooled cluster sample for each mosquito species. However, one may reasonably expect (and at least one study has demonstrated) that the prevalence of the marker in the different species can be correlated at the site level. Though beyond the scope of this paper, positive  correlations between species at the cluster level will likely reduce the statistical power of surveys to detect differences at the population level.

***Should I elaborate two-sided tests?***

Alternatively we can use existing sample size calculation methods for binary tested data and then multiply by the design effect to appropriately adjust the total sample size. A common method of calculating power for when using simple random sampling and individual testing is the arcsine transformation methods (eg.[reference]). In this method power is approximated as:
$$
z(p) \approx 2\left(\arcsin(\sqrt{\theta_0})-\arcsin(\sqrt{\theta_a})\right) \sqrt{n} - z(c)
$$
and for given level of power the minimum sample size is
$$
n \geq \left[\frac{z(c) + z(p)}{2\left(\arcsin(\sqrt{\theta_0})-\arcsin(\sqrt{\theta_a})\right)}\right]^2.
$$
To estimate the total number of units required for a survey with units tested in pools and either simple random or cluster sampling, we simply multiply by the above sample size by design effect. In the cluster sampling case the minimum of replicates is then calculated as
$$ {\}
\begin{align*}
J  &\geq  \frac{D(\theta_0, \rho, \underline{N},\underline{s},\psi, \varphi)}{\underline{N}\cdot\underline{s}} \left[\frac{z(c) + z(p)}{2\left(\arcsin(\sqrt{\theta_0})-\arcsin(\sqrt{\theta_a})\right)}\right]^2 \\
& \equiv I(\theta_0|1,1,\varphi,\psi)  \ \Sigma(\theta_0, \rho|\underline{N}, \underline{s}, \varphi, \psi)_{1,1} \left[\frac{z(c) + z(p)}{2\left(\arcsin(\sqrt{\theta_0})-\arcsin(\sqrt{\theta_a})\right)}\right]^2.
\end{align*}
$$



The figures below compare the two approaches for a range of parameters. Generally the agreement is good when the target and design prevalence are both low which is the typical case for pool-tested surveys, however the computed minimum sample size is typically larger when using the design 

![](.\Figures\SampleSizebeta.png)

***Figure[]*** Minimum sample size to achieve 80% power and 95% that true prevalence is greater/less than 1%, using a one-sided test and pooled testing. The top row is the case for simple random sampling (or cluster sampling where correlation at sites is known to be 0). Other rows are for cluster sampling with very low and modest correlation but where correlation needs to be estimated. Different line types indicate the size of pools, while each columns indicates the number of units collected per site. The colours indicate how sample sizes are calculates. Blue uses the arcsine method and assumes simple random sampling and individual tests. Red multiplies the arcsine method by the design effect. Green uses the Fisher information matrix directly to determine sample sizes. Correlation model assumed beta distributed site prevalence.

![](.\Figures\SampleSizediscrete.png)

***Figure[]*** Minimum sample size to achieve 80% power and 95% that true prevalence is greater/less than 1%, using a one-sided test and pooled testing. Same as previous figure but correlation model assumes the simple discrete distribution of site prevalence.

[Rather than optimising the unit cost of information as we have in the previous sections, we could also choose a design which meets minimum power and confidence requirements while minimising cost direct. This turns out to have a slightly different optimum designs in numerical experiments, but I haven't played around with this too much]

### Surveys with random catch sizes

Surveys using pooled testing do not always use fixed sample sizes. In molecular xenomonitoring surveys, it is common instead to place traps to catch disease vectors for set period of time and analyse however many vectors were caught at the end of the time period. Consequently in these surveys neither $\underline{N}$ nor $\underline{s}$ are predetermined but random. They may often be some predetermined rules about how the units caught are split into pools. Example rules might be to divide the units, however many are caught, into the pools of size 25, with any remainder assigned to an additional, smaller pool. Another example might be to divide the units into four pools of approximately equal size. To provide a notation for possible rules, let the number of units sampled at each site $j$, $V_j$, be iid random variables with probability mass function $f_V(v|\underline{\xi})$ defined by (known) parameters $\underline{\xi}$. Let $r_N$ and $r_s$ be functions $\mathbb{N} \to \mathbb{N}^\infty$ describing the pooling rules such that $\underline{N}_j = r_N(V_j)$ and $\underline{s}_j = r_s(V_j)$. Noting that each unit is placed in at most one group but some rules may choose to discard some units (say if the catch size is much larger that expected) the pooling rules must satisfy $r_s(v)\cdot r_N(v) \leq v$. We assume the same pooled rules apply for each sampling location and that number and size of the pools are deterministic given $v$; however we assume that units are randomly assigned to pools to make up the required size. We also assume catch sizes at a location are independent of the testing results of individual units.

Under this setup the Fisher information from a site is the Fisher information averaged over possible catch sizes. [I probably don't have to write the below up since it seems very general, obvious, and not specific to pooled testing in any way, but it might be helpful to spell it out?]

First note that the likelihood of the parameters of interest given a catch size $v$ and testing results $\underline{y}$ (and pooling rules and other known parameters) is simply
$$
L(\theta, \rho|\underline{y},v;r_N,r_s,\underline{\xi},\psi, \varphi) = f_V(v|\underline{\xi}) L(\theta, \rho|\underline{y},r_N(v),r_s(v); \psi, \varphi).
$$
Consequently the derivatives the likelihood are also simply
$$
L_\theta(\theta, \rho|\underline{y},v) = f_V(v|\underline{\xi}) L_\theta(\theta, \rho|\underline{y})\\
L_\rho(\theta, \rho|\underline{y},v) = f_V(v|\underline{\xi}) L_\rho(\theta, \rho|\underline{y})\\
$$
where have suppressed the dependence on pooling rules, catch distribution parameters, and test sensitivity and specificity for notational compactness. Recall that $\Xi_{\underline{N}}$ is the set of all possible outcomes $\underline{y}$ given pool numbers $\underline{N}$. Then the Fisher information from a cluster survey with pooled testing and random catch sizes is
$$
\begin{align*}
I(\theta, \rho|\xi,r_N,r_s,\psi, \varphi)
&=\sum_{v = 0}^{\infty}\sum_{\underline{y} \in \Xi_{r_N(v)}}
L(\theta, \rho|\underline{y},v)^{-1}
\begin{bmatrix}
 L_\theta(\theta, \rho|\underline{y},v) ^2 &
 L_\theta(\theta, \rho|\underline{y},v) L_\rho(\theta, \rho|\underline{y},v)\\
 L_\theta(\theta, \rho|\underline{y},v) L_\rho(\theta, \rho|\underline{y},v) &
 L_\rho(\theta, \rho|\underline{y},v) ^2\\
\end{bmatrix}\\
&=\sum_{v = 0}^{\infty}\sum_{\underline{y} \in \Xi_{r_N(v)}}
 f_V(v|\xi) L(\theta, \rho|\underline{y})^{-1}
\begin{bmatrix}
 L_\theta(\theta, \rho|\underline{y}) ^2 &
 L_\theta(\theta, \rho|\underline{y}) L_\rho(\theta, \rho|\underline{y})\\
 L_\theta(\theta, \rho|\underline{y}) L_\rho(\theta, \rho|\underline{y}) &
 L_\rho(\theta, \rho|\underline{y}) ^2\\
\end{bmatrix} \\
&= \sum_{v = 0}^{\infty} f_V(v|\xi)
  I(\theta, \rho|r_N(v), r_s(v),\psi, \varphi)\\
&= E_V[I(\theta, \rho|r_N(V), r_s(V),\psi, \varphi)].
 \end{align*}
$$

[I have code which is set up to do these calculations for negative binomially distributed catch sizes. The calculations aren't particularly hard numerically since you just do finite sums until the sum converges. I'm not we would really want to go into it any more than this here. Early experiments suggest that if Fisher information calculations based on the *expected* catch tend to overestimate Fisher information based on random catches in the cluster randomised setup, i.e. $E[I(\theta, \rho| r_N(V),r_s(V))] \leq I(\theta, \rho|r_N(E[V]), r_s(E[V]))$. One could probably prove this (or give conditions when this happens) if one could show that (or when) Fisher information increases sub-linearly with catch size. This should be an equality in the case with individual testing and simple random sampling]

### Cluster surveys — multiple levels of clustering

Many MX surveys designs adopt designs with multiple levels of clustering. For instance sampling my proceed by first selecting random areas, then selecting random trapping sites from within these areas. One might analyse the results of such surveys with spatially correlated random effects for each sampling site. However, in analogy to the results in previous sections, we might also use multiple levels of nested random effects.

Let $V_{ijk}$ be the state (0 or 1) of unit $i$ selected from sampling site $j$ within sampling area $k$. In analogy to before we say that two units from different areas and sites are independent, however there is some degree of correlation between units from the same site and another degree of correlation from different sites but the same area:

$$
Corr(V_{ijk},V_{i'j'k'}) = \begin{cases}
1, & i = i', j = j', k = k' \\
\rho_1, & i \neq i', j = j', k = k' \\
\rho_2, & i \neq i', j \neq j', k = k' \\
0, &  k \neq k' \\
\end{cases}
$$

As before there are many different correlation structures which have this property, but as before we consider a general approach to specify the joint distribution of $\{V_{ijk}\}$ where we model the prevalence as variable across sites and areas, with units being independent conditioned on the site.  Let $\Theta_{jk}$ represent the prevalence at site $j$ in area $k$ 
$$
V_{ijk}|\Theta_{jk} \sim Bern(\Theta_{jk}),
$$
where $\{V_{ijk}|\Theta_{jk}\}$ are independent.

We model the prevalence in each area $k$ $\{\Theta_{k}\}$ as i.i.d. random variables with mean $\theta$ and support on $[0,1]$ and density $f_\Theta (p)$. We model the prevalence in each site in each location $\{\Theta_{jk}\}$ as identically distributed random variables with mean $\theta$ and conditional mean $E[\Theta_{jk}|\Theta_k] = \Theta_k$, support on $[0,1]$, and conditional densities $g(p|\Theta_k)$. We assume that prevalence at each site $\Theta_{jk}$ is independent from prevalence at sites in distinct areas, and assume that two site prevalences from the same area are independent when conditioned on area level prevalence.

Then, if $n_{jk}$ is the number of units from site $j$ in area $k$ then
$$
\sum_{i=1}^{n_{jk}} V_{ijk} |\Theta_{jk} \sim Bin(\Theta_{jk},n_{jk}).
$$
This set up preserves $E[V_{ijk}] = \theta$. Units from different sites and areas are independent and therefore uncorrelated. Correlation between two distinct units from the same site (and therefore same area), $\rho_1$, is
$$
\begin{align*}
\rho_1 := Corr(V_{ijk},V_{i'jk}) &= \frac{E[V_{ijk}V_{i'jk}] - E[V_{ijk}][V_{i'jk}]}{\sqrt{Var[V_{ijk}]Var[V_{i'jk}]}} \\
&=  \frac{P(V_{ijk} = V_{i'jk} = 1) - E[V_{ijk}]^2}{Var[V_{ijk}]} \\
&=  \frac{E[\Theta_{jk}^2] - E[\Theta_{jk}]^2}{Var[V_{ijk}]} \\
&= \frac{Var[\Theta_{jk}]}{\theta (1-\theta)}.
\end{align*}
$$

Correlation between two distinct units from the same area $k$ but distinct sites $j$ and $j'$, $\rho_2$, is
$$
\begin{align*}
\rho_2 := Corr(V_{ijk},V_{i'j'k}) &= \frac{E[V_{ijk}V_{i'j'k}] - E[V_{ijk}][V_{i'j'k}]}{\sqrt{Var[V_{ijk}]Var[V_{i'j'k}]}} \\
&=  \frac{P(V_{ijk} = V_{i'j'k} = 1) - E[V_{ijk}]^2}{Var[V_{ijk}]} \\
&=  \frac{E[\Theta_{k}^2] - E[\Theta_{k}]^2}{Var[V_{ijk}]} \\
&= \frac{Var[\Theta_{k}]}{\theta (1-\theta)}.
\end{align*}
$$

$$
\begin{align*}
P(V_{ijk} = V_{i'j'k} = 1)

& = \int_0^1 f_\Theta(p) E[\Theta_{jk}|\Theta_k = p] E[\Theta_{j'k}|\Theta_k = p] dp \\
& = \int_0^1 f_\Theta(p) E[\Theta_{jk}|\Theta_k = p]^2 dp \\
& = \int_0^1 f_\Theta(p) p^2 dp \\
& = E[\Theta_k^2]
\end{align*}
$$

Note the while the above setup has the same degree of correlation $\rho_1$ between units from a single site for every pair of sites across all areas, if you condition on area-level prevalence, $\Theta_k$ then intra-site correlations become zero, and inter-site correlations depend on $\Theta_k$
$$
Corr(V_{ijk}, V_{i'j'k})|\Theta_k =
\begin{cases}
1, & i = i', j = j'\\
\frac{Var[\Theta_{jk}|\Theta_k]}{\Theta_k(1-\Theta_k)}, & i \neq i', j = j' \\
0, &j \neq j'\\

\end{cases}
$$
We can use any two-parameter families with support on $[0,1]$ to model the distribution of prevalences across sites. The most common way to analyse data from nested cluster surveys like this would be to use a generalised linear model with two levels of random effects which are normally distributed on the link scale. If $\sigma_1$ and $\sigma_2$ are the standard deviations of the site-level and area-level random effects, and $h$ is the link function, then in this framework

$h^{-1}(\Theta_{jk}) \sim Normal(\mu, \sigma_1^2 + \sigma_2^2)$

However the distribution of $\Theta_k$ is more complicated. It does **not** follow

$h^{-1}(\Theta_k) \sim Normal(\mu, \sigma_2^2)$

as one might expect. To see this note that this would **not** preserve $E[\Theta_k] = E[\Theta_{jk}]$.

Instead if have $W_k \sim Normal(0, \sigma_2^2)$ i.i.d. and $W_{jk} \sim Normal(0, \sigma^2_1)$ i.i.d. then

$$
\begin{align*}
\Theta_{jk}
& \sim h(\mu + W_k + W_{jk})\\
\Theta_k
& \sim  E[h(\mu + W_k + W_{jk})|W_k]\\
&   =   \int_{-\infty}^\infty h(w) f_N(w|\mu + W_k, \sigma_1^2) dw
\end{align*}
$$

This can be generalised to three or more levels of clustering.

### Alternative parameterisations for Fisher information calculations

The complicated form of the correlation structures presented in the previous section is unfortunate and perhaps could be avoided by parameterising our distributions (and calculating Fisher information) with an alternate parameterisation.

For this class of random effect models considered in the previous section, the commonly used STATA software calculates an ICC of the following form
$$
ICC_1 = \frac{\sigma_1^2 + \sigma_2^2}{\gamma + \sigma_1^2 + \sigma_2^2}
$$

$$
ICC_2 = \frac{\sigma_2^2}{\gamma + \sigma_1^2 + \sigma_2^2}
$$

where $\gamma$ is the variance of standard logistic random variable (i.e. $\pi^2/3$). (See Goldstein et al 2002)

The advantage of this version of the ICC, in addition to ease of computation, is independence of $\theta$. Therefore in a mixed effect model, the ICCs are the same for every combination of fixed effects. The major downside is that while their names would suggest that ICCs should approximate $\rho_1$ and $\rho_2$, whenever prevalence is close to 0 (or 1) the ICCs are potentially orders of magnitude larger than the latter. Nevertheless, if the 'correlation coefficients' reported in the literature are going to be of the type calculated by STATA, and these are going to be used as a starting point for sample size calculations, then it is worth considering parameterising our Fisher information calculations in terms of $\theta$ and $ICC$.

This can be done fairly straightforwardly as we have in the previous sections with a change of variables.

## Designing surveys for the detection of disease markers

[Note that I actually wrote out this section first, but it is — I think — less interesting. I have used slightly different notation in this section that in the previous and rely more heavily on the language and examples of molecular xenomonitoring rather than using more generic language. I will correct once I know the right direction/format for this section.]

We consider the case where the aim the survey is to correctly identify whether locations do or do not have units that are positive for the disease marker. There are four outcomes of a survey

* The marker is absent in the population ($\theta = 0$) and no pools/groups return a positive result (true absence)
* the marker is present in the population ($\theta >0$)  and at least one pool/group returns a positive result (true presence)
* the marker is absent in the population ($\theta =0$) and at least one pool/group returns a positive result (false presence; type I error)
* the marker is present in the population ($\theta >0$)  and no pools/groups return a positive result (false absence; type II error)

We want to identify survey designs that keep the probability of false absence (type II error) and false presence (type I error) below acceptable limits while minimising cost. To calculate type I and type II errors we will need to consider the characteristics of the test and the survey design.

### Type II error using a simple random survey with a fixed sample size

A type II error can happen in two ways. Though the marker is present in the population, by chance the sampled units may all be negative. Alternatively, the sample may include one or more units carrying the marker of interest, but an imperfectly sensitive test might return a negative result for the pools/groups containing the positive units. Note that an imperfectly specific test acts to *decrease* type I error. A pool/group without the marker may return a positive result correctly indicating the *presence* of the marker in the *population*, though incorrectly identifying the presence of the marker in *group/pool*.

To calculate type I error for detection of a maker from a simple random sample, we reiterate two assumptions about the test characteristic:

* for pools/groups with at least one vector with the marker, the test sensitivity ($\varphi$) does not depend on the number vectors with the marker or the total size of the pool/group
* false negative test and false positive test results occur at random and independently across pools/groups

With these assumptions we can show that the probability that a single pool/group of size $s$ returns a positive results is:

$(1-\varphi - \psi)(1-\theta)^s + \varphi$

and type II error rate for a simple random survey with $N$ pools/groups with $s$ vectors per group is

$\epsilon_1(\theta, s, N, \varphi, \psi) = (1 - (1-\varphi - \psi)(1-\theta)^s - \varphi)^N$.

The type II error for simple random surveys with a mix of $K$ different pool sizes is similar. If the survey involves $\underline{N}$ pools of sizes $\underline{s}$ for a total of $N_{tot} = \sum_{k=1}^KN_k$ pools and $n = \underline{N}\cdot\underline{s}$ units, then the type II error rate is

$\epsilon_1(\theta, \underline{s}, \underline{N}, \varphi, \psi) = \prod_{k=1}^K(1 - (1-\varphi - \psi)(1-\theta)^{s_k} - \varphi)^{N_k}$.

Note that if the test is perfectly sensitive ($\varphi = 1$) then this reduces to

$\epsilon_1(\theta, \underline{s}, \underline{N}, 1, \psi) = \prod_{k =1}^K\psi^{N_k}(1-\theta)^{s_kN_k} = \psi^{N_{tot}}(1-\theta)^{n}$

and furthermore if the test is also perfectly specific this further reduces to 

$\epsilon_1(\theta, \underline{s}, \underline{N}, 1,1) = (1-\theta)^{n}$.

In other words, for a given number of units to be tested with a perfect test, the type-II error does not depend on whether units are tested in pools/groups or individually. Similarly if a perfectly sensitive test used, the type II error depends only on the total number units and pools/groups and not how the units are distributed amongst pools/groups. 

### Type II error using a simple random survey with a random sample size

In the previous section we calculated the type II errors assuming the number of groups ($\underline{N}$) and the number of vectors per group ($\underline{s}$) were fixed by design. This might happen if the survey design involves continuing sampling until a fixed number of units, $n = \underline{N}\cdot\underline{s}$, have been sampled. However, in many practical MX designs involving traps, the total number of vectors captured at a site is not known ahead of time. Rather, traps are set out for a fixed length of time $t$ and the catch total catch of units is a random variable, as are the number of size of pools/groups. If the size of pools/groups are fixed to be mostly to a single size ($s$), then the number of groups is $N_{tot} \approx V/s$. The exact relationship between catch size and number of pools, depends on what is done with vectors left over after distribution into groups of size $s$; however, we use the above approximation for mathematical convenience. We model the number of vectors trapped in time $t$ as a negative binomial distributed random variable with mean $t \mu$ and dispersion $tk$. Then the type II error is approximately

$$
\begin{align*}
E[\epsilon_1(\theta, \underline{s}, \underline{N}, \varphi, \psi)] &\approx E[\epsilon_1(\theta, s, V/s, \varphi, \psi)] \\
&= \sum_{x = 0}^\infty NB(x|t\mu, tk)\epsilon_1(\theta, s, x/s, \varphi, \psi) \\
&=  \left(\dfrac{k}{\mu + k - \epsilon_1(\theta, s, 1/s, \varphi, \psi)\mu}\right)^{t\mu},
\end{align*}
$$

where $NB$ is the probability mass function of a negative binomial distribution.

For a perfect test, the distribution of units between pools/groups and the number of pools does not matter than the type II error is exactly

$$
\begin{align*}
E[\epsilon_1(\theta, \underline{s}, \underline{N}, 1, 1)] 
&= \sum_{x = 0}^\infty NB(x|t\mu, tk) (1-\theta)^x \\
&=  \left(\dfrac{k}{k + \theta \mu}\right)^{t\mu}.
\end{align*}
$$



### Type II error using cluster surveys

The previous two sections we have assumed simple random sample from the population. However, many surveys including MX surveys, often involve a hierarchical or clustered sample design. We consider hierarchical survey designs with two or three levels. In a two-level survey, $M$ locations (e.g. households, buildings, breeding sites) are chosen by simple random sampling from a list of possible sampling sites and the vector populations at these sites are sampled. In a three level survey the area to be surveyed is split up into units (e.g. villages), a number ($M_2$) of these areas are selected randomly and within each of these areas a number ($M_1$) of locations are chosen to sample the vector populations. An example of a two-level survey might be to select households by simple random sampling from a list of households, and then place vector traps near selected households. An example of a three-level survey might be to might be to select villages by simple random sampling from a list of villages, then within each selected village select households by simple random sampling from a list of households, and then place traps near selected households.

If the prevalence of the marker of interest is clustered within sampling units (e.g. around households in a two-level survey or in villages and households in a three-level survey), a clustered survey will generally have a lower probability of detection than a simple random survey of vectors. To calculate the probability of detection for a survey, we need to have information (or make assumptions) about the kind and degree of clustering that exists at each level of sampling. We will model clustering in a manner analogous to random-effect logistic regression model, as this is the most common means of estimating modelling hierarchical survey data. This allows estimates of clustering from previous surveys to be used in the design of future surveys.

For a two-level survey we model the prevalence of the marker at each locations $\theta_1,...,\theta_M$ as identical and independently distributed logit-normal variables with mean $\theta$ and variance on the logit-scale of $\sigma^2$. i.e.

$logit(\theta_m) = log\left(\frac{\theta_m}{1- \theta_m}\right)\sim Normal(\theta', \sigma^2)$

with $\theta'$ such that $E[\theta_m] = \theta$. This is analogous to a random-intercept logistic regression model with normal random effects with variance $\sigma^2$.

In a two-level survey a type II error would occur when none of the detections at any of the sampling locations return a positive result even though the prevalence  $\theta>0$. If the size of the pools/groups ($s$) and the number of pools/groups per location ($N$) is fixed by design, the type II error probability is

$\left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma^2) \epsilon_1(f(x), s, N, \varphi, \psi)dx \right]^M $,

where $f$ is the inverse-logit (expit) function and $Normal$ is the density of the normal distribution.

If instead the time spent sampling at each location is fixed and we model the vector catch sizes at each location, $C_m$ as independent negative binomial random variables each with mean $t \mu$ and dispersion $tk$, then the number of pools/groups at each location is $N_m \approx C_m/s$, and the type II error probability is approximately
$$
\begin{align*}
\left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma^2) E[\epsilon_1(f(x), n, N_m, \varphi, \psi)] dx \right]^M & \approx 
\left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma^2) E[\epsilon_1(f(x), n, V_m/n, \varphi, \psi)] dx \right]^M \\
&=  \left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma^2)\left(\dfrac{k}{\mu + k - \epsilon_1(f(x), n, 1/n, \varphi, \psi)\mu}\right)^{t\mu} dx \right]^M.
\end{align*}
$$

The above integrals do not have neat analytical forms, however they can be approximated easily with numerical integration techniques. Note that type-II error decreases as $\sigma^2$ decreases, i.e. as the degree of clustering decreases. In the limiting case where the marker is not clustered at all ($\sigma^2 = 0$) the type II error probability for a survey with fixed sample size reduces to

$\epsilon_1(\theta, s, N, \varphi, \psi)^M  = \epsilon_1(\theta, s, NM, \varphi, \psi) $,

i.e. the type-I error probability for a simple random survey with $N\times M$ groups/pools. Similarly the type II error probability for survey with a fixed sampling time reduces to

$\left[\left(\dfrac{k}{\mu + k - \epsilon_1(\theta, s, 1/s, \varphi, \psi)\mu}\right)^{t\mu}\right]^M = \left(\dfrac{k}{\mu + k - \epsilon_1(\theta, s, 1/s, \varphi, \psi)\mu}\right)^{tM \mu}$

i.e. the type-II error probability of a simple random survey with a sampling time of $t\times M$.

We take a similar approach for a three-level hierarchical sampling design, with two sets of random effects. If $\sigma_1$ and $\sigma_2$ are the standard deviations at the two levels (e.g. household and village), then for a survey with a fixed number ($N$) of pools/groups at each of the $M_1 \times M_2$ sites has type I error is
$$
\left[\int_{-\infty}^{\infty}Normal(x|\theta',\sigma_2^2) \left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma_1^2)\epsilon_1(f(x), n, N, \varphi, \psi) dx \right]^{M_1}\right]^{M_2}.
$$

If the sample sizes are random, but sampling time at each of the $M_1\times M_2$ sites is fixed at $t$, then the type I error is approximately

$$
\left[\int_{-\infty}^{\infty}Normal(x|\theta',\sigma_2^2) \left[\int_{-\infty}^{\infty} Normal(x|\theta',\sigma_1^2)\left(\dfrac{k}{\mu + k - \epsilon_1(f(x), n, 1/n, \varphi, \psi)\mu}\right)^{t\mu} dx \right]^{M_1}\right]^{M_2}.
$$

### Type I error

A type I error is only possible if the test procedure has imperfect specificity, i.e. has some probability greater than zero of returning a positive test result on a pool/group of vectors that do not have the marker of interest. We make two assumptions:

* the test specificity ($\psi$) does not depend on the number of vectors per pool/group (s)
* false positive test results occur at random and independently across pools/groups

With these two assumptions it follows that the probability of getting at least one positive pool/group when the marker is absent ($\theta = 0$), depends only on the specificity ($\psi$) and the total number of groups/pools being tested: $(1-\psi)^{N_{tot}}$.

When sample sizes are fixed by design, for a simple random of sample of vectors $N_{tot} = N$, for a two-level cluster survey $N_{tot} = M N$, while for a three-level cluster survey $N_{tot} = M_1 M_2 N$. When sampling time/effort is fixed at $t$ but sample size is random, then $N_{tot}$ is also random. As before we model catch size at a site, $C$, as a negative binomial random variable with mean $t \mu$ and dispersion $tk$. For a two-level cluster survey $N_{tot} \approx \frac{1}{s} \sum_{m=1}^M C_m$. The sum over $C_m$ (the total catch size at all sites) is also negative binomially distributed, but with mean $tM\mu$ and dispersion $tMk$. Similarly the total catch size in a three-level cluster survey is negative binomially distributed with mean $tM_1M_2\mu$ and dispersion $tM_1M_2k$. So if $t_{tot}$ is the total sampling time across all locations in a survey — i.e. $t$ for a simple random survey, $tM$ for a cluster survey, and $tM_1M_2$ for a three-level cluster survey — then  type I error is approximately
$$
1 - \left(\dfrac{k}{\mu + k - \psi^{\frac{1}{n}}}\right)^{t_{tot}k}.
$$
For all the above designs, increasing specificity decreases the probability of a type I error while increasing the number of pools/groups increases the probability of a type II error. Therefore for a given sampling effort (i.e. either given $t_{tot}$ or $N_{tot}$), using larger pools/groups decreases the probability of a type I error. For survey designs with random catch sizes, the probability of a type I error increases as the dispersion parameters $k$ decreases i.e. as the variance of the catch size increases.
