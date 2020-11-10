#'---
#' title: "PCS inference for linear regression"
#' output: html_document
#' ---

#+ setup, message=FALSE
library(tidyselect)
library(ggplot2)
set.seed(193)

#' ## Setup
#' The functions below calculate classical and PCS p-values in a simple linear 
#' regression setting given feature matrix `x` and response vector `y` 
#' without intercept and with unit variance.
#+ inference_functions
z.test <- function(y, x){
  # Generates classical p-value from feature matrix x and response vector y
  return(pnorm(-sum(x * y) / sqrt(sum(x ^ 2))))
}

pcs.test <- function(y, x){
  # Generates PCS p-value from feature matrix x and response vector y
  y_test <- y[1:(n / 2)]
  x_test <- x[1:(n / 2)]
  
  y_train <- y[-(1:(n / 2))]
  x_train <- x[-(1:(n / 2))]
  
  c_train <- max(sum(x_train * y_train) / sum(x_train ^ 2), 0)
  
  if (sum(y_test ^ 2) < sum((y_test - c_train * x_test) ^ 2)) {
    return(1)
  } else{
    y_null <- rnorm(n / 2)
    B <- 1000
    pval_B <- numeric(B)
    for(b in 1:B){
      indB <- sample(1:(n / 2), n / 2, replace = T)
      pval_B[b] <- sum(x_test[indB] * y_test[indB]) <= sum(x_test[indB] * y_null[indB]) 
    }
    return( (sum(pval_B) + 1) / B)
  }
}

#' The paramters below are used across all simulation settings
#+ global <- params
n <- 1000 # number of samples
N <- 1000 # number of simulation replicates
pval <- numeric(N) #classical p-value
pval_pcs <- numeric(N) #PCS p-value

#' ## Models
#' #### Correctly specified H1
#' Below, we simulate under a correctly specified linear regression model with
#' coefficient `c`.
#+ models_1
c <- 4
for(i in 1:N){
  x <- runif(n)
  x <- x/sqrt(sum(x^2))
  
  y <- rnorm(n, c * x) #correctly specified
  pval[i] <- z.test(y,x)
  pval_pcs[i] <- pcs.test(y,x)
}

pval_correct <- pval
pval_pcs_correct <- pval_pcs


#' #### Correctly specified H0
#' Below, we simulate under a correctly specified null model i.e. no association
#' between `x` and `y`.
#+ models_2
for(i in 1:N){
  y <- rnorm(n) #correctly specified
  pval[i] <- z.test(y,x)
  pval_pcs[i] <- pcs.test(y,x)
}

pval_correct_null <- pval
pval_pcs_correct_null <- pval_pcs


#' #### Misspecified H0
#' Below, we run inference under a misspecified null. Responses `y` are
#' generated as Gaussian mean = 0,  standard deviation = $\sigma$ random variables i.e. no
#' assocation between `x` and `y`. P-values are computed relative to a Gaussian mean = 0,
#' standard deviation = 1 distribution.
#+ models_3
Sigma <- c(2, 10, 100)

pval_var <- vector('list', length(Sigma))
pval_pcs_var <- vector('list', length(Sigma))
for(j in 1:length(Sigma)){
  for(i in 1:N){
    x <- runif(n)
    x <- x/sqrt(sum(x^2))
    
    y <- rnorm(n, 0, Sigma[j]) #large variance
    pval[i] <- z.test(y,x)
    pval_pcs[i] <- pcs.test(y,x)
  }
  pval_var[[j]] <- pval
  pval_pcs_var[[j]] <- pval_pcs
}

#' #### Misspecified H1
#' Below, we run inference under a misspecified alternative. The majority $(1-\epsilon)$ of
#' responses `y` are generated under a simple linear model with coefficient `c`.
#' A small proportion $(\epsilon)$ of responses are generated using a different
#' coefficient `-Cvar`. 
#+ models_4

c <- 10
Cvar <- c(100, 10000)
Eps <- rbind(c(0.15, 0.09, 0.04),
             c(0.01, 0.001, 0.0002))

pval_mix <- vector('list', length(Eps))
pval_pcs_mix <- vector('list', length(Eps))

for(l in 1:length(Cvar)){
  for(j in 1:ncol(Eps)){
    for(i in 1:N){
      x <- runif(n)
      x <- x/sqrt(sum(x^2))
      
      ind_mix <- rbinom(n, 1, Eps[l,j] ) == 1
      y <- rnorm(n, c * x) 
      y[ind_mix] <- rnorm(sum(ind_mix), - Cvar[l] * x[ind_mix] )
      
      pval[i] <- z.test(y,x)
      pval_pcs[i] <- pcs.test(y,x)
      
    }
    pval_mix[[(l - 1) * ncol(Eps) + j]] <- pval
    pval_pcs_mix[[(l - 1) * ncol(Eps) + j]] <- pval_pcs
  }
}


#' #### Details on correctly specified $H_1$ case
#' 
#'Assume that the true responses are drawn independently from a probabilistic model such that 
#'\begin{align*}
#'  y_i \; \overset{ind.}{\sim} \; N(c x_i,1) \text{ for some } c > 0
#'\end{align*}
#'and analog for $\tilde{y}_i$ (with the same constant $c$). 
#'Then it is easy to
#'show that the classical p-value in a simple linear regression model is distributed as
#' \begin{align*}
#' \text{ p-value } = \Phi(-c + N(0,1)),
#' \end{align*}
#' such that for large values of $c$ it holds approximately that
#' \begin{align}
#' -\log(\text{ p-value }) \approx c^2.
#' \end{align}
#' Next, we compare this with the PCS p-value. When both, the sample size $n$
#' and the constant $c$, are sufficiently large, it holds true that $\hat{c}
#' \approx \sum_{i = 1}^n x_i y_i \approx c$ and thus, the prediction screening
#' is unlikely to yield a PCS p-value of one. Therefore, for the PCS p-value we
#' find that 
#' \begin{align*}
#' \text{PCS p-value } \approx P\left( \sum_{i = 1}^n x_i \left( y_i - {y}_{0,i} \right) < 0 \right) = P\left( c + N(0,2) > 0 \right) = \Phi\left(\frac{- c}{\sqrt{2}} \right)
#' \end{align*}
#' and thus, for large values of $c$ it holds approximately that
#' \begin{align}
#' -\log(\text{ PCS p-value }) \approx \frac{c^2}{2}.
#' \end{align}
#' Comparing the p-values above,  we see that for a correctly specified $H_1$
#' model, the PCS p-value has a decreased detection power.
#' From the above analysis we see that this is because exploring 
#' the full empirical distribution of the test statistic $T(Y)$ 
#' yields an additional variance term, which decreases the power.
#' 
#' #### Details on misspecified $H_0$ case
#' Assume that the true responses are drawn from a probabilistic model such that 
#' \begin{align*}
#' y_i \overset{i.i.d.}{\sim} \mathcal{N}(0,\sigma^2) \text{ for some } \sigma^2  > 1
#' \end{align*} 
#' and analog for $\tilde{y}_i$ (with the same $\sigma^2$).
#' Then for the normal p-value we obtain
#' \begin{align*}
#' \text{ p-value } = \Phi(- \mathcal{N}(0,\sigma^2)) \Longrightarrow \text{Bernoulli}(0.5) \quad \text{ as } \sigma^2 \to \infty.
#' \end{align*}
#' This implies that for any significance level $\alpha \in (0,1)$ and sufficiently large $\sigma^2$, 
#' the classical p-value rejects the null model $H_0$ with probability approximately $50\%$.
#' The PCS p-value, in contrast, after prediction screening, can be approximated as
#' \begin{align*}
#' \text{ PCS p-value } = P\left(\sum_{i = 1}^n x_i (y_i - y_{i,0}) < 0 \right) = P\left(\mathcal{N}(0, \sigma^2 + 1) < 0 \right)
#' \end{align*}
#' and hence, the PCS p-value follows a sub-uniform distribution, as expected under the null.
#' 
#' #### Details on misspecified $H_1$ case
#' Finally, we consider responses that are drawn from a mixture distribution such that
#' \begin{align*} 
#' y_i \overset{ind.}{\sim} (1 - \epsilon) \; \mathcal{N}(c \; x_i, 1) + \epsilon \;  \mathcal{N}(- c^\prime \; x_i, 1),
#' \end{align*}
#' for two constants $c, c^\prime > 0$ and some small $\epsilon > 0$.
#' This means that, on average, a $(1 - \epsilon)$ fraction of all observations follow 
#' the correctly specified alternative model $H_1$ with coefficient $c > 0$. 
#' However, on average an $\epsilon$ fraction of the observations are outliers, with a negative 
#' regression coefficient $- c^\prime < 0$.
#' 

#'
#' ## Results
#' The results below compare classical and PCS p-values under correctly
#' specified and misspecified models. In summary, we find that the classical
#' p-value outperforms PCS p-value when the data generative models are specified 
#' correctly. However, the PCS p-value is more robust to settings where the 
#' data generating process is misspecified.
#'
#' #### Correct H0: 
#' When data are generated under the null,
#' both the classical p-value and the PCS p-value follow a sub-uniform
#' distribution.  While the classical p-value exactly follows a uniform
#' distribution, the distribution of the PCS p-value is sub-uniform due to the
#' prediction screening step. For example, while the classical p-value is
#' smaller than $0.5$ in $50\%$ of cases, the PCS p-value is smaller than $0.5$
#' in only $\sim 40\%$ of cases. However, for smaller quantiles (that are of
#' dominant interest in practice) the different between the PCS p-value and the
#' classical p-value distributions is almost negligible under the exact null
#' distribution.
#'
#' #### Correct H1: 
#' When data are generated under the alternative, we find that the
#' PCS p-value has less power than the classical p-value. For example, at a
#' nominal level of $\alpha = 0.05$ the classical p-value correctly rejects the
#' null hypothesis in $99\%$ of cases, but the PCS p-value rejects in only
#' $64\%$ of cases. Similarly, at nominal level $\alpha = 0.1$ the classical
#' p-value rejects in $99.7\%$ of cases, but the PCS p-value in only $74.9\%$ of
#' cases. Intuitively, the reason for this is that exploring the full
#' distribution of the test statistic $T(Y)$ via the bootstrap sampling results
#' in an additional variance term, which leads to a slightly weaker detection
#' power when the observed responses exactly follow the specified alternative
#' model $H_1$.
#'
#' #### Misspecified H0: 
#' When data are generated under a misspecified version of
#' the null, we find that the classical p-value yields severe false positives.
#' For example, when $\sigma=10$ at a nominal level of $\alpha = 0.05$ the classical p-value
#' falsely rejects the null hypothesis in $42\%$ of cases. In contrast, the PCS
#' p-value, via the prediction screening as well as the bootstrap sampling, is
#' robust to this misspecification. E.g., at nominal level $\alpha = 0.05$ the
#' PCS p-value falsely rejects the null hypothesis in only $5\%$ of cases, as is
#' expected under the null. Similarly, at nominal level $\alpha = 0.01$ the PCS
#' p-value rejects in $\%1$ of cases, but the classical p-value still rejects in
#' as much as $40\%$ of cases.
#'
#' #### Misspecified H1:
#' First, consider the situation where $c^\prime$ is of moderate size.
#' On the one hand, for small values $\epsilon$ the average regression 
#' coefficient $(1 - \epsilon) c - \epsilon c^\prime$ is still positive, although 
#' the alternative model is misspecified. 
#' In that case, one observes that the clasical p-value still has a higher detection power than the 
#' PCS p-value, analog as in the correctly specified case.
#' On the other hand, for moderately sized values of $\epsilon$ 
#' the average regression coefficient $(1 - \epsilon) c - \epsilon c^\prime$ is negative 
#' and there is also a relatively large fraction of individual observations 
#' with a negative regression coefficient.
#' In that case, both the clasical p-value and the PCS p-value follow a sub-uniform distribution, 
#' i.e., the null hypothesis will not be rejected on average.
#' 
#' Second, consider the situation where $c^\prime$ is very large.
#' For small value of $\epsilon$, we find that the PCS p-value 
#' has significantly larger power than the clasical p-value.
#' This corresponds to a situation where there are very few, but extreme outliers.
#' Note that in this case on 
#' average for a $(1-(1-\epsilon)^{1000})$ fraction
#' (= probability that a sample of $n = 1,000$ data points contains at least one outlier) 
#' of the Monte Carlo runs the test statistic is dominated 
#' by the outliers with very small regression coefficient  $- c^\prime$, which leads to 
#' a classical p-value of approximately $1$.
#' On the other hand, the PCS p-value explores the full distribution 
#' of the test statistic via the bootstrap sub-sampling. 
#' Because the outliers are very rare, most of the bootstrap sub-samples will not contain an outlier 
#' and therefore, the PCS p-value obtains a higher detection power compared to the classical p-value, 
#' whenever $c^\prime$ is sufficiently large (i.e., a strong model misspecification).

#'  For example, when $\epsilon = 0.001$, $c^\prime = 10000$, and $c = 10$,
#'  at a nominal level of $\alpha = 0.05$ the
#' classical p-value correctly rejects the null hypothesis only in $62\%$ of
#' cases, but the PCS p-value rejects in $69\%$ of cases.
#' 
#+ results, echo=FALSE
# Aggregate results for correctly specified models
df_correct <- data.frame(pval = c(pval_correct, pval_pcs_correct), 
                         type = c(rep("zTest", N), rep("PCS", N)))
df_correct_null <- data.frame(pval = c(pval_correct_null, pval_pcs_correct_null), 
                              type = c(rep("zTest", N), rep("PCS", N)))

#Aggregate results for misspecified models
df_var <- data.frame(pval = c(pval_var[[2]], pval_pcs_var[[2]]), 
                       type = c(rep("zTest", N), rep("PCS", N)))
df_mix <- data.frame(pval = c(pval_mix[[5]], pval_pcs_mix[[5]]), 
                       type = c(rep("zTest", N), rep("PCS", N)))


# all parameter settings
df_var_all <- data.frame(pval = c(unlist(pval_var), unlist(pval_pcs_var)), 
                     type = c(rep("ztest", N * length(Sigma)), 
                              rep("PCS", N * length(Sigma))),
                     sigma = c(rep(rep(as.character(Sigma), each = N), 2)))

Eps_name <- t(matrix(paste0("epsilon", 1:ncol(Eps)), ncol = length(Cvar), nrow = ncol(Eps)))

df_mix_all <-  data.frame(pval = c( unlist(pval_mix) , unlist(pval_pcs_mix)), 
                  type = c( rep( "zTest", N * ncol(Eps) * length(Cvar) ), 
                            rep( "PCS", N * ncol(Eps) * length(Cvar))),
                  epsilon = rep( rep( as.character( t(Eps_name) ), each = N), 2 ),
                  constant = rep( rep(paste( "c'  =", Cvar )  , each = N * ncol(Eps) ), 2 ))



# Aggregate results across all models for plotting
df <- list(df_correct_null,
           df_correct, 
           df_var, df_mix,
           df_var_all, df_mix_all)

title <- c("Correct H0",
           "Correct H1 (c = 4)", 
           expression(paste("Misspecified H0 (", sigma, " = 10 )" )), 
           expression(paste("Misspecified H1 ( c' = 10000 and ", epsilon, " = 0.001 )")) ,
           "Misspecified H0 - multiple parameter settings", 
           "Misspecified H1 - multiple paramter settings")

#+ figures, echo=FALSE, fig.width=12, fig.height=8
# Generate plots from each model setting

for(i in 1:length(df)){
  if (i <= 4) {
    pl <- ggplot(df[[i]], aes(pval, linetype = type)) 
  }
  if (i == 5) {
    pl <- ggplot(df[[i]], aes(pval, linetype=type, colour=sigma)) + 
      labs(color  = expression(sigma)) 
  }
  if (i == 6) {
    pl <- ggplot(df[[i]], aes(pval, linetype = type, colour = epsilon)) + 
      labs(color  = expression(epsilon)) + 
      facet_grid( cols = vars(constant)) 
  }
  pl <- pl +
    stat_ecdf(size = 1) +
    geom_line(aes(x=pval, y=pval), colour="darkgray", size = 1)  +
    labs(x="t", y ="P( p-value \u2264 t )", title = title[i])

  plot(pl)
}
