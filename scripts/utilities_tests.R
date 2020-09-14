###############################################################################
## general functions           
###############################################################################

fitTree <- function(pheno_ft, geno_ft, ind.it, cp_max = 0.01, rpart_method = "anova", all = TRUE){
  # Fit CART model
  # Input: pheno_ft numeric vector with phenotype values
  #        geno_ft  numeric matrix with genotype values
  #        ind.it   vector with column indexes of geno_ft that should be used for the CART model
  #        cp_max   initial complexity parameter 
  #        rpart_method method that should be used for rpart, default is "anova" 
  #        all  a logical value indicating whether all cp_max should be dynamically decreased such that all features appear in CART model
  library(rpart)
  cp <- cp_max
  if(all){
    valid <- FALSE
    while(!valid){
      gen <- as.matrix(geno_ft[, ind.it])
      colnames(gen) <- colnames(geno_ft)[ind.it]
      tree <- rpart(y ~ ., data = data.frame(y = pheno_ft, gen), method = rpart_method,
                    control = rpart.control(cp = cp, maxsurrogate = 0, xval = 0), model = T)
      valid <- (length(unique(tree$frame$var))-1 == length(ind.it)) | cp < 1e-3
      cp <- cp/1.1
    }
  }else{
    gen <- as.matrix(geno_ft[, ind.it])
    colnames(gen) <- colnames(geno_ft)[ind.it]
    tree <- rpart(y ~ ., data = data.frame(y = pheno_ft, gen), method = rpart_method,
                  control = rpart.control(cp = cp, maxsurrogate = 0, xval = 0), model = T)
  }
  return(tree)
}


gen_sub <- function(ind.it, geno){
  #extract ind.it columns of geno
  gen <- as.matrix(geno[, ind.it])
  colnames(gen) <- colnames(geno)[ind.it]
  return(gen)
}

mean_pattern <- function(x, pattern){
  #return means of x within pattern groups
  ans <- numeric(length(x))
  for(i in levels(pattern)){
    ind <- which(pattern == i)
    if(length(ind) > 0){
      ans[ind] <- mean(x[ind])
    }
  }
  return(ans)
}

###############################################################################
## tests for gene expression features 
###############################################################################

glmTestPredixMulti <- function(geno, geno.train, pheno, pheno.train, ind.int, fam = 'binomial', pca, return.fitted = F) {
  # Likelihood ratio test for gene-interaction
  # Test logistic or linear regression model with all interaction terms (NULL)
  # vs. logistic or linear regression model with all but the highest order interaction term
  # Input: geno numeric matrix with imputed gene expression data
  #        geno.train training data to compute prediction accuracy
  #        pheno numeric vector with phenotype values
  #        pheno.train training data to compute prediction accuracy
  #        ind.int vector with indexes of interaction features corresponding to columns in geno and geno.train
  #        fam family of glm, either 'binomial' or 'gaussian'
  #        pca numeric matrix with pca components of individulas are added as covariates in the model
  #        return.fitted a logical value indicating whether fitted models should be returned
  # Returns: p-value (numeric value between 0 and 1) with several attributes.
  library(stringr)
  library(glmnet)
  library(lmtest)
  library(PRROC)
  
  xint <- geno[,unlist(ind.int)]
  xint <- as.matrix(xint)
  
  xint.train <- geno.train[,unlist(ind.int)]
  xint.train <- as.matrix(xint.train)
  
  #interaction term to be tested
  inter <- apply(xint, 1, prod)
  inter.train <- apply(xint.train, 1, prod)
  
  #lower order interactions
  comb <- expand.grid(rep(list(c(0,1)), length(ind.int)))
  comb <- comb[rowSums(comb) > 0 & rowSums(comb) < length(ind.int),]
  
  xint <- apply(comb, 1, function(x) apply(xint, 1, function(y)  prod(y[x == 1])))
  xint.train <- apply(comb, 1, function(x) apply(xint.train, 1, function(y)  prod(y[x == 1])))
  
  
  d <- data.frame(y=pheno, x = xint, pca = pca, inter = inter)
  d.train <- data.frame(y=pheno.train, x = xint.train, inter = inter.train)
  
  # fit NULL models
  flin <- glm(y ~ ., data=d, family=fam)
  flin.train <- glm(y ~ ., data=d.train, family=fam)
  
  # coefficient and effect size
  if(is.na(tail(flin$coefficients,1))){
    co <- NA
    es <- NA
  }else{
    co <- tail(as.numeric(coef(flin)), n = 1)
    es <- co * mean(apply(xint, 1, prod))
  }
  
  # fit alternative models
  flin2 <- glm(y ~ . - inter, data = d, family = fam)
  flin2.train <- glm(y ~ . - inter, data = d.train, family = fam)
  
  # evaluate likelihood ratio test
  flr <- lrtest(flin, flin2)
  pv <- flr$`Pr(>Chisq)`[2]
  chisq <- flr$Chisq[2]
  
  # prediction accuracy on test data
  predA <- predict(flin.train, newdata = d, type = "response")
  predH <- predict(flin2.train, newdata = d, type = "response")
  
  # mse of null and alternative model
  mse1 <- mean((pheno - predA)^2)
  mse2 <- mean((pheno - predH)^2)
  
  # auc of null and alternative model
  auc1 <- roc.curve(predA[pheno == 1], predA[pheno == 0])$auc
  auc2 <- roc.curve(predH[pheno == 1], predH[pheno == 0])$auc
  
  # difference in mean prediction
  mDiff1 <- mean(predA[pheno == 1]) - mean(predA[pheno == 0])
  mDiff2 <- mean(predH[pheno == 1]) - mean(predH[pheno == 0])
  
  # log-likelihood
  logL_single <- sum(dbinom(pheno, size = 1, prob =  predH, log = T))
  logL_inter <- sum(dbinom(pheno, size = 1, prob =  predA, log = T))
  
  # log-likelihood ratio
  loglr <- logL_inter - logL_single
  
  attr(pv, "lr") <- exp(loglr)
  attr(pv, "log_lr") <- loglr
  attr(pv, "logLA") <- logL_inter
  attr(pv, "logLH") <- logL_single
  attr(pv, "mseA") <- mse1
  attr(pv, "mseH") <- mse2
  attr(pv, "aucA") <- auc1
  attr(pv, "aucH") <- auc2
  attr(pv, "mDiffA") <- mDiff1
  attr(pv, "mDiffH") <- mDiff2
  attr(pv, "effS") <- es
  attr(pv, "coef") <- co
  attr(pv, "chisq") <- chisq
  
  if(return.fitted){
    attr(pv, "fittedH") <- flin2$fitted.values
    attr(pv, "fittedA") <- flin$fitted.values
    attr(pv, "pS") <- predH
    attr(pv, "pI") <- predA
  }
  return(pv)
}



pcsTestPredixCART <- function(geno, pheno, ind.int, geno.train, pheno.train, return.fitted = F, return.tree = F) {
  # PCS p-value with CART model
  # Test P(Y | A,B) =  CART(A) + CART(B) model (NULL)
  # vs. P(Y | A,B) = CART(A,B) (epistasis)
  # Input: geno matrix with continuous gene expression data
  #        pheno binary vector
  #        ind.int vector with indeces of interaction features corresponding to columns in geno
  #        fam family of glm, either 'binomial' or 'gaussian'
  #        geno.train, pheno.train training data that is used to fit the trees
  #        pca matrix with pca components of individulas added as covariates in the model
  #        return.fitted logical whether fitted p vectors should be returned
  #        return.tree logical whether fitted trees should be returned
  # Returns: p-value (numeric value between 0 and 1) with several attributes.
  library(stringr)
  library(glmnet)
  library(rpart)
  library(lmtest)
  library(ranger)
  library(pracma)
  library(gtools)
  library(PRROC)
  
  
  rpart_method = "anova"
  
  geno.all <- geno.train
  pheno.all <- pheno.train
  
  xint <- geno[,unlist(ind.int)]
  xint <- as.matrix(xint)
  
  xint.train <- geno.train[,unlist(ind.int)]
  xint.train <- as.matrix(xint.train)
  
  #interaction term to be tested
  tree <- fitTree(pheno_ft = pheno.train, geno_ft = geno.train, ind.it = ind.int, rpart_method = rpart_method)
  cp_max <- tree$control$cp
  p.inter <- predict(tree, newdata = data.frame(y = pheno, xint))
  all_features <- (length(unique(tree$frame$var))-1 == length(ind.int))
  
  #lower order interactions
  comb <- expand.grid(rep(list(c(0,1)), length(ind.int)))
  comb <- comb[rowSums(comb) > 0 & rowSums(comb) < length(ind.int),]
  
  #fit individual trees with backfitting
  res <- pheno.all - mean(pheno.all)
  tree_comb <- vector('list', nrow(comb))
  tree_comb_pred <- vector('list', nrow(comb))
  diff_pred <- Inf
  max_iter <- 10
  j <- 0
  
  while(j < max_iter & diff_pred > 0.01){
    j <- j + 1
    tree_comb_pred_old <- tree_comb_pred
    perm_c <- gtools::permute(1:nrow(comb))
    for(i in 1:nrow(comb)){
      i_p <- perm_c[i]
      if(!is.null(tree_comb_pred[[i_p]])){
        res <- res + tree_comb_pred[[i_p]]
      }
      
      #backfitting step
      tree_comb[[i_p]] <- fitTree(pheno_ft = res, geno_ft = geno.all, ind.it = ind.int[comb[i_p, ] == 1], cp_max = cp_max, rpart_method = rpart_method)
      tree_comb_pred[[i_p]] <- predict(tree_comb[[i_p]], newdata = data.frame(gen_sub(ind.int[comb[i_p,] == 1], geno.all)))
      
      #mean centering of estimated function
      tree_comb_pred[[i_p]] <- tree_comb_pred[[i_p]] - mean(tree_comb_pred[[i_p]])
      
      res <- res - tree_comb_pred[[i_p]] 
    }
    
    if(!is.null(tree_comb_pred_old[[1]])){
      diff_pred <- max(sapply(1:nrow(comb), function(x) sqrt(mean((tree_comb_pred[[x]] - tree_comb_pred_old[[x]])^2))))
    }

  }
  
  
  # compute predictions on test set
  for(i in 1:nrow(comb)){
    tree_comb_pred[[i]] <- predict(tree_comb[[i]], newdata = data.frame(gen_sub(ind.int[comb[i,] == 1], geno))) 
    
    #mean centering of estimated function
    tree_comb_pred[[i]] <- tree_comb_pred[[i]] - mean(tree_comb_pred[[i]])
  }
  
  p.single <- mean(pheno.all) + Reduce('+', tree_comb_pred)
  
  ep <- 1e-10
  p.inter <- pmax(ep, pmin(1 - ep, p.inter))
  p.single <- pmax(ep, pmin(1 - ep, p.single))
  
  
  y <- pheno
  n <- length(y)
  weights <- log(p.inter) - log(1 - p.inter) - log(p.single) + log(1 - p.single)
  delta <- (p.single - y) * weights
  
  delta_mean <- mean(delta)
  delta_sd <- sd(delta)
  
  stat <- sqrt(n) * delta_mean / delta_sd
  pv <- pnorm(stat)
  
  if(!all_features){
    pv <- 1
  }
  
  lr_increase <- sum(dbinom(pheno, size = 1, prob =  p.inter, log = T)) > sum(dbinom(pheno, size = 1, prob =  p.single, log = T))
  if(!lr_increase){
    pv <- 1
  }
  
  
  # mse of NULL and alternative model
  mse1 <- mean((pheno - p.inter)^2)
  mse2 <- mean((pheno - p.single)^2)
  
  # auc of roc curve of NULL and alternative model
  roc_curve_inter <- roc.curve(p.inter[pheno == 1], p.inter[pheno == 0], curve = T)
  auc1 <- roc_curve_inter$auc
  roc_curve_single <- roc.curve(p.single[pheno == 1], p.single[pheno == 0], curve = T)
  auc2 <- roc_curve_single$auc
  
  # difference in mean prediction
  mDiff1 <- mean(p.inter[pheno == 1]) - mean(p.inter[pheno == 0])
  mDiff2 <- mean(p.single[pheno == 1]) - mean(p.single[pheno == 0])
  
  # log-likelihood
  logL_single <- sum(dbinom(pheno, size = 1, prob =  p.single, log = T))
  logL_inter <- sum(dbinom(pheno, size = 1, prob =  p.inter, log = T))
  
  # log-likelihood-ratio statistic
  loglr <- logL_inter - logL_single
  
  attr(pv, "lr") <- exp(loglr)
  attr(pv, "log_lr") <- loglr
  attr(pv, "logLA") <- logL_inter
  attr(pv, "logLH") <- logL_single
  attr(pv, "mseA") <- mse1
  attr(pv, "mseH") <- mse2
  attr(pv, "aucA") <- auc1
  attr(pv, "aucH") <- auc2
  attr(pv, "mDiffA") <- mDiff1
  attr(pv, "mDiffH") <- mDiff2

  attr(pv, "stat") <- stat
  attr(pv, "delta_mean") <- delta_mean
  attr(pv, "delta_sd") <- delta_sd
  
  if(return.tree){
    attr(pv, "treeA") <- tree
    attr(pv, "treeH") <- tree_comb
    attr(pv, "meanPheno") <- mean(pheno.all) 
  }
  
  if(return.fitted){
    attr(pv, "pS") <- p.single
    attr(pv, "pI") <- p.inter
    attr(pv, "pS_comb") <- tree_comb_pred
    attr(pv, "delta") <- delta
  }
  
  return(pv)
  
}


###############################################################################
## tests for SNP features            
###############################################################################

glmTestSnp <- function(geno, geno.train, pheno, pheno.train, ind.int, fam = 'binomial', model = 'full', pca, return.fitted = F) {
  # LR test for gene-interaction
  # Test logistic/linear regression model with all interaction terms (NULL)
  # vs. logistic/linear regression model with all but the higherst order interaction term
  # Input: geno numeric matrix with snp data taking values 0,1,2
  #        geno.train training data to evaluate prediction accuracy
  #        pheno numeric vector with phenotype values
  #        pheno.train training data to evaluate prediction accuracy
  #        ind.int vector with indeces of interaction features corresponding to columns in geno
  #        fam family of glm, either 'binomial' or 'gaussian'
  #        model one of 'full', 'dom', or 'rec' indicating whether 0,1,2 levels of SNPs should be merged into two levels.
  #              'full' keeps the three levels 0,1,2; 'dom' considers dominant model and merges 1,2 into one level; 
  #              'rec' considers recessive model and merges 0,1 into one level
  #        pca matrix with pca components of individulas added as covariates in the model
  #        return.fitted a logical value indicating whether fitted models should be returned
  # Returns: p-value (numeric value between 0 and 1) with several attributes.
  library(stringr)
  library(glmnet)
  library(lmtest)
  library(PRROC)
  
  # Construct dummy variables for 0,1,2 valued SNPs as in 
  #(Cordell, "Epistasis: what it means, what is doesn't man, and statistical methods to detect it in humans', Human Molecular Genetics, 2002, Vol. 11, No. 20)
  
  if(model == 'full'){
    xVal <- function(g) c(-1, 0, 1)[g+1]
    zVal <- function(g) c(-0.5, 0.5, -0.5)[g+1]
    
    xint <- geno[,unlist(ind.int)]
    xint <- cbind(apply(xint, c(1,2), xVal), apply(xint, c(1,2), zVal) )
    
    xint.train <- geno.train[,unlist(ind.int)]
    xint.train <- cbind(apply(xint.train, c(1,2), xVal), apply(xint.train, c(1,2), zVal) )
  }
  if(model == 'dom'){
    xVal <- function(g) c( 0,1,1)[g+1]
    
    xint <- geno[,unlist(ind.int)]
    xint <- apply(xint, c(1,2), xVal)
    
    xint.train <- geno.train[,unlist(ind.int)]
    xint.train <- apply(xint.train, c(1,2), xVal)
  }
  if(model == 'rec'){
    xVal <- function(g) c( 0,0,1)[g+1]
    
    xint <- geno[,unlist(ind.int)]
    xint <- apply(xint, c(1,2), xVal)
    
    xint.train <- geno.train[,unlist(ind.int)]
    xint.train <- apply(xint.train, c(1,2), xVal)
  }
  
  xint <- as.matrix(xint)
  xint.train <- as.matrix(xint.train)
  
  
  if(model == 'full'){
    #interaction term(s) to be tested
    comb <- expand.grid(rep(list(c(0,1)), length(ind.int)))
    inter <- apply(comb, 1, function(c) apply( xint[, as.numeric(seq(1, 2 * length(ind.int), 2) + c)], 1, prod))
    inter.train <- apply(comb, 1, function(c) apply( xint.train[, as.numeric(seq(1, 2 * length(ind.int), 2) + c)], 1, prod))
    
    
    #lower order interactions
    comb <- expand.grid(rep(list(c(0,1)), length(ind.int)))
    comb <- comb[rowSums(comb) > 0 & rowSums(comb) < length(ind.int),]
    comb <- t(do.call(cbind, lapply(1:nrow(comb), 
                                   function(c) apply(expand.grid(rep(list(c(1,2)), sum(comb[c,]))), 1, 
                                                     function(cc) comb_inout(comb[c,], cc)))))
    xint <- apply(comb, 1, function(x) apply(xint, 1, function(y)  prod(y[x == 1])))
    xint.train <- apply(comb, 1, function(x) apply(xint.train, 1, function(y)  prod(y[x == 1])))
    
    
  }else{
    inter <- matrix(apply(xint, 1, prod), nrow = nrow(geno))
    inter.train <- matrix(apply(xint.train, 1, prod), nrow = nrow(geno.train))
    
    #lower order interactions
    comb <- expand.grid(rep(list(c(0,1)), length(ind.int)))
    comb <- comb[rowSums(comb) > 0 & rowSums(comb) < length(ind.int),]
    xint <- apply(comb, 1, function(x) apply(xint, 1, function(y)  prod(y[x == 1])))
    
    xint.train <- apply(comb, 1, function(x) apply(xint.train, 1, function(y)  prod(y[x == 1])))
    
  }
  
  d <- data.frame(y = pheno, x = xint, pca = pca, inter = inter)
  d.train <- data.frame(y = pheno.train, x = xint.train, inter = inter.train)
  
  # fit null model (no epistasis)
  flin <- glm(y ~ ., data = d, family = fam, x = T)
  flin.train <- glm(y ~ ., data = d.train, family = fam)
  
  co <- tail(coef(flin), n = ncol(inter))
  es <- co * tail(colMeans(flin$x), n = ncol(inter))
  
  d <- data.frame(y = pheno, x = xint, pca = pca)
  d.train <- data.frame(y = pheno.train, x = xint.train)
  
  # fit alternative model (epistasis)
  flin2 <- glm(y ~ ., data = d, family = fam)
  flin2.train <- glm(y ~ ., data = d.train, family = fam)
  
  # evaluate likelihood ratio test
  flr <- lrtest(flin, flin2)
  pv <- flr$`Pr(>Chisq)`[2]
  chisq <- flr$Chisq[2]
  
  # predicted response
  predH <- predict(flin2.train, newdata = data.frame(y = pheno, x = xint, inter = inter), type = "response")
  predA <- predict(flin.train, newdata = data.frame(y = pheno, x = xint, inter = inter), type = "response")
  
  # mse
  mse1 <- mean((predA - pheno)^2)
  mse2 <- mean((predH - pheno)^2)
  
  # auc of roc curve
  auc1 <- roc.curve(predA[pheno == 1], predA[pheno == 0])$auc
  auc2 <- roc.curve(predH[pheno == 1], predH[pheno == 0])$auc
  
  # difference in means
  mDiff1 <- mean(predA[pheno == 1]) - mean(predA[pheno == 0])
  mDiff2 <- mean(predH[pheno == 1]) - mean(predH[pheno == 0])
  
  # log likelihood
  logL_single <- sum(dbinom(pheno, size = 1, prob =  predH, log = T))
  logL_inter <- sum(dbinom(pheno, size = 1, prob =  predA, log = T))
  
  # log-likelihood-ratio
  loglr <- logL_inter - logL_single
  
  
  attr(pv, "lr") <- exp(loglr)
  attr(pv, "log_lr") <- loglr
  attr(pv, "logLA") <- logL_inter
  attr(pv, "logLH") <- logL_single
  attr(pv, "mseA") <- mse1
  attr(pv, "mseH") <- mse2
  attr(pv, "aucA") <- auc1
  attr(pv, "aucH") <- auc2
  attr(pv, "mDiffA") <- mDiff1
  attr(pv, "mDiffH") <- mDiff2
  attr(pv, "effS") <- es
  attr(pv, "coef") <- co
  attr(pv, "chisq") <- chisq
  if(return.fitted){
    attr(pv, "fittedH") <- flin2
    attr(pv, "fittedA") <- flin
    attr(pv, "pS") <- predH
    attr(pv, "pI") <- predA
  }

  return(pv)
}




pcsTestSnp <- function(geno, pheno,ind.int, geno.train, pheno.train, model = 'full', return.fitted = F) {
  # PCS p-values for SNP data for binary phenotypes
  # Input: geno numeric matrix with snp data taking values 0,1,2
  #        geno.train training data to evaluate prediction accuracy
  #        pheno binary vector with phenotype values
  #        pheno.train training data to evaluate prediction accuracy
  #        ind.int vector with indeces of interaction features corresponding to columns in geno
  #        model one of 'full', 'dom', or 'rec' indicating whether 0,1,2 levels of SNPs should be merged into two levels.
  #              'full' keeps the three levels 0,1,2; 'dom' considers dominant model and merges 1,2 into one level; 
  #              'rec' considers recessive model and merges 0,1 into one level
  #        return.fitted a logical value indicating whether fitted models should be returned
  # Returns: p-value (numeric value between 0 and 1) with several attributes.
  library(stringr)
  library(glmnet)
  library(lmtest)
  library(PRROC)
  library(pracma)
  library(gtools)
  
  xint <- geno[,unlist(ind.int)]
  xint.train <- geno.train[,unlist(ind.int)]

  if(model == 'dom'){
    xVal <- function(g) c( 0,1,1)[g+1]
    xint <- apply(xint, c(1,2), xVal)
    xint.train <- apply(xint.train, c(1,2), xVal)
  }
  if(model == 'rec'){
    xVal <- function(g) c( 0,0,1)[g+1]
    xint <- apply(xint, c(1,2), xVal)
    xint.train <- apply(xint.train, c(1,2), xVal)
  }
  
  xint <- as.matrix(xint)
  xint.train <- as.matrix(xint.train)
  
  xint.all <- xint.train
  pheno.all <- pheno.train

  #fit individual probabilities with backfitting
  res <- pheno.all - mean(pheno.all)
  fit_comb <- vector('list', length(ind.int))
  pattern_comb <- vector('list', length(ind.int))
  for(i in 1 : length(ind.int)){
    pattern_comb[[i]] <- as.factor(apply(xint.all, 1, function(x) paste(x[-i], collapse = "")))
  }
  
  diff_pred <- Inf
  max_iter <- 10
  j <- 0
  
  while(j < max_iter & diff_pred > 0.01){
    j <- j + 1
    fit_comb_old <- fit_comb
    perm_c <- gtools::permute(1:length(ind.int))
    for(i in 1: length(ind.int) ){
      i_p <- perm_c[i]
      
      if(!is.null(fit_comb[[i_p]])){
        res <- res + fit_comb[[i_p]]
      }
      
      #backfitting step
      fit_comb[[i_p]] <- mean_pattern(res, pattern_comb[[i_p]])
      #mean centering of estimated function
      fit_comb[[i_p]] <- fit_comb[[i_p]] - mean(fit_comb[[i_p]])
      
      res <- res - fit_comb[[i_p]] 
    }
    
    if(!is.null(fit_comb_old[[1]])){
      diff_pred <- max(sapply(1:length(ind.int), function(x) sqrt(mean((fit_comb[[x]] - fit_comb_old[[x]])^2))))
    }
    
  }
  
  for(i in 1 : length(ind.int)){
    pattern <- as.factor(apply(xint, 1, function(x) paste(x[-i], collapse = "")))
    fit <- numeric(length(pheno))
    for(j in levels(pattern)){
      if(sum(pattern_comb[[i]] == j) > 0){
        fit[pattern == j] <- fit_comb[[i]][pattern_comb[[i]] == j][1]
      }
    }
    fit_comb[[i]] <- fit
  }
  
  p.single <- mean(pheno.all) + Reduce('+', fit_comb)
  
  
  #fit interaction term
  pattern_train <- as.factor(apply(xint.train, 1, function(x) paste(x, collapse = "")))
  fit_int <- mean_pattern(pheno.train, pattern_train)
  
  pattern <- as.factor(apply(xint, 1, function(x) paste(x, collapse = "")))
  fit <- numeric(length(pheno))
  for(j in levels(pattern)){
    ind <- which(pattern == j)
    if(length(ind) > 0 & length(which(pattern_train == j)) > 0){
      fit[ind] <- fit_int[which(pattern_train == j)][1]
    }
  }
  p.inter <- fit
  
  ep <- 1e-10
  p.inter <- pmax(ep, pmin(1 - ep, p.inter))
  p.single <- pmax(ep, pmin(1 - ep, p.single))
  
  y <- pheno
  n <- length(y)
  weights <- log(p.inter) - log(1 - p.inter) - log(p.single) + log(1 - p.single)
  delta <- (p.single - y) * weights
  
  delta_mean <- mean(delta)
  delta_sd <- sd(delta)
  
  stat <- sqrt(n) * delta_mean / delta_sd
  pv <- pnorm(stat)
  
  lr_increase <- sum(dbinom(pheno, size = 1, prob =  p.inter, log = T)) > sum(dbinom(pheno, size = 1, prob =  p.single, log = T))
  if(!lr_increase){
    pv <- 1
  }
  
  mse1 <- mean((pheno - p.inter)^2)
  mse2 <- mean((pheno - p.single)^2)
  
  auc1 <- roc.curve(p.inter[pheno == 1], p.inter[pheno == 0], curve = T)$auc
  auc2 <- roc.curve(p.single[pheno == 1], p.single[pheno == 0], curve = T)$auc
  
  mDiff1 <- mean(p.inter[pheno == 1]) - mean(p.inter[pheno == 0])
  mDiff2 <- mean(p.single[pheno == 1]) - mean(p.single[pheno == 0])
  
  logL_single <- sum(dbinom(pheno, size = 1, prob =  p.single, log = T))
  logL_inter <- sum(dbinom(pheno, size = 1, prob =  p.inter, log = T))
  
  loglr <- logL_inter - logL_single
  
  attr(pv, "lr") <- exp(loglr)
  attr(pv, "log_lr") <- loglr
  attr(pv, "logLA") <- logL_inter
  attr(pv, "logLH") <- logL_single
  attr(pv, "mseA") <- mse1
  attr(pv, "mseH") <- mse2
  attr(pv, "aucA") <- auc1
  attr(pv, "aucH") <- auc2
  attr(pv, "mDiffA") <- mDiff1
  attr(pv, "mDiffH") <- mDiff2
  attr(pv, "stat") <- stat
  attr(pv, "delta_mean") <- delta_mean
  attr(pv, "delta_sd") <- delta_sd
  
  if(return.fitted){
    attr(pv, "fitA") <- fit
    attr(pv, "fitH") <- fit_comb
    attr(pv, "pS") <- p.single
    attr(pv, "pI") <- p.inter
  }
  
  return(pv)
}
