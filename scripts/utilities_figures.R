plotCurve <- function(x, leg=names(x), type='pr', main=NULL, ngrid=500) {
  # Plot precision PR/ROC curves from PRROC using ggplot
  # args:
  #   x: a list of PR/ROC curves as returned by pr.curve/roc.curve.
  #   leg: Labels for each list entry of x
  #   type: one of pr/roc, specifying whether to plot PR/ROC curve.
  #   main: title of the plot
  #   ngrid: size of grid (recall values) to plot curve over
  require(ggplot2)
  require(dplyr)
  
  # Group curves by iRF type
  d <- lapply(1:length(x), function(i) {
    out <- data.frame(y=x[[i]]$curve[,2], x=x[[i]]$curve[,1], id=i)
    ngrid <- min(ngrid, nrow(out))
    id <- floor(seq(1, nrow(out), length.out=ngrid))
    return(out[id,])
  })
  
  # Get auc-info for each curve and set axis labels
  if (type == 'pr') {
    auc <- round(sapply(x, function(z) z$auc.integral), 2)
    x.lab <- 'Recall'
    y.lab <- 'Precision'
  } else if (type == 'roc') {
    auc <- round(sapply(x, function(z) z$auc), 2)
    x.lab <- 'FPR'
    y.lab <- 'TPR'
  } else {
    stop('type must be one of: pr, roc')
  }
  
  leg <- paste0(leg, ' (AUC = ', auc, ')')
  p <- do.call(rbind, d) %>% mutate(type=leg[id]) %>%
    ggplot(aes(x=x, y=y, col=type)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(main) +
    scale_color_discrete(name="") +
    theme_bw() +
    theme(legend.position=c(0.5, 0.15)) +
    theme(legend.title=element_blank()) +
    ylim(0:1)
  
  return(adjustAxes(p))
}

adjustAxes <- function(p) {
  require(ggplot2)
  p <- p + theme(plot.title=element_text(size=20),
                 axis.text=element_text(size=18),
                 axis.title=element_text(size=20),
                 strip.text=element_text(size=18,
                                         face="bold"),
                 legend.title=element_text(size=20),
                 legend.text=element_text(size=18)
  )
  return(p)
}





plotResponse <- function(pval.pcs, data = NULL, range1, range2){
  # Plot response surface of PCS p-value CART models for the null model (additive) and the alternative model (interacton)
  # args:
  # pval.pcs PCS p-value output from function pcsTestPredixCART (see "utilities_tests.R") with return.tree = TRUE
  # data (optional) should be a list with entried geno and pheno
                   # geno should be a matrix with two columns giving the genotype data
                   # pheno should be a vector with the phenotype data of same length as nrow(geno)
  # range1 numeric vector of length 2 giving the plotted range for the first feature (when missing the range from the input data is used)
  # range2 same as range1, but for the second feature
  # mean_pheno 
  library(tidyverse)
  library(RColorBrewer)
  library(ggpubr)
  
  colorScale <- scale_color_gradientn( colours = c(brewer.pal(n = 9, name = "Blues")[9:1],  
                                                   brewer.pal(n = 9, name = "Reds")), limits = c(0, 1) )
  colorScale_s <- scale_color_gradient(low = "#2171B5", high = "#CB181D", limits = c(0, 1) )
  colorScale_fill <- scale_color_gradientn( colours = c(brewer.pal(n = 9, name = "Blues")[8:1],  
                                                        brewer.pal(n = 9, name = "Reds")[1:8]), limits = c(0, 1), aesthetics = "fill")
  
  
  if(!is.null(data)){
    geno.sub <- data$geno
    pheno <- data$pheno
  }
  
  pv <- pval.pcs

  mean_pheno <- attr(pv, "meanPheno")
  tree <- attr(pv, "treeA")
  tree_comb <- attr(pv, "treeH")
  
  gnames <- attr(tree$terms, "term.labels")
  
  if(is.null(data)){
    if(missing(range1) | missing(range2)){
      stop("When no data is specified via geno.sub, then both, range1 and rang2, need to bespecified.")
    }
  }else{
    if(missing(range1)){
      range1 <- range(geno.sub[,1])
      
    }
    if(missing(range2)){
      range2 <- range(geno.sub[,2])
      
    }
  }

  
  grid1 <- seq(range1[1], range1[2], length.out = 200)
  grid2 <- seq(range2[1], range2[2], length.out = 200)
  
  n = length(grid1)^2
  geno_null <- matrix(nrow = n, ncol = 2)
  geno_null <- expand.grid(grid1, grid2)
  colnames(geno_null) <- gnames
  
  inter <- predict(tree, newdata = data.frame(geno_null))
  

  single <- lapply(1:length(tree_comb), function(x) predict(tree_comb[[x]],newdata = geno_null))
  
  poSi <- 2
  
  
  geno_null$y <- inter
  plotAB_tree <- ggplot(geno_null, aes_string(x = colnames(geno_null)[1], y = colnames(geno_null)[2], col= "y")) + 
    geom_point(size = poSi, shape = 15) +
    colorScale +
    ggtitle(paste0('CART(A,B) with p-value = 10^(-', round(-log10(as.numeric(pv)),0), ')')) +
    xlab(paste(colnames(geno_null)[1]) )+
    ylab(paste(colnames(geno_null)[2])) + 
    scale_fill_continuous(breaks = c(0, 1))
  
  if(length(single) > 1){
    geno_null$y <- pmin(1, pmax(0, do.call("+", single)  + mean_pheno))
  }else{
    geno_null$y <- pmin(1, pmax(0, do.call("+", single) ))
  }
  plotAB_sum <- ggplot(geno_null, aes_string(x = colnames(geno_null)[1], y = colnames(geno_null)[2], col= "y")) + 
    geom_point(size = poSi, shape = 15) +
    colorScale  +
    ggtitle('CART(A) + CART(B)') +
    xlab(paste(colnames(geno_null)[1]) )+
    ylab(paste(colnames(geno_null)[2])) + 
    scale_fill_continuous(breaks = c(0, 1)) 
  
  
  if(!is.null(data)){
    ind.random <- sample(1:nrow(geno.sub), nrow(geno.sub))
    data <- data.frame(geno.sub[ind.random, ], y = pheno[ind.random])
    
    plotScatter_all <- ggplot(data, aes_string(x = colnames(geno.sub)[1], y = colnames(geno.sub)[2], z= "y")) +
      colorScale_fill +
      stat_summary_hex(bins = 25, show.legend = F) +
      ggtitle(paste('data smoothed'))
    
    plot_combined <- ggarrange(plotAB_sum, plotAB_tree, plotScatter_all,
                               ncol = 3, nrow = 1, common.legend = T)
    
  }else{
    plot_combined <- ggarrange(plotAB_sum, plotAB_tree,
                               ncol = 2, nrow = 1, common.legend = T)
  }
  
  
  return(plot_combined)
  
}





plotPvalues <- function(pval.info, n, logL.full){
  # Plot summary figure with p-values (first column), stability score (second column), and prediction error of null and alternative model (third column)
  # args:
  # pval.info a data.frame with the following entries:
  # genes character vector with gene names of the interaction, genes separated by ' , '
  # pval numeric vector with p-values (between 0 and 1)
  # stab (optional) numeric vector with stability scores (between 0 and 1)
  # logLA numeric vector with log-likelihood for alternative model
  # logLH numeric vector with log-likelihood for null model
  # n (optional) single integer to rescale log-likelihood
  # logL.full (optimal) numeric value with log-likelihood of full model
  
  library(rpart.plot)
  library(tidyverse)
  library(RColorBrewer)
  library(ggpubr)
  library(patchwork)
  
  if(missing(n)){
    n <- 1
    namePA <- "Prediction Error \n - log P(Y|p)"
  }else{
    namePA <- "Prediction Error \n - log P(Y|p) / n "
  }
  
  pval.info <- pval.info %>%
    dplyr::mutate(order = str_count(genes, ",") + 1) %>% 
    dplyr::mutate(genes = paste0(genes,' (p-value = ', signif(pval,2), ')')) %>%
    dplyr::arrange(order, desc(pval) ) %>%
    dplyr::mutate(genes = fct_inorder(genes))
    
  
  tx_s <- 40
  #col_bar <- c("#E6AB02", "#FFF2AE")
  col_bar <- c("#009E73", "#999999")
  
  text_bar <- c("epistasis", "NULL")
  
  plotPV <-  pval.info %>% dplyr::select(genes, pval, order) %>%
    gather(key = type, value = pValue, -genes, -order) %>%  
    dplyr::mutate(pValue = -log10(pValue)) %>%
    ggplot(aes(y = genes, x=pValue, fill=type)) +
    geom_col(position="dodge", width = 0.5) + 
    scale_fill_manual(values = col_bar)+
    xlab("p-value (-log10)") +
    ylab("") + 
    theme(legend.title=element_blank(),
          text = element_text(size = tx_s),
          legend.position = "none",
          strip.text.y = element_blank()) +
    facet_grid(rows = vars(order), space = "free", scales = "free") 
  
  if('stab' %in% colnames(pval.info)){
    plotStab <- pval.info %>% dplyr::select(genes, stab, order) %>%
      gather(key = type, value = stability, -genes, -order) %>%
      ggplot(aes(y = genes, x=stability, fill=type))  +
      geom_col(position="dodge", width = 0.5) + 
      xlab("Stability") +
      ylab("") +
      coord_cartesian(xlim = c(0, 1))+ 
      scale_fill_manual(values = col_bar)+
      theme(legend.title=element_blank(),
            text = element_text(size = tx_s),
            legend.position = "none",
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.y = element_blank()) +
      facet_grid(rows = vars(order), space = "free", scales = "free") 
    
  }
    
  logL_min <- pval.info %>% dplyr::select(starts_with("logL")) %>% min 
  #logL_max <- pval.info %>% dplyr::select(starts_with("logL")) %>% max + 20
  logL_max <- 0
  
  data2 <- pval.info %>% dplyr::select(genes, logLA, order) %>%
    gather(key = type, value = prediction, -genes, -order) %>%
    mutate(prediction = - prediction / n)
  
  # plotPA <- pval.info %>% dplyr::select(genes, logLH, order) %>%
  #   gather(key = type, value = prediction, -genes, -order) %>%
  #   dplyr::mutate(prediction = - prediction / n) %>%
  #   ggplot(aes(y = genes, x=prediction, fill=type)) +
  #   geom_col(position="dodge", width = 0.5) + 
  #   geom_col(position="dodge", width = 0.5, data = data2, aes(y = genes, x=prediction, fill=type)) + 
  #   xlab(namePA) +
  #   ylab("") +
  #   scale_fill_manual(values = col_bar, labels = text_bar)+
  #   coord_cartesian(xlim = c(- logL_max / n, - logL_min / n))+ 
  #   theme(legend.title=element_blank(),
  #         text = element_text(size = tx_s, hjust = 0.2),
  #         axis.title.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.ticks.y=element_blank(),
  #         legend.text = element_text(size = tx_s),
  #         legend.key.height=unit(4,"line"),
  #         strip.text.y = element_text(angle = 0)) +
  #   facet_grid(rows = vars(order), space = "free", scales = "free") 
  
  
  
  plotPA <- pval.info %>% dplyr::select(genes, logLH, logLA, order) %>%
    gather(key = type, value = prediction, -genes, -order) %>%
    dplyr::mutate(prediction = - prediction / n) %>%
    ggplot(aes(y = genes, x=prediction, fill=type)) +
    geom_col(position="dodge", width = 0.5) + 
    xlab(namePA) +
    ylab("") +
    scale_fill_manual(values = col_bar, labels = text_bar)+
    coord_cartesian(xlim = c(- logL_max / n, - logL_min / n))+ 
    theme(legend.title=element_blank(),
          text = element_text(size = tx_s, hjust = 0.2),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text = element_text(size = tx_s),
          legend.key.height=unit(4,"line"),
          strip.text.y = element_text(angle = 0)) +
    facet_grid(rows = vars(order), space = "free", scales = "free") 
  
  
  if(!missing(logL.full) & !missing(n)){
    plotPA <- plotPA +
      geom_vline(xintercept = - logL.full/n)
  }
  
  if('stab' %in% colnames(pval.info)){
    plot_combined <- plotPV + plotStab + plotPA  + plot_layout(nrow = 1)
  }else{
    plot_combined <- plotPV + plotPA  + plot_layout(nrow = 1)
  }
  return(plot_combined)
}

