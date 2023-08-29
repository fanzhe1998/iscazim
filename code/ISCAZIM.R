## ISCAZIM
# Title: Benchmarking statistical correlation methods for zero-inflated microbiome-metabolome association analysis
# Version: 0.0.1
# Authors: Zhe Fan (fanzhe0308@163.com)
# Description: 
# Date: 2023/8/27

#####################################################################################
#                              Instructions
# The code was adapted on the function 'zeroinfl' from 'pscl' package,
# 'bootstrapMI' from 'maigesPack' package, 'mine' from 'minerva' package.
# Please see these packages for more details.
#####################################################################################
require(pscl)
require(psych)
require(compositions)
require(maigesPack)
require(minerva)

MIC <- function(x,y,R=100,...) {
  #MIC
  mic<-mine(x=x,y=y,...)
  
  #For calculate p-value de MIC
  if (! is.null(R)) {
    R <- floor(R)
    if (R < 1) R <- 100
  } else {
    R <- 100
  }
  Rep<-as.data.frame(rep(as.data.frame(y),R))
  Rep2<-as.data.frame(apply(Rep,2,sample))
  permic<-matrix(NA,nrow=R,ncol=7)
  colnames(permic)<-c("MIC","MAS","MEV","MCN","MIC-R2", "GMIC","TIC")
  for (i in 1:R){
    p<-mine(x=x,y=Rep2[,i],...)
    permic[i,1:7]<-c(p$MIC,p$MAS,p$MEV,p$MCN,p$`MIC-R2`,p$GMIC,p$TIC)
  }
  
  permic<-as.data.frame(permic)
  pvalor<-nrow(permic[which(permic$MIC>=mic$MIC),])/nrow(permic)
  
  mic2<-as.data.frame(cbind(mic$MIC,pvalor))
  colnames(mic2)<-c("MIC","p-value")
  return(mic2)
}

winsorization <- function (x, minval = NULL, maxval = NULL, winsor.qt = c(0, 0.97), 
          na.rm = FALSE, type = 7) 
{
  if (is.null(minval) || is.null(maxval)) {
    xq <- quantile(x = x, probs = winsor.qt, na.rm = na.rm, type = type)
    if (is.null(minval)) 
      minval <- xq[1L]
    if (is.null(maxval)) 
      maxval <- round(xq[2L],0)
  }
  x[x < minval] <- minval
  x[x > maxval] <- maxval
  return(x)
}

Univariate_ISCAZIM <- function(data, x, y, winsor = TRUE, winsor.qt = c(0,0.97), rep = 200){
  ZIR <- length(which(data[,y] == 0))/length(data[,y])
  
  # Winsorization
  if (winsor == TRUE) {
    otu.win <- winsorization(data[,y])
  } else {
    otu.win <- data[,y]
  }
  
  if(ZIR >= 0.2){
    zi <- zeroinfl(formula = otu ~ predictor, data = data.frame(otu = otu.win, predictor = data[,x]), dist = "negbin")
    linear_cofig <- ifelse(rownames(summary(zi)[["coefficients"]][["count"]])[2] == "Log(theta)", NA,
                        summary(zi)[["coefficients"]][["count"]][2,4])
  } else {
    pea <- cor.test(~ otu + predictor, data = data.frame(otu = otu.win, predictor = data[,x]), method = "pearson")
    linear_cofig <- pear[["p.value"]]
  }
  

  if(linear_cofig < 0.05){
    if(ZIR >= 0.2){
      zi <- zeroinfl(formula = otu ~ predictor, data = data.frame(otu = otu.win, predictor = data[,x]), dist = "negbin")
      cor_pvalue <- ifelse(rownames(summary(zl)[["coefficients"]][["count"]])[2] == "Log(theta)", NA,
                           summary(zl)[["coefficients"]][["count"]][2,4])
      cor_R <- ifelse(rownames(summary(zl)[["coefficients"]][["count"]])[2] == "Log(theta)", NA,
                         summary(zl)[["coefficients"]][["count"]][2,1])
      cor_pvalue <- data.frame(otu = colnames(data)[y], 
                               predictor = colnames(data)[x], 
                               linear = "Linear",
                               ZIR = paste(ZIR*100, "%", sep = ""),
                               method = "ZINB",
                               R = cor_R,
                               pvalue = cor_pvalue)
    } else {
      pea <- cor.test(~ otu + predictor, data = data.frame(otu = otu.win, predictor = data[,x]), method = "pearson")
      cor_R <- pear[["estimate"]][["cor"]]
      cor_pvalue <- pear[["p.value"]]
      cor_pvalue <- data.frame(otu = colnames(data)[y], 
                               predictor = colnames(data)[x], 
                               linear = "Linear",
                               ZIR = paste(ZIR*100, "%", sep = ""),
                               method = "Pearson",
                               R = cor_R,
                               pvalue = cor_pvalue)
    }
  } else {
    if(length(which(data[,y] == 0))/length(data[,y]) < 0.4){
      mic <- MIC(data[,y], data[,x])
      cor_R <- mic$"MIC"
      cor_pvalue <- mic$'p-value'
      cor_pvalue <- data.frame(otu = colnames(data)[y], 
                               predictor = colnames(data)[x], 
                               linear = "Nonlinear",
                               ZIR = paste(ZIR*100, "%", sep = ""),
                               method = "ZINB",
                               R = cor_R,
                               pvalue = cor_pvalue)
    } else {
      cor_pvalue <- bootstrapMI(data[,y], data[,x], bRep = rep, ret = "p-value")
      cor_R <- MI(data[,y], data[,x])
      cor_pvalue <- data.frame(otu = colnames(data)[y], 
                               predictor = colnames(data)[x], 
                               linear = "Nonlinear",
                               ZIR = paste(ZIR*100, "%", sep = ""),
                               method = "MI",
                               R = cor_R,
                               pvalue = cor_pvalue)
    }
  }
  return(cor_pvalue)
}

ISCAZIM <- function(data, x, y, winsor = TRUE, winsor.qt = c(0,0.97), rep = 200,n_cores = 1){
# Initializes the parallel environment:
  myCluster = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(myCluster)
  
  res <- foreach(y = c(1:length(y)),.combine='rbind') %dopar% Univariate_ISCAZIM()
    
  parallel::stopCluster(myCluster)
  res$FDR <- p.adjust(res$pvalue, method = "fdr")
  return(res)
}
