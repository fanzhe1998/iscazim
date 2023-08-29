MIC <- function(x,y,R=100,...) {
  #MIC
  mic<-minerva::mine(x=x,y=y,...)
  
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
    p<-minerva::mine(x=x,y=Rep2[,i],...)
    permic[i,1:7]<-c(p$MIC,p$MAS,p$MEV,p$MCN,p$`MIC-R2`,p$GMIC,p$TIC)
  }
  
  permic<-as.data.frame(permic)
  pvalor<-nrow(permic[which(permic$MIC>=mic$MIC),])/nrow(permic)
  
  mic2<-as.data.frame(cbind(mic$MIC,pvalor))
  colnames(mic2)<-c("MIC","p-value")
  return(mic2)
}