#' Dichotomize Data
#'
#' Dichotomization of continuous variables is frequently criticized for its cost of losing information. In this work, we showed that, for variables of bimodal Gaussian mixtures, we can properly preserve more information than commonly used splits by applying our proposed algorithm to obtain optimized cut points.
#'
#' @param X A numeric vector or data frame. Categorical variables and missing data are unacceptable.
#' @export
#' @examples
#' \dontrun{
#'   #simulation data
#'   n1=100;n2=100
#'   m1=0;m2=4;s1=2;s2=1
#'   X1=rnorm(n1,m1,s1)
#'   X2=rnorm(n2,m2,s2)
#'   X=c(X1,X2)
#'
#'   dicho.cor.x(X)
#' }



dicho.cor.x=function(X){
  dicho_output=matrix(data=NA,nrow =1,ncol=6)
  rho=matrix(data=NA,nrow =1,ncol=5)
  N=length(X)
  xx=data.frame(X)

  ########### use mclust ###########
  #mod <- Mclust(xx,verbose=F)
  #mod
  #plot.Mclust(mod, what = "BIC")
  par(mfrow=c(1,2))
  BIC=mclust::mclustBIC(xx)
  #ICL=mclustICL(xx)
  mod2=mclust::Mclust(xx, G = 2, modelNames = "VVE", x=BIC,verbose=FALSE)
  #drmod=MclustDR(mod2, lambda = 1)
  #win.graph()
  #plot(drmod)
  a=summary(mod2,parameters = TRUE)
  m1=as.numeric(a$mean[1])
  dicho_output[,1]=m1
  m2=as.numeric(a$mean[2])
  dicho_output[,2]=m2

  v1=as.numeric(a$variance[1])
  dicho_output[,3]=v1
  s1=sqrt(v1)
  v2=as.numeric(a$variance[2])
  dicho_output[,4]=v2
  s2=sqrt(v2)
  n=plyr::count(a$classification,1)["freq"]
  r=n[1,1]/N
  dicho_output[,5]=r

  a=abs(m2-m1)/(s1+s2)
  dicho_output[,6]=a


  o=function(h){
    l1=r*(1-r)*(m1-m2)^2
    l2=r*(v1)+(1-r)*(v2)
    l3=r*(2*sqrt(v1)*dnorm((h-m1)/sqrt(v1),0,1))
    l4=(1-r)*(2*sqrt(v2)*dnorm((h-m2)/sqrt(v2),0,1))
    l5=r*(1-2*pnorm((h-m1)/sqrt(v1),0,1))
    l6=(1-r)*(1-2*pnorm((h-m2)/sqrt(v2),0,1))
    mx=r*m1+(1-r)*m2

    -(l3+l4+(m1-mx)*l5+(m2-mx)*l6)/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
  }

  #win.graph()
  #plot(o,m1,m2)
  #abline(v=m1,col="red")
  #abline(v=m2,col="blue")
  #o(0)
  opt=optimize(o,c(m1,m2))
  #opt
  ########### h way ###########
  h1=median(X)
  h2=as.numeric(quantile(X,r))
  h3=(s2*m1+s1*m2)/(s1+s2)
  h4=opt$minimum
  h5=mean(X)
  #corr_Xd_X=abs(c$objective)


  ############ dichotomize by h_hat ###########
  Xd.h1=NULL
  Xd.h2=NULL
  Xd.h3=NULL
  Xd.h4=NULL
  Xd.h5=NULL
  for (k in 1:N) {if(X[k]>=h1){Xd.h1[k]=1}else{Xd.h1[k]=-1}}
  for (k in 1:N) {if(X[k]>=h2){Xd.h2[k]=1}else{Xd.h2[k]=-1}}
  for (k in 1:N) {if(X[k]>=h3){Xd.h3[k]=1}else{Xd.h3[k]=-1}}
  for (k in 1:N) {if(X[k]>=h4){Xd.h4[k]=1}else{Xd.h4[k]=-1}}
  for (k in 1:N) {if(X[k]>=h5){Xd.h5[k]=1}else{Xd.h5[k]=-1}}
  rho[,1]=cor(Xd.h1,X)
  rho[,2]=cor(Xd.h2,X)
  rho[,3]=cor(Xd.h3,X)
  rho[,4]=cor(Xd.h4,X)
  rho[,5]=cor(Xd.h5,X)


  rho.cm=colMeans(rho)
  va.cm=colMeans(dicho_output)
  output1=cbind(va.cm[1],va.cm[2],va.cm[3],va.cm[4],va.cm[5],va.cm[6])
  colnames(output1) = c("mean1","mean2","variance1","variance2","gamma","a")

  ht=huxtable::hux(
    method= c('meadian', 'gamma.percentile', '(s2*m1+s1*m2)/(s1+s2)','optimize','mean'),
    h_hat= c(h1,h2,h3,h4,h5),
    corr_Xd.X= c(rho.cm[1],rho.cm[2], rho.cm[3],rho.cm[4],rho.cm[5]),
    add_colnames = TRUE
  )
  huxtable::bottom_border(ht)[1,]  <- 0.4

  #output
  round(output1,4)
  listSample=list(variable=output1,corr_Xd.X=ht)
  listSample
}
