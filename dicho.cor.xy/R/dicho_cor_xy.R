#' Dichotomize Data
#'
#' Dichotomization of continuous variables is frequently criticized for its cost of losing information. In this work, we showed that, for variables of bimodal Gaussian mixtures, we can properly preserve more information than commonly used splits by applying our proposed algorithm to obtain optimized cut points.
#'
#' @param  X A numeric vector or data frame. Categorical variables and missing data are unacceptable.
#' @param  Y A numeric vector or data frame. Categorical variables and missing data are unacceptable.
#' @export
#' @examples
#' \dontrun{
#'   #simulation data
#'   n1=100;n2=100
#'   m1=0;m2=4;s1=2;s2=1
#'   X1=rnorm(n1,m1,s1)
#'   X2=rnorm(n2,m2,s2)
#'   X=c(X1,X2)
#'   B1=0.1;B2=0.15
#'   eps=rnorm(n1+n2,sd=2)
#'   Y1=B1*X1+eps[1:n1]
#'   Y2=B2*X2+eps[n1+1:n2]
#'   Y=c(Y1,Y2)
#'
#'   dicho.cor.xy(X,Y)
#' }

dicho.cor.xy=function(X,Y){
  dicho_output=matrix(data=NA,nrow =10,ncol=8)
  rho_Xd.X=matrix(data=NA,nrow =1,ncol=5)
  rho_Xd.Y=matrix(data=NA,nrow =1,ncol=5)
  N=length(X)
  xx=data.frame(X)
  ##########piecewise#######
  fit=lm(Y ~ X)
  segmented.fit=segmented::segmented(fit, seg.Z = ~X, psi=mean(X))
  piecewise=summary(segmented.fit)
  B1=piecewise$coefficients[2,1]
  B2=piecewise$coefficients[3,1]+piecewise$coefficients[2,1]


  ########### use mclust ###########
  #mod <- Mclust(xx)
  #mod
  #plot.Mclust(mod, what = "BIC")
  par(mfrow=c(1,2))
  BIC=mclust::mclustBIC(xx)
  #ICL=mclustICL(xx)
  mod2=mclust::Mclust(xx, G = 2, modelNames = "VVE", x=BIC,verbose=F)
  #drmod=MclustDR(mod2, lambda = 1)
  #win.graph()
  #plot(drmod)
  a=summary(mod2,parameters = T)
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
    C1=(l3+m1*((1-r)*l5-r*l6))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    C2=(l4+m2*(r*l6-(1-r)*l5))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    l1_=r*(1-r)*(B1*m1-B2*m2)^2
    l2_=r*B1^2*v1+(1-r)*B2^2*v2
    (B1*C1+B2*C2)*sqrt(l1+l2)/sqrt(l1_+l2_+1^2)
  }
  opt=optimize(o,c(m1,m2))
  if(opt$objective<0){o=function(h){
    l1=r*(1-r)*(m1-m2)^2
    l2=r*(v1)+(1-r)*(v2)
    l3=r*(2*sqrt(v1)*dnorm((h-m1)/sqrt(v1),0,1))
    l4=(1-r)*(2*sqrt(v2)*dnorm((h-m2)/sqrt(v2),0,1))
    l5=r*(1-2*pnorm((h-m1)/sqrt(v1),0,1))
    l6=(1-r)*(1-2*pnorm((h-m2)/sqrt(v2),0,1))
    mx=r*m1+(1-r)*m2
    C1=(l3+m1*((1-r)*l5-r*l6))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    C2=(l4+m2*(r*l6-(1-r)*l5))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    l1_=r*(1-r)*(B1*m1-B2*m2)^2
    l2_=r*B1^2*v1+(1-r)*B2^2*v2
    (B1*C1+B2*C2)*sqrt(l1+l2)/sqrt(l1_+l2_+1^2)
  }}else{o=function(h){
    l1=r*(1-r)*(m1-m2)^2
    l2=r*(v1)+(1-r)*(v2)
    l3=r*(2*sqrt(v1)*dnorm((h-m1)/sqrt(v1),0,1))
    l4=(1-r)*(2*sqrt(v2)*dnorm((h-m2)/sqrt(v2),0,1))
    l5=r*(1-2*pnorm((h-m1)/sqrt(v1),0,1))
    l6=(1-r)*(1-2*pnorm((h-m2)/sqrt(v2),0,1))
    mx=r*m1+(1-r)*m2
    C1=(l3+m1*((1-r)*l5-r*l6))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    C2=(l4+m2*(r*l6-(1-r)*l5))/(sqrt(1-(l5+l6)^2)*sqrt(l1+l2))
    l1_=r*(1-r)*(B1*m1-B2*m2)^2
    l2_=r*B1^2*v1+(1-r)*B2^2*v2
    -(B1*C1+B2*C2)*sqrt(l1+l2)/sqrt(l1_+l2_+1^2)
  }}
  opt=optimize(o,c(m1,m2))
  ########### h way ###########
  h1=median(X)
  h2=as.numeric(quantile(X,r))
  h3=(s2*m1+s1*m2)/(s1+s2)
  h4=opt$minimum
  h5=mean(X)

  ############ dichotomize by h ###########
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
  rho_Xd.X[,1]=cor(Xd.h1,X)
  rho_Xd.X[,2]=cor(Xd.h2,X)
  rho_Xd.X[,3]=cor(Xd.h3,X)
  rho_Xd.X[,4]=cor(Xd.h4,X)
  rho_Xd.X[,5]=cor(Xd.h5,X)

  rho_Xd.Y[,1]=cor(Xd.h1,Y)
  rho_Xd.Y[,2]=cor(Xd.h2,Y)
  rho_Xd.Y[,3]=cor(Xd.h3,Y)
  rho_Xd.Y[,4]=cor(Xd.h4,Y)
  rho_Xd.Y[,5]=cor(Xd.h5,Y)

  rho.x.cm=colMeans(rho_Xd.X)
  rho.y.cm=colMeans(rho_Xd.Y)
  va.cm=colMeans(dicho_output)
  output1=cbind(va.cm[1],va.cm[2],va.cm[3],va.cm[4],va.cm[5],va.cm[6],B1,B2)
  colnames(output1) = c("mean1","mean2","variance1","variance2","gamma","a","beta1","beta2")

  ht=huxtable::hux(
    method= c('meadian', 'gamma.percentile', '(s2*m1+s1*m2)/(s1+s2)','optimize','mean'),
    h_hat= c(h1,h2,h3,h4,h5),
    corr_Xd.X= c(rho.x.cm[1],rho.x.cm[2], rho.x.cm[3],rho.x.cm[4],rho.x.cm[5]),
    add_colnames = TRUE
  )
  huxtable::bottom_border(ht)[1,]  <- 0.4

  ht_=huxtable::hux(
    method= c('meadian', 'gamma.percentile', '(s2*m1+s1*m2)/(s1+s2)','optimize','mean'),
    h_hat= c(h1,h2,h3,h4,h5),
    corr_Xd.Y= c(rho.y.cm[1],rho.y.cm[2],rho.y.cm[3],rho.y.cm[4],rho.y.cm[5]),
    add_colnames = TRUE
  )
  huxtable::bottom_border(ht_)[1,]  <- 0.4

  #output
  round(output1,4)
  listSample=list(variable=output1,corr_Xd.X=ht,corr_Xd.Y=ht_)
  listSample
}
