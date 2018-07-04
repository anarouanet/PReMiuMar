kernel='Quadratic'
LArray<- riskProfileObj$LArray
LMeans <- matrix(0,ncol=ifelse(kernel=="SQexponential",3,4),nrow=dim(LArray)[2])
#LMean <- c()

for(i in 1:ifelse(kernel=="SQexponential",3,4)){
  LMeans[,i]<-exp(apply(LArray[,,i],2,mean,trim=0.005))
#  LMean[i]<-sum(LMeans*clusterSizes)/sum(clusterSizes)
}

times=seq(min(riskProfileObj$riskProfClusObj$clusObjRunInfoObj$longMat$time),max(riskProfileObj$riskProfClusObj$clusObjRunInfoObj$longMat$time),by=0.1)
nT=length(times)
mu=rep(0,nT)
Sigma_quad1=matrix(0,nT,nT)
Sigma_quad2=matrix(0,nT,nT)
Sigma_quad3=matrix(0,nT,nT)


g=2

for(i in 1:nT){
  for(j in 1:nT){
    Sigma_quad1[i,j]= (LMeans[g-1,1]+LMeans[g-1,2]*(times[i]-LMeans[g-1,4])*(times[j]-LMeans[g-1,4]))^2+ifelse(i!=j,0,LMeans[g-1,3])
    Sigma_quad2[i,j]= (LMeans[g,1]+LMeans[g,2]*(times[i]-LMeans[g,4])*(times[j]-LMeans[g,4]))^2+ifelse(i!=j,0,LMeans[g,3])
    Sigma_quad3[i,j]= (LMeans[g+1,1]+LMeans[g+1,2]*(times[i]-LMeans[g+1,4])*(times[j]-LMeans[g+1,4]))^2+ifelse(i!=j,0,LMeans[g+1,3])
  }
}

library(MASS)
y_quad1=mvrnorm(100,mu,Sigma_quad1)
y_quad2=mvrnorm(100,mu,Sigma_quad2)
y_quad3=mvrnorm(100,mu,Sigma_quad3)

#plot(y~times,type='l',ylim=c(min(y),max(y)))
par(mfrow = c(1,3))
plot(y_quad1[1,]~times,type='l',ylim=c(min(y_quad1),max(y_quad1)),ylab=c(sigma_b[1],sigma_v[1]))
for(i in 1:10)
  lines(y_quad1[i,]~times)

plot(y_quad2[1,]~times,type='l',ylim=c(min(y_quad1),max(y_quad1)),ylab=c(sigma_b[1],sigma_v[2]))
for(i in 1:10)
  lines(y_quad2[i,]~times)


plot(y_quad3[1,]~times,type='l',ylim=c(min(y_quad1),max(y_quad1)),ylab=c(sigma_b[2],sigma_v[1]))
for(i in 1:10)
  lines(y_quad3[i,]~times)
sigma_b<-c(0.05,0.5)
sigma_v<-c(-3,-2)
l_lin <-c(1.4,2.07)
sigma_e <- c(-4, -4)

sigma_b<-exp(sigma_b)
sigma_v<-exp(sigma_v)
l_lin <-exp(l_lin)
sigma_e <-exp(sigma_e)

nT=length(times)
mu=rep(0,nT)
Sigma_exp=matrix(0,nT,nT)
Sigma_lin=matrix(0,nT,nT)
Sigma_quad1=matrix(0,nT,nT)
Sigma_quad2=matrix(0,nT,nT)
Sigma_quad3=matrix(0,nT,nT)
library(MASS)
lexp=35
l2=20


for(i in 1:nT){
  for(j in 1:nT){
    Sigma_exp[i,j]=l2*exp(-(times[i]-times[j])^2/(2*lexp))
    Sigma_lin[i,j]=sigma_b[1]^2+sigma_v[1]^2*(times[i]-l_lin[1])*(times[j]-l_lin[1])
    Sigma_quad1[i,j]= (sigma_b[1]+sigma_v[1]*(times[i]-l_lin[1])*(times[j]-l_lin[1]))^2+ifelse(i!=j,0,sigma_e[1])
    Sigma_quad2[i,j]= (sigma_b[1]+sigma_v[1]*(times[i]-l_lin[1])*(times[j]-l_lin[1]))^2+ifelse(i!=j,0,sigma_e[1])
    Sigma_quad3[i,j]= (sigma_b[2]+sigma_v[2]*(times[i]-l_lin[2])*(times[j]-l_lin[2]))^2+ifelse(i!=j,0,sigma_e[1])
  }
}

y_lin=mvrnorm(100,mu,Sigma_lin)
y_exp=mvrnorm(100,mu,Sigma_exp)
y_quad1=mvrnorm(100,mu,Sigma_quad1)
y_quad2=mvrnorm(1,mu,Sigma_quad2)
y_quad3=mvrnorm(1,mu,Sigma_quad3)

#plot(y~times,type='l',ylim=c(min(y),max(y)))
par(mfrow=c(1,1))
plot(y_quad2~times,type='l',ylim=c(min(c(y_quad2,y_quad3)),max(c(y_quad2,y_quad3))))
lines(y_quad3~times,col=2)
