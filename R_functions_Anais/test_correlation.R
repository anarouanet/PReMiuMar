

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



old=1

if(old==1){
  times=seq(50,100,by=1)
  sigma_b<-c(0,0.1)
  sigma_v<-c(0.05,0.11)
  l_lin=log(65)

  sigma_b<-c(0,0)
  sigma_v<-c(-3,-2)
  l_lin <-c(log(65),log(80))
  sigma_e <- c(-4, -4)

  #linear
  sigma_b<-c(0,0)
  sigma_v<-c(-10,-10)
  l_lin <-c(log(65),log(80))
  sigma_e <- c(-4, -4)
}else{
  times=seq(0,12,by=1)
  sigma_b<-c(0,0.1^2)
  sigma_v<-c(0.05^2,0.11^2)
  l_lin=5

  sigma_b<-c(0,0)
  sigma_v<-c(-3,-2)
  l_lin <-c(1.09,2.07)
  sigma_e <- c(-4, -4)
}



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
y_quad2=mvrnorm(100,mu,Sigma_quad2)
y_quad3=mvrnorm(100,mu,Sigma_quad3)

#plot(y~times,type='l',ylim=c(min(y),max(y)))
par(mfrow=c(1,1))
plot(y_quad2~times,type='l',ylim=c(min(c(y_quad2,y_quad3)),max(c(y_quad2,y_quad3))))
lines(y_quad3~times)

a=rlnorm(100, meanlog = 0, sdlog = 1)
hist(a)
mean(a)

lin=0
if(lin==1){
  par(mfrow = c(1,3))
  plot(y_exp[1,]~times,type='l',ylim=c(min(y_exp),max(y_exp)))
  for(i in 1:10)
    lines(y_exp[i,]~times)

  plot(y_lin[1,]~times,type='l',ylim=c(min(y_lin),max(y_lin)))
  for(i in 1:10)
    lines(y_lin[i,]~times)


  plot(y_quad[1,]~times,type='l',ylim=c(min(y_quad),max(y_quad)))
  for(i in 1:10)
    lines(y_quad[i,]~times)

}else{
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
}



dat<-y_quad[,1]
for(i in 1:dim(y_quad)[1])
  dat<- c(dat,y_quad[,i])

data<-matrix(0,nrow=length(dat),ncol=2,dimnames=c('x','y'))
colnames(data)<-c("x","y")
data$y<-t(dat)
data$x<-rep(times,dim(y_quad)[1])

head(data)

m <- gp(y_quad[,1] ~ rbf('x'), data = data, family = gaussian)

# predict from it
pred_df <- data.frame(x = seq(min(df$x), max(df$x), len = 500))
lambda <- predict(m, pred_df, type = 'response')

# plot the predicted rate parameter, the true model and the data
plot(lambda ~ pred_df$x, type = 'l', lwd = 2, ylim = range(y))
lines(exp(f(pred_df$x)) ~ pred_df$x, lty = 2)
points(y ~ x, data = df)
