## This R program is used to estimate the performance of 1 se rule on estimating standard error of cross validation error
library(MASS)
library(boot)
library(glmnet)
library(matrixStats)
library(foreach)
library (ModelMetrics)
library(doParallel)
nslots <- Sys.getenv( "SLURM_CPUS_ON_NODE" )
print( nslots )
co<- makeCluster(as.numeric(nslots)/2)
registerDoParallel( cores = co )


setwd(".../")
##1. Generate data set
#parameters settings
n_rep<- 1000
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.1,0.5,0.9) #correlation between
beta.value<- 1.0 #constant beta value
q.sd<- c(10,100,200) #number of non-zero coefficients
qp.sd<- c(0.1,0.5,10/11) #the ratio between the number of non-zero coefficients to total number of coefficients

##generate scenarios for independent variables
scenarios_sd<- expand.grid(q.sd=q.sd,qp.sd=qp.sd,n=n)
##generate scenarios for variables with correlation
scenarios_sd2<- expand.grid(q.sd=q.sd,qp.sd=qp.sd,n=n,rho=rho)

#################################################
#constant coefficient and independent predictors
X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set 

Y.1<- function(Q,Np,N,Beta){
  z<- Beta*X.1(Np,N)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients



#work session
cv.SE.lasso<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.lasso<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SE.ridge<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.ridge<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
SEtoSD.lasso1<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
SEtoSD.ridge1<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
mse.min_lasso<- matrix(0,nrow=4,ncol=n_rep)
mse.min_ridge<- matrix(0,nrow=4,ncol=n_rep)
se.1se_lasso<- matrix(0,nrow=4,ncol=n_rep)
se.1se_ridge<- matrix(0,nrow=4,ncol=n_rep)
a<- list()
fit<- list()


#with constant coefficients but without correlation

for(i in 1:nrow(scenarios_sd)){
q0<- scenarios_sd[i,1]
p0<- q0+scenarios_sd[i,2]
n0<- scenarios_sd[i,3]

a<- foreach(j = 1:n_rep,.packages = 'doSNOW')%dopar%{
  x0<-X.1(p0,n0) 
  y0<- Y.1(q0,p0,n0,1)
  
  foreach(h = 1:4,.packages='glmnet')%do%{
    fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
    fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
    
    i.min_lasso <- which(fit.lasso$lambda == fit.lasso$lambda.min)
    fit$lmin<- fit.lasso$cvm[i.min_lasso]
    
    i.min_ridge <- which(fit.ridge$lambda == fit.ridge$lambda.min)
    fit$rmin<- fit.ridge$cvm[i.min_ridge]
    
    fit$lse<- fit.lasso$cvsd[i.min_lasso]
    
    fit$rse<- fit.ridge$cvsd[i.min_ridge]
    return(fit)
  }
}
for(j in 1:n_rep){
  for(h in 1:4){
    mse.min_lasso[h,j]<-a[[j]][[h]]$lmin
    mse.min_ridge[h,j]<-a[[j]][[h]]$rmin
    se.1se_lasso[h,j]<-a[[j]][[h]]$lse
    se.1se_ridge[h,j]<-a[[j]][[h]]$rse
  }
}

cv.SE.lasso[,i]<- rowMeans(se.1se_lasso,na.rm=TRUE)
cv.SE.ridge[,i]<- rowMeans(se.1se_ridge,na.rm=TRUE)
cv.SD.lasso[,i]<- rowSds(mse.min_lasso,na.rm=TRUE)
cv.SD.ridge[,i]<- rowSds(mse.min_ridge,na.rm=TRUE)
SEtoSD.lasso[,i]<- cv.SE.lasso[,i]/cv.SD.lasso[,i]
SEtoSD.ridge[,i]<- cv.SE.ridge[,i]/cv.SD.ridge[,i]
}
write.matrix(SEtoSD.lasso1,file="setosd1_lasso.csv")
write.matrix(SEtoSD.ridge1,file="setosd1_ridge.csv")


###################################################################
#decaying coefficient and independent predictors
X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}
Y.2<- function(Q,Np,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X.1(Np,N)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with decaying coefficients without correlation


#with decaying coefficients but without correlation
cv.SE.lasso2<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.lasso2<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SE.ridge2<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.ridge2<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
SEtoSD.lasso2<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
SEtoSD.ridge2<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
mse.min_lasso2<- matrix(0,nrow=4,ncol=n_rep)
mse.min_ridge2<- matrix(0,nrow=4,ncol=n_rep)
se.1se_lasso2<- matrix(0,nrow=4,ncol=n_rep)
se.1se_ridge2<- matrix(0,nrow=4,ncol=n_rep)

for (i in 1:nrow(scenarios_sd)){
  q0<- scenarios_sd[i,1]
  p0<- q0/scenarios_sd[i,2]
  n0<- scenarios_sd[i,3]
  
  for (j in 1:n_rep){
    x0<-X.1(p0,n0)
    y0<- Y.2(q0,p0,n0,1)
    
    for (h in 1:4){
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
      fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
      
      i.min_lasso <- which(fit.lasso$lambda == fit.lasso$lambda.min)
      mse.min_lasso2[h,j] <- fit.lasso$cvm[i.min_lasso]
      
      i.min_ridge <- which(fit.ridge$lambda == fit.ridge$lambda.min)
      mse.min_ridge2[h,j] <- fit.ridge$cvm[i.min_ridge]
      
      se.1se_lasso2[h,j] <- fit.lasso$cvsd[i.min_lasso]
      
      se.1se_ridge2[h,j] <- fit.ridge$cvsd[i.min_ridge]
    }
  }
  cv.SE.lasso2[,i]<- rowMeans(se.1se_lasso2)
  cv.SE.ridge2[,i]<- rowMeans(se.1se_ridge2)
  cv.SD.lasso2[,i]<- rowSds(mse.min_lasso2)
  cv.SD.ridge2[,i]<- rowSds(mse.min_ridge2)
  SEtoSD.lasso2[,i]<- cv.SE.lasso2[,i]/cv.SD.lasso2[,i]
  SEtoSD.ridge2[,i]<- cv.SE.ridge2[,i]/cv.SD.ridge2[,i]
}
write.matrix(SEtoSD.lasso2,file="setosd2_lasso.csv")
write.matrix(SEtoSD.ridge2,file="setosd2_ridge.csv")

############################################
#constant coefficient and compound symmetry correlation structure
X.3<- function(Np,N,RHO){
  CS<- CSgenerate(Np,RHO)
  mvrnorm(N,rep(1,Np),CS)
}#This is to generate train set with compound symmetry covariance
Y.5<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients and compound symmetry correlation structure


#with constant coefficients and compound simmetry covariance structure
cv.SE.lasso5<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SD.lasso5<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SE.ridge5<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SD.ridge5<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
SEtoSD.lasso5<- matrix(0,ncol= nrow(scenarios_sd2),nrow=4)
SEtoSD.ridge5<- matrix(0,ncol= nrow(scenarios_sd2),nrow=4)
mse.min_lasso5<- matrix(0,nrow=4,ncol=n_rep)
mse.min_ridge5<- matrix(0,nrow=4,ncol=n_rep)
se.1se_lasso5<- matrix(0,nrow=4,ncol=n_rep)
se.1se_ridge5<- matrix(0,nrow=4,ncol=n_rep)

for (i in 1:nrow(scenario_sd2)){
  q0<- scenarios_sd2[i,1]
  p0<- q0/scenarios_sd2[i,2]
  n0<- scenarios_sd2[i,3]
  rho0<- scenarios_sd2[i,4]
  a<- foreach(j = 1:n_rep,.packages = c('doSNOW','MixMatrix','MASS'))%dopar%{
    x0<-X.3(p0,n0,rho0)
    y0<- Y.5(q0,x0,n0,1)
    
  foreach(h = 1:4,.packages='glmnet')%do%{
    fitting <- FALSE
    
    while(!fitting) {
      
      tmp <- tryCatch({
        
        fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
        fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
        fitting <- TRUE
        
      },
      error = function(e) {
        x0<-X.3(p0,n0,rho0)
        y0<- Y.5(q0,x0,n0,1)
      })
    }
    
    i.min_lasso <- which(fit.lasso$lambda == fit.lasso$lambda.min)
    fit$lmin<- fit.lasso$cvm[i.min_lasso]
    
    i.min_ridge <- which(fit.ridge$lambda == fit.ridge$lambda.min)
    fit$rmin<- fit.ridge$cvm[i.min_ridge]
    
    fit$lse<- fit.lasso$cvsd[i.min_lasso]
    
    fit$rse<- fit.ridge$cvsd[i.min_ridge]
    return(fit)
  }
}

for(j in 1:n_rep){
  for(h in 1:4){
    mse.min_lasso5[h,j]<-a[[j]][[h]]$lmin
    mse.min_ridge5[h,j]<-a[[j]][[h]]$rmin
    se.1se_lasso5[h,j]<-a[[j]][[h]]$lse
    se.1se_ridge5[h,j]<-a[[j]][[h]]$rse
  }
}
  cv.SE.lasso5[,i]<- rowMeans(se.1se_lasso5)
  cv.SE.ridge5[,i]<- rowMeans(se.1se_ridge5)
  cv.SD.lasso5[,i]<- rowSds(mse.min_lasso5)
  cv.SD.ridge5[,i]<- rowSds(mse.min_ridge5)
  SDtoSE.lasso5[,i]<- cv.SE.lasso5[,i]/cv.SD.lasso5[,i]
  SDtoSE.ridge5[,i]<- cv.SE.ridge5[,i]/cv.SD.ridge5[,i]
}
write.matrix(SEtoSD.lasso5,file="setosd3_lasso.csv")
write.matrix(SEtoSD.ridge5,file="setosd3_ridge.csv")


#####################################################
#decaying coefficient and compound symmetry correlation structure
X.3<- function(Np,N,RHO){
  CS<- CSgenerate(Np,RHO)
  mvrnorm(N,rep(1,Np),CS)
}#This is to generate train set with compound symmetry covariance
Y.6<- function(Q,Np,N,RHO,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X.3(Np,N,RHO)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with decaying coefficients and compound symmetry correlation


#with decaying coefficients and compound simmetry covariance structure
cv.SE.lasso6<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SD.lasso6<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SE.ridge6<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
cv.SD.ridge6<- matrix(0,ncol = nrow(scenarios_sd2),nrow=4)
SEtoSD.lasso6<- matrix(0,ncol= nrow(scenarios_sd2),nrow=4)
SEtoSD.ridge6<- matrix(0,ncol= nrow(scenarios_sd2),nrow=4)
mse.min_lasso6<- matrix(0,nrow=4,ncol=n_rep)
mse.min_ridge6<- matrix(0,nrow=4,ncol=n_rep)
se.1se_lasso6<- matrix(0,nrow=4,ncol=n_rep)
se.1se_ridge6<- matrix(0,nrow=4,ncol=n_rep)
a<- list()
fit<- list()

for (i in 1:nrow(scenarios_sd2)){
  q0<- scenarios_sd2[i,1]
  p0<- q0/scenarios_sd2[i,2]
  n0<- scenarios_sd2[i,3]
  rho0<- scenarios_sd2[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = 'doSNOW')%dopar%{
    x0<-X.3(p0,n0,rho0)
    y0<- Y.6(q0,x0,n0,rho0,1)
    
    foreach(h = 1:4,.packages='glmnet')%do%{
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
          fitting <- TRUE
          
        },
        error = function(e) {
          x0<-X.3(p0,n0,rho0)
          y0<- Y.6(q0,x0,n0,rho0,1)
        })
      }
      
      i.min_lasso <- which(fit.lasso$lambda == fit.lasso$lambda.min)
      fit$lmin<- fit.lasso$cvm[i.min_lasso]
      
      i.min_ridge <- which(fit.ridge$lambda == fit.ridge$lambda.min)
      fit$rmin<- fit.ridge$cvm[i.min_ridge]
      
      fit$lse<- fit.lasso$cvsd[i.min_lasso]
      
      fit$rse<- fit.ridge$cvsd[i.min_ridge]
      return(fit)
    }
  }
  
  for(j in 1:n_rep){
    for(h in 1:4){
      mse.min_lasso6[h,j]<-a[[j]][[h]]$lmin
      mse.min_ridge6[h,j]<-a[[j]][[h]]$rmin
      se.1se_lasso6[h,j]<-a[[j]][[h]]$lse
      se.1se_ridge6[h,j]<-a[[j]][[h]]$rse
    }
  }
  cv.SE.lasso6[,i]<- rowMeans(se.1se_lasso6)
  cv.SE.ridge6[,i]<- rowMeans(se.1se_ridge6)
  cv.SD.lasso6[,i]<- rowSds(mse.min_lasso6)
  cv.SD.ridge6[,i]<- rowSds(mse.min_ridge6)
  SEtoSD.lasso6[,i]<- cv.SE.lasso6[,i]/cv.SD.lasso6[,i]
  SEtoSD.ridge6[,i]<- cv.SE.ridge6[,i]/cv.SD.ridge6[,i]
}
write.matrix(SEtoSD.lasso6,file="setosd4_lasso.csv")
write.matrix(SEtoSD.ridge6,file="setosd4_ridge.csv")
