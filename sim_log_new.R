## This R program is used to estimate the performance of 1 se rule in cross validation when tuning logistic model.
##Both Lasso and Ridge would be tested. The estimation would contain 3 main parts: i) standard error of cross validation error;
## ii) accuracy of variable selection; iii) prediction error.
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
registerDoParallel( cores = 6 )


setwd("C:/Users/joeyt/Desktop/research_project/simulation")
##1. Generate data set
#parameters settings
n_rep<- 1000
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.01,0.1,0.5,0.9) #correlation between
beta.value<- c(1.0,1.5,2.0,3.0) #constant beta value

X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set 
X.1.test<- function(Np,N){
  matrix(rnorm(N*Np*0.3),nrow=0.3*N,ncol=Np)
} 
Y.1<- function(Q,Np,N,Beta){
  z<- Beta*X.1(Np,N)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients
Y.1.test<- function(Q,Np,N,Beta){
  z<- Beta*X.1.test(Np,N)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(0.3*N,1,pr)
}

##Maybe to set specific data sets for each part is a better idea. The data set I set above may not make sense to every part of our test. 
##The ratio between parameters may not be comparable.

##true standard error for cross validation error vs.  standard error ignoring the correlation among groups
#In this part, we would test under following conditions:i) K={4,5,10,20};ii)n={100,1000};iii)sigma={0.1,1,2};iv)rho={0,0.1,0.5,0.9};
#v)the value of constant beta={0.1,0.5,1,2};vi)the correlation structure is AR(1) and compound symmetry covariance.


q.sd<- c(5) #the value of q for standard deviation test
zq.sd<- c(1,2,3,4)#the ratio between the number of non-zero coefficients to total number of coefficients


scenarios_sd<- expand.grid(q.sd=q.sd,zq.sd=zq.sd,n=n,beta.value=beta.value)

#work session
cv.SE.lasso<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.lasso<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SE.ridge<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
cv.SD.ridge<- matrix(0,ncol = nrow(scenarios_sd),nrow=4)
SEtoSD.lasso<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
SEtoSD.ridge<- matrix(0,ncol= nrow(scenarios_sd),nrow=4)
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
beta0<- scenarios_sd[i,4]

a<- foreach(j = 1:n_rep,.packages = 'doSNOW')%dopar%{
  x0<-X.1(p0,n0) 
  y0<- Y.1(q0,p0,n0,beta0)
  
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

SEtoSD.lasso
SEtoSD.ridge

write.matrix(SEtoSD.lasso,file="sd_new_lasso.csv")
write.matrix(SEtoSD.ridge,file="sd_new_ridge.csv")
