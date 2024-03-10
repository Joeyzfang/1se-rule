## This R program is used to estimate the performance of 1 se rule on estimating variable selection.
library(boot)
library(glmnet)
library(matrixStats)
library(MixMatrix)
library(foreach)
library (ModelMetrics)
library(doParallel)
registerDoParallel(cores= )
##1. Generate data set
#parameters settings
n_rep<- 1000
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.1,0.9) #correlation between
##Variable selection (only for lasso)
q.vs<- c(10,100,200)    #This is the vector of potential value of q for variable selection
qp.vs<- c(0.1,0.5,10/11) #In this section, qp.vs refers to the difference between q and p which is equal to p - q.

scenarios_vs<- expand.grid(q.vs=q.vs,qp.vs=qp.vs,n=n,beta.value=beta.value)
scenarios_vs2<- expand.grid(q.vs=q.vs,qp.vs=qp.vs,n=n,beta.value=beta.value,rho=rho)
###############################################################################
X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set

Y.1<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients


fit<- list()
a<- list()

#with constant coefficients but without correlation

trco.min<- matrix(0,nrow=4,ncol=n_rep)
trco.1se<- matrix(0,nrow=4,ncol=n_rep)
faco.min<- matrix(0,nrow=4,ncol=n_rep)
faco.1se<- matrix(0,nrow=4,ncol=n_rep)
trze.min<- matrix(0,nrow=4,ncol=n_rep)
trze.1se<- matrix(0,nrow=4,ncol=n_rep)
faze.min<- matrix(0,nrow=4,ncol=n_rep)
faze.1se<- matrix(0,nrow=4,ncol=n_rep)

for(i in 1:nrow(scenarios_vs)){
  q0<- scenarios_vs[i,1]
  p0<- q0/scenarios_vs[i,2]
  n0<- scenarios_vs[i,3]
  beta0<- scenarios_vs[i,4]
  a<- foreach (j = 1:n_rep,.packages='doSNOW')%dopar%{
    x0<-X.1(p0,n0)
    y0<- Y.1(q0,x0,n0,beta0)
    
    foreach (h = 1:4,.packages='glmnet')%do%{
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fitting <- TRUE
          
        },error = function(e) {
          x0<-X.1(p0,n0)
          y0<- Y.1(q0,x0,n0,beta0)
        })
      }
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1)
      coef.min<- coef(fit.lasso,fit.lasso$lambda.min)
      coef.1se<- coef(fit.lasso,fit.lasso$lambda.1se)
      coef.min<- coef.min[-1,]
      coef.1se<- coef.1se[-1,]
      true.correct_min<- length(which(which(coef.min!=0)< (q0+1)))
      true.correct_1se<- length(which(which(coef.1se!=0)< (q0+1)))
      fake.correct_min<- length(which(coef.min!=0))-true.correct_min
      fake.correct_1se<- length(which(coef.1se!=0))-true.correct_1se
      fake.zero_min<- q0-true.correct_min
      fake.zero_1se<- q0-true.correct_1se
      true.zero_min<- length(which(coef.min==0))-fake.zero_min
      true.zero_1se<- length(which(coef.1se==0))-fake.zero_1se
      fit$trco.min<- true.correct_min/q0
      fit$trco.1se<- true.correct_1se/q0
      fit$faco.min<- fake.correct_min/q0
      fit$faco.1se<- fake.correct_1se/q0
      fit$trze.min<- true.zero_min/q0
      fit$trze.1se<- true.zero_1se/q0
      fit$faze.min<- fake.zero_min/q0
      fit$faze.1se<- fake.zero_1se/q0
      
      
      return(fit)
    }
  }
  for (j in 1:n_rep) {
    for(h in 1:4){
      trco.min[h,j]<-a[[j]][[h]]$trco.min
      trco.1se[h,j]<-a[[j]][[h]]$trco.1se
      faco.min[h,j]<-a[[j]][[h]]$faco.min
      faco.1se[h,j]<-a[[j]][[h]]$faco.1se
      trze.min[h,j]<-a[[j]][[h]]$trze.min
      trze.1se[h,j]<-a[[j]][[h]]$trze.1se
      faze.min[h,j]<-a[[j]][[h]]$faze.min
      faze.1se[h,j]<-a[[j]][[h]]$faze.1se
    }
  }
  write.matrix(trco.min,file=paste0("vs.trco.min1_",i,".csv"))
  write.matrix(trco.1se,file=paste0("vs.trco.1se1_",i,".csv"))
  write.matrix(faco.min,file=paste0("vs.faco.min1_",i,".csv"))
  write.matrix(faco.1se,file=paste0("vs.faco.1se1_",i,".csv"))
  write.matrix(trze.min,file=paste0("vs.trze.min1_",i,".csv"))
  write.matrix(trze.1se,file=paste0("vs.trze.1se1_",i,".csv"))
  write.matrix(faze.min,file=paste0("vs.faze.min1_",i,".csv"))
  write.matrix(faze.1se,file=paste0("vs.faze.1se1_",i,".csv"))
}


#####################################################################
X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set
Y.2<- function(Q,Np,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X.1(Np,N)[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with decaying coefficients without correlation


fit<- list()
a<- list()

#with decaying coefficients but without correlation
trco.min<- matrix(0,nrow=4,ncol=n_rep)
trco.1se<- matrix(0,nrow=4,ncol=n_rep)
faco.min<- matrix(0,nrow=4,ncol=n_rep)
faco.1se<- matrix(0,nrow=4,ncol=n_rep)
trze.min<- matrix(0,nrow=4,ncol=n_rep)
trze.1se<- matrix(0,nrow=4,ncol=n_rep)
faze.min<- matrix(0,nrow=4,ncol=n_rep)
faze.1se<- matrix(0,nrow=4,ncol=n_rep)


for (i in 1:nrow(scenarios_vs)){
  q0<- scenarios_vs[i,1]
  p0<- q0+scenarios_vs[i,2]
  n0<- scenarios_vs[i,3]
  beta0<- scenarios_vs[i,4]
  a<- foreach (j = 1:n_rep,.packages='doSNOW')%dopar%{
    x0<-X.1(p0,n0)
    y0<- Y.2(q0,p0,n0,beta0)
    
    foreach (h = 1:4,.packages='glmnet')%do%{
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1)
      coef.min<- coef(fit.lasso,fit.lasso$lambda.min)
      coef.1se<- coef(fit.lasso,fit.lasso$lambda.1se)
      coef.min<- coef.min[-1,]
      coef.1se<- coef.1se[-1,]
      true.correct_min<- length(which(which(coef.min!=0)< (q0+1)))
      true.correct_1se<- length(which(which(coef.1se!=0)< (q0+1)))
      fake.correct_min<- length(which(coef.min!=0))-true.correct_min
      fake.correct_1se<- length(which(coef.1se!=0))-true.correct_1se
      fake.zero_min<- q0-true.correct_min
      fake.zero_1se<- q0-true.correct_1se
      true.zero_min<- length(which(coef.min==0))-fake.zero_min
      true.zero_1se<- length(which(coef.1se==0))-fake.zero_1se
      fit$trco.min<- true.correct_min/q0
      fit$trco.1se<- true.correct_1se/q0
      fit$faco.min<- fake.correct_min/q0
      fit$faco.1se<- fake.correct_1se/q0
      fit$trze.min<- true.zero_min/q0
      fit$trze.1se<- true.zero_1se/q0
      fit$faze.min<- fake.zero_min/q0
      fit$faze.1se<- fake.zero_1se/q0
      
      return(fit)
    }
  }
  for (j in 1:n_rep) {
    for(h in 1:4){
      trco.min[h,j]<-a[[j]][[h]]$trco.min
      trco.1se[h,j]<-a[[j]][[h]]$trco.1se
      faco.min[h,j]<-a[[j]][[h]]$faco.min
      faco.1se[h,j]<-a[[j]][[h]]$faco.1se
      trze.min[h,j]<-a[[j]][[h]]$trze.min
      trze.1se[h,j]<-a[[j]][[h]]$trze.1se
      faze.min[h,j]<-a[[j]][[h]]$faze.min
      faze.1se[h,j]<-a[[j]][[h]]$faze.1se
    }
  }
  write.matrix(trco.min,file=paste0("vs.trco.min2_",i,".csv"))
  write.matrix(trco.1se,file=paste0("vs.trco.1se2_",i,".csv"))
  write.matrix(faco.min,file=paste0("vs.faco.min2_",i,".csv"))
  write.matrix(faco.1se,file=paste0("vs.faco.1se2_",i,".csv"))
  write.matrix(trze.min,file=paste0("vs.trze.min2_",i,".csv"))
  write.matrix(trze.1se,file=paste0("vs.trze.1se2_",i,".csv"))
  write.matrix(faze.min,file=paste0("vs.faze.min2_",i,".csv"))
  write.matrix(faze.1se,file=paste0("vs.faze.1se2_",i,".csv"))
}


#######################################################################################
X.3<- function(Np,N,RHO){
  CS<- CSgenerate(Np,RHO)
  mvrnorm(N,rep(1,Np),CS)
}#This is to generate train set with compound symmetry covariance
Y.5<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients and compound symmetry correlation structure

q.vs<- c(10,100,200)    #This is the vector of potential value of q for variable selection
qp.vs<- c(0.1,0.5,10/11) #In this section, qp.vs refers to the difference between q and p which is equal to p - q.
#generate scenarios for variables with correlation 
scenarios_vs2<- expand.grid(q.vs=q.vs,qp.vs=qp.vs,n=n,beta.value=beta.value,rho=rho)

fit<- list()
a<- list()

#with constant coefficients and compound symmetry correlation

trco.min<- matrix(0,nrow=4,ncol=n_rep)
trco.1se<- matrix(0,nrow=4,ncol=n_rep)
faco.min<- matrix(0,nrow=4,ncol=n_rep)
faco.1se<- matrix(0,nrow=4,ncol=n_rep)
trze.min<- matrix(0,nrow=4,ncol=n_rep)
trze.1se<- matrix(0,nrow=4,ncol=n_rep)
faze.min<- matrix(0,nrow=4,ncol=n_rep)
faze.1se<- matrix(0,nrow=4,ncol=n_rep)


for (i in 1:nrow(scenarios_vs2)){
  q0<- scenarios_vs2[i,1]
  p0<- q0/scenarios_vs2[i,2]
  n0<- scenarios_vs2[i,3]
  beta0<- scenarios_vs2[i,4]
  rho0<- scenarios_vs2[i,5]
  a<- foreach (j = 1:n_rep,.packages=c('ModelMetrics','MixMatrix','glmnet','stats','doSNOW','MASS'))%dopar%{
    x0<-X.3(p0,n0,rho0)
    y0<- Y.5(q0,x0,n0,beta0)
    
    foreach (h = 1:4,.packages='glmnet')%do%{
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1)
      coef.min<- coef(fit.lasso,fit.lasso$lambda.min)
      coef.1se<- coef(fit.lasso,fit.lasso$lambda.1se)
      coef.min<- coef.min[-1,]
      coef.1se<- coef.1se[-1,]
      true.correct_min<- length(which(which(coef.min!=0)< (q0+1)))
      true.correct_1se<- length(which(which(coef.1se!=0)< (q0+1)))
      fake.correct_min<- length(which(coef.min!=0))-true.correct_min
      fake.correct_1se<- length(which(coef.1se!=0))-true.correct_1se
      fake.zero_min<- q0-true.correct_min
      fake.zero_1se<- q0-true.correct_1se
      true.zero_min<- length(which(coef.min==0))-fake.zero_min
      true.zero_1se<- length(which(coef.1se==0))-fake.zero_1se
      fit$trco.min<- true.correct_min/q0
      fit$trco.1se<- true.correct_1se/q0
      fit$faco.min<- fake.correct_min/q0
      fit$faco.1se<- fake.correct_1se/q0
      fit$trze.min<- true.zero_min/q0
      fit$trze.1se<- true.zero_1se/q0
      fit$faze.min<- fake.zero_min/q0
      fit$faze.1se<- fake.zero_1se/q0
      
      return(fit) 
    }
  }
  for (j in 1:n_rep) {
    for(h in 1:4){
      trco.min[h,j]<-a[[j]][[h]]$trco.min
      trco.1se[h,j]<-a[[j]][[h]]$trco.1se
      faco.min[h,j]<-a[[j]][[h]]$faco.min
      faco.1se[h,j]<-a[[j]][[h]]$faco.1se
      trze.min[h,j]<-a[[j]][[h]]$trze.min
      trze.1se[h,j]<-a[[j]][[h]]$trze.1se
      faze.min[h,j]<-a[[j]][[h]]$faze.min
      faze.1se[h,j]<-a[[j]][[h]]$faze.1se
    }
  }
  write.matrix(trco.min,file=paste0("vs.trco.min3_",i,".csv"))
  write.matrix(trco.1se,file=paste0("vs.trco.1se3_",i,".csv"))
  write.matrix(faco.min,file=paste0("vs.faco.min3_",i,".csv"))
  write.matrix(faco.1se,file=paste0("vs.faco.1se3_",i,".csv"))
  write.matrix(trze.min,file=paste0("vs.trze.min3_",i,".csv"))
  write.matrix(trze.1se,file=paste0("vs.trze.1se3_",i,".csv"))
  write.matrix(faze.min,file=paste0("vs.faze.min3_",i,".csv"))
  write.matrix(faze.1se,file=paste0("vs.faze.1se3_",i,".csv"))
  
}


###################################################################################
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

fit<- list()
a<- list()

trco.min<- matrix(0,nrow=4,ncol=n_rep)
trco.1se<- matrix(0,nrow=4,ncol=n_rep)
faco.min<- matrix(0,nrow=4,ncol=n_rep)
faco.1se<- matrix(0,nrow=4,ncol=n_rep)
trze.min<- matrix(0,nrow=4,ncol=n_rep)
trze.1se<- matrix(0,nrow=4,ncol=n_rep)
faze.min<- matrix(0,nrow=4,ncol=n_rep)
faze.1se<- matrix(0,nrow=4,ncol=n_rep)

for (i in 1:nrow(scenarios_vs2)){
  q0<- scenarios_vs2[i,1]
  p0<- q0+scenarios_vs2[i,2]
  n0<- scenarios_vs2[i,3]
  beta0<- scenarios_vs2[i,4]
  rho0<- scenarios_vs2[i,5]
  a<- foreach (j = 1:n_rep,.packages='doSNOW')%dopar%{
    x0<-X.3(p0,n0,rho0)
    y0<- Y.6(q0,p0,n0,rho0,beta0)
    
    foreach (h = 1:4,.packages='glmnet')%do%{
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1)
      coef.min<- coef(fit.lasso,fit.lasso$lambda.min)
      coef.1se<- coef(fit.lasso,fit.lasso$lambda.1se)
      coef.min<- coef.min[-1,]
      coef.1se<- coef.1se[-1,]
      true.correct_min<- length(which(which(coef.min!=0)< (q0+1)))
      true.correct_1se<- length(which(which(coef.1se!=0)< (q0+1)))
      fake.correct_min<- length(which(coef.min!=0))-true.correct_min
      fake.correct_1se<- length(which(coef.1se!=0))-true.correct_1se
      fake.zero_min<- q0-true.correct_min
      fake.zero_1se<- q0-true.correct_1se
      true.zero_min<- length(which(coef.min==0))-fake.zero_min
      true.zero_1se<- length(which(coef.1se==0))-fake.zero_1se
      fit$trco.min<- true.correct_min/q0
      fit$trco.1se<- true.correct_1se/q0
      fit$faco.min<- fake.correct_min/q0
      fit$faco.1se<- fake.correct_1se/q0
      fit$trze.min<- true.zero_min/q0
      fit$trze.1se<- true.zero_1se/q0
      fit$faze.min<- fake.zero_min/q0
      fit$faze.1se<- fake.zero_1se/q0
      
      return(fit)
    }
  }
  for (j in 1:n_rep) {
    for(h in 1:4){
      trco.min[h,j]<-a[[j]][[h]]$trco.min
      trco.1se[h,j]<-a[[j]][[h]]$trco.1se
      faco.min[h,j]<-a[[j]][[h]]$faco.min
      faco.1se[h,j]<-a[[j]][[h]]$faco.1se
      trze.min[h,j]<-a[[j]][[h]]$trze.min
      trze.1se[h,j]<-a[[j]][[h]]$trze.1se
      faze.min[h,j]<-a[[j]][[h]]$faze.min
      faze.1se[h,j]<-a[[j]][[h]]$faze.1se
    }
  }
  write.matrix(trco.min,file=paste0("vs.trco.min4_",i,".csv"))
  write.matrix(trco.1se,file=paste0("vs.trco.1se4_",i,".csv"))
  write.matrix(faco.min,file=paste0("vs.faco.min4_",i,".csv"))
  write.matrix(faco.1se,file=paste0("vs.faco.1se4_",i,".csv"))
  write.matrix(trze.min,file=paste0("vs.trze.min4_",i,".csv"))
  write.matrix(trze.1se,file=paste0("vs.trze.1se4_",i,".csv"))
  write.matrix(faze.min,file=paste0("vs.faze.min4_",i,".csv"))
  write.matrix(faze.1se,file=paste0("vs.faze.1se4_",i,".csv"))
}

