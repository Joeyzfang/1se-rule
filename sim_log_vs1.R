
library(boot)
library(glmnet)
library(matrixStats)
library(MixMatrix)
library(foreach)
library (ModelMetrics)
library(doParallel)
registerDoParallel(cores=6)
##1. Generate data set
#parameters settings
n_rep<- 1000
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.1,0.7,0.9) #correlation between
beta.value<- c(0.1,0.5,1.0,2.0,3.0) #constant beta value

X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set
x1<- X.1(10,100)
Y.1<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients

##Variable selection (only for lasso)
q.vs<- c(10,100,1000)    #This is the vector of potential value of q for variable selection
qp.vs<- c(0.1,0.5,10/11) #In this section, qp.vs refers to the difference between q and p which is equal to p - q.
scenarios_vs<- expand.grid(q.vs=q.vs,qp.vs=qp.vs,n=n,beta.value=beta.value)
fit<- list()
a<- list()

#with constant coefficients but without correlation
DP.1<- matrix(0,nrow=4,ncol=nrow(scenarios_vs))
trco.min<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
trco.1se<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
faco.min<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
faco.1se<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
trze.min<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
trze.1se<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
faze.min<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
faze.1se<- array(0,dim = c(nrow(scenarios_vs),4,n_rep))
diff.P.1<- array(0,dim = c(nrow(scenarios_vs),4,n_rep)) #change my mind from ratio to difference, ratio could raise some issue about mean

for(i in 10:18){
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
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
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
      trco.min[i,h,j]<-a[[j]][[h]]$trco.min
      trco.1se[i,h,j]<-a[[j]][[h]]$trco.1se
      faco.min[i,h,j]<-a[[j]][[h]]$faco.min
      faco.1se[i,h,j]<-a[[j]][[h]]$faco.1se
      trze.min[i,h,j]<-a[[j]][[h]]$trze.min
      trze.1se[i,h,j]<-a[[j]][[h]]$trze.1se
      faze.min[i,h,j]<-a[[j]][[h]]$faze.min
      faze.1se[i,h,j]<-a[[j]][[h]]$faze.1se
    }
  }
}
trco.min
write.matrix(trpo.min,file="trpo.min1_1.csv")
write.matrix(trpo.1se,file="trpo.1se1_1.csv")

write.matrix(fapo.min,file="fapo.min1_1.csv")
write.matrix(fapo.1se,file="fapo.1se1_1.csv")

write.matrix(trne.min,file="trne.min1_1.csv")
write.matrix(trne.1se,file="trne.1se1_1.csv")

write.matrix(fane.min,file="fane.min1_1.csv")
write.matrix(fane.1se,file="fane.1se1_1.csv")

