library(MASS)
library(boot)
library(glmnet)
library(matrixStats)
library(MixMatrix)
library (ModelMetrics)
library(stats)
library(foreach)
library(doParallel)
registerDoParallel(cores=6)

n_rep<- 1000
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.01,0.1,0.5,0.9) #correlation between
beta.value<- c(0.2,0.4,0.6,0.8,1.0,1.5,2.0,3.0) #constant beta value
q.pre<- c(10,100,200) #the value of q for standard deviation test
qp.pre<- c(0.1,0.5,10/11) #the ratio between the number of non-zero coefficients to total number of coefficients

scenarios_pre<- expand.grid(q.pre=q.pre,qp.pre=qp.pre,n=n,beta.value=beta.value)
scenarios_pre2<- expand.grid(q.pre=q.pre,qp.pre=qp.pre,n=n,beta.value=beta.value,rho=rho)

scenarios_full<- expand.grid(q=q.pre,qp=qp.pre,n=n,beta.value=beta.value,K=K)
scenarios_full2<- expand.grid(q=q.pre,qp=qp.pre,n=n,beta.value=beta.value,rho=rho,K=K)
############################################################
#with constant coefficients but without correlation

X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set 
X.1.test<- function(Np,N){
  matrix(rnorm(N*Np*0.3),nrow=0.3*N,ncol=Np)
}

Y.1<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients
Y.1.test<- function(Q,X,N,Beta){
  z<- Beta*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(0.3*N,1,pr)
}

min_better.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
min_better.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
se_better.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
se_better.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
same.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
same.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)

br.lasso.min<- matrix(0,nrow=4,ncol=n_rep)
br.lasso.1se<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.min<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.1se<- matrix(0,nrow=4,ncol=n_rep)
ratio.lasso<- matrix(0,nrow=4,ncol=n_rep)
ratio.ridge<- matrix(0,nrow=4,ncol=n_rep)


a<- list()
fit<- list()


for (i in 1:nrow(scenarios_pre)){
  q0<- scenarios_pre[i,1]
  p0<- q0/scenarios_pre[i,2]
  n0<- scenarios_pre[i,3]
  beta0<- scenarios_pre[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = c('ModelMetrics','glmnet','MixMatrix','MASS','doSNOW'))%dopar%{
    x0<-X.1(p0,n0)
    x0.test<- X.1.test(p0,n0)
    y0<- Y.1(q0,x0,n0,beta0)
    y0.test<- Y.1.test(q0,x0.test,n0,beta0)
    foreach(h = 1:4)%do%{
      
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
          fitting <- TRUE
          
        },
        error = function(e) {
          x0<-X.1(p0,n0)
          x0.test<- X.1.test(p0,n0)
          y0<- Y.1(q0,x0,n0,beta0)
          y0.test<- Y.1.test(q0,x0.test,n0,beta0)
          
        })
      }
      
      pre.lasso.min<- predict(fit.lasso,s = "lambda.min",newx = x0.test,type = "response")
      pre.lasso.1se<- predict(fit.lasso,s = "lambda.1se",newx = x0.test,type = "response")
      pre.ridge.min<- predict(fit.ridge,s = "lambda.min",newx = x0.test,type = "response")
      pre.ridge.1se<- predict(fit.ridge,s = "lambda.1se",newx = x0.test,type = "response")
      
      fit$lmin<- brier(y0.test,pre.lasso.min)
      fit$lse<- brier(y0.test,pre.lasso.1se)
      fit$rmin<- brier(y0.test,pre.ridge.min)
      fit$rse<- brier(y0.test,pre.ridge.1se)
      
      return(fit)
    }
  }
  for(j in 1:n_rep){
    for(h in 1:4){
      br.lasso.min[h,j]<-a[[j]][[h]]$lmin
      br.lasso.1se[h,j]<-a[[j]][[h]]$lse
      br.ridge.min[h,j]<-a[[j]][[h]]$rmin
      br.ridge.1se[h,j]<-a[[j]][[h]]$rse
      
      ratio.lasso[h,j]<- br.lasso.min[h,j]/br.lasso.1se[h,j]
      ratio.ridge[h,j]<- br.ridge.min[h,j]/br.ridge.1se[h,j]
    }
  }
  
  for (h in 1:4) {
    min_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]<1])/n_rep
    min_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]<1])/n_rep
    
    se_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    se_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
    
    same.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]==1])/n_rep
    same.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]==1])/n_rep
  }
}

min_better_lasso1<- cbind(scenarios_full,c(min_better.lasso))
min_better_ridge1<- cbind(scenarios_full,c(min_better.ridge))

se_better_lasso1<- cbind(scenarios_full,c(se_better.lasso))
se_better_ridge1<- cbind(scenarios_full,c(se_better.ridge))

same_lasso1<- cbind(scenarios_full,c(same.lasso))
same_ridge1<- cbind(scenarios_full,c(same.ridge))

con_nocor<- list(min_better_lasso1=min_better_lasso1,min_better_ridge1=min_better_ridge1,
                 se_better_lasso1=se_better_lasso1,se_better_ridge1=se_better_ridge1,
                 same_lasso1=same_lasso1,same_ridge1=same_ridge1)
saveRDS(con_nocor,file = "br_con_nocor.RDS")


############################################################
#with decaying coefficients but without correlation

X.1<- function(Np,N){
  matrix(rnorm(N*Np),nrow=N,ncol=Np)
}         #Np refers to the number of total predictors,N refers to the sample size, this function is for generating training set 
X.1.test<- function(Np,N){
  matrix(rnorm(N*Np*0.3),nrow=0.3*N,ncol=Np)
}

Y.2<- function(Q,X,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with decaying coefficients without correlation
Y.2.test<- function(Q,X,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(0.3*N,1,pr)
}

min_better.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
min_better.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
se_better.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
se_better.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
same.lasso<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)
same.ridge<- matrix(0,nrow = nrow(scenarios_pre),ncol=4)

br.lasso.min<- matrix(0,nrow=4,ncol=n_rep)
br.lasso.1se<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.min<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.1se<- matrix(0,nrow=4,ncol=n_rep)
ratio.lasso<- matrix(0,nrow=4,ncol=n_rep)
ratio.ridge<- matrix(0,nrow=4,ncol=n_rep)


a<- list()
fit<- list()


for (i in 1:nrow(scenarios_pre)){
  q0<- scenarios_pre[i,1]
  p0<- q0/scenarios_pre[i,2]
  n0<- scenarios_pre[i,3]
  beta0<- scenarios_pre[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = c('ModelMetrics','glmnet','MixMatrix','MASS','doSNOW'))%dopar%{
    x0<-X.1(p0,n0)
    x0.test<- X.1.test(p0,n0)
    y0<- Y.2(q0,x0,n0,beta0)
    y0.test<- Y.2.test(q0,x0.test,n0,beta0)
    foreach(h = 1:4)%do%{
      
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
          fitting <- TRUE
          
        },
        error = function(e) {
          x0<-X.1(p0,n0)
          x0.test<- X.1.test(p0,n0)
          y0<- Y.2(q0,x0,n0,beta0)
          y0.test<- Y.2.test(q0,x0.test,n0,beta0)
          
        })
      }
      
      pre.lasso.min<- predict(fit.lasso,s = "lambda.min",newx = x0.test,type = "response")
      pre.lasso.1se<- predict(fit.lasso,s = "lambda.1se",newx = x0.test,type = "response")
      pre.ridge.min<- predict(fit.ridge,s = "lambda.min",newx = x0.test,type = "response")
      pre.ridge.1se<- predict(fit.ridge,s = "lambda.1se",newx = x0.test,type = "response")
      
      fit$lmin<- brier(y0.test,pre.lasso.min)
      fit$lse<- brier(y0.test,pre.lasso.1se)
      fit$rmin<- brier(y0.test,pre.ridge.min)
      fit$rse<- brier(y0.test,pre.ridge.1se)
      
      return(fit)
    }
  }
  for(j in 1:n_rep){
    for(h in 1:4){
      br.lasso.min[h,j]<-a[[j]][[h]]$lmin
      br.lasso.1se[h,j]<-a[[j]][[h]]$lse
      br.ridge.min[h,j]<-a[[j]][[h]]$rmin
      br.ridge.1se[h,j]<-a[[j]][[h]]$rse
      
      ratio.lasso[h,j]<- br.lasso.min[h,j]/br.lasso.1se[h,j]
      ratio.ridge[h,j]<- br.ridge.min[h,j]/br.ridge.1se[h,j]
    }
  }
  
  for (h in 1:4) {
    min_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]<1])/n_rep
    min_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]<1])/n_rep
    
    se_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    se_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
    
    same.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]==1])/n_rep
    same.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]==1])/n_rep
  }
}

min_better_lasso2<- cbind(scenarios_full,c(min_better.lasso))
min_better_ridge2<- cbind(scenarios_full,c(min_better.ridge))

se_better_lasso2<- cbind(scenarios_full,c(se_better.lasso))
se_better_ridge2<- cbind(scenarios_full,c(se_better.ridge))

same_lasso2<- cbind(scenarios_full,c(same.lasso))
same_ridge2<- cbind(scenarios_full,c(same.ridge))


dec_nocor<- list(min_better_lasso2=min_better_lasso2,min_better_ridge2=min_better_ridge2,
                 se_better_lasso2=se_better_lasso2,se_better_ridge2=se_better_ridge2,
                 same_lasso2=same_lasso2,same_ridge2=same_ridge2)
saveRDS(dec_nocor,file = "br_dec_nocor.RDS")


############################################################
#with constant coefficients with CS correlation

X.3<- function(Np,N,RHO){
 Omega <- matrix(RHO, ncol = Np, nrow = Np)
 diag(Omega)<- 1
 mvrnorm(N,rep(1,Np),Omega)
}#This is to generate train set with compound symmetry covariance
X.3.test<- function(Np,N,RHO){
Omega <- matrix(RHO, ncol = Np, nrow = Np)
diag(Omega)<- 1
mvrnorm(0.3*N,rep(1,Np),Omega)
}
Y.5<- function(Q,X,N,Beta){
z<- Beta*X[,1:Q]
pr<- 1/(1+exp(-z))
rbinom(N,1,pr)
}#This is to generate the outcome with constant coefficients and compound symmetry correlation structure
Y.5.test<- function(Q,X,N,Beta){
z<- Beta*X[,1:Q]
pr<- 1/(1+exp(-z))
rbinom(0.3*N,1,pr)
}


min_better.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
min_better.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
se_better.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
se_better.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
same.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
same.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)

br.lasso.min<- matrix(0,nrow=4,ncol=n_rep)
br.lasso.1se<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.min<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.1se<- matrix(0,nrow=4,ncol=n_rep)
ratio.lasso<- matrix(0,nrow=4,ncol=n_rep)
ratio.ridge<- matrix(0,nrow=4,ncol=n_rep)


a<- list()
fit<- list()


for (i in 1:nrow(scenarios_pre2)){
  q0<- scenarios_pre2[i,1]
  p0<- q0/scenarios_pre2[i,2]
  n0<- scenarios_pre2[i,3]
  beta0<- scenarios_pre2[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = c('ModelMetrics','glmnet','MixMatrix','MASS','doSNOW'))%dopar%{
    x0<-X.3(p0,n0)
    x0.test<- X.3.test(p0,n0)
    y0<- Y.5(q0,x0,n0,beta0)
    y0.test<- Y.5.test(q0,x0.test,n0,beta0)
    foreach(h = 1:4)%do%{
      
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
          fitting <- TRUE
          
        },
        error = function(e) {
          x0<-X.3(p0,n0)
          x0.test<- X.3.test(p0,n0)
          y0<- Y.5(q0,x0,n0,beta0)
          y0.test<- Y.5.test(q0,x0.test,n0,beta0)
          
        })
      }
      
      pre.lasso.min<- predict(fit.lasso,s = "lambda.min",newx = x0.test,type = "response")
      pre.lasso.1se<- predict(fit.lasso,s = "lambda.1se",newx = x0.test,type = "response")
      pre.ridge.min<- predict(fit.ridge,s = "lambda.min",newx = x0.test,type = "response")
      pre.ridge.1se<- predict(fit.ridge,s = "lambda.1se",newx = x0.test,type = "response")
      
      fit$lmin<- brier(y0.test,pre.lasso.min)
      fit$lse<- brier(y0.test,pre.lasso.1se)
      fit$rmin<- brier(y0.test,pre.ridge.min)
      fit$rse<- brier(y0.test,pre.ridge.1se)
      
      return(fit)
    }
  }
  for(j in 1:n_rep){
    for(h in 1:4){
      br.lasso.min[h,j]<-a[[j]][[h]]$lmin
      br.lasso.1se[h,j]<-a[[j]][[h]]$lse
      br.ridge.min[h,j]<-a[[j]][[h]]$rmin
      br.ridge.1se[h,j]<-a[[j]][[h]]$rse
      
      ratio.lasso[h,j]<- br.lasso.min[h,j]/br.lasso.1se[h,j]
      ratio.ridge[h,j]<- br.ridge.min[h,j]/br.ridge.1se[h,j]
    }
  }
  
  for (h in 1:4) {
    min_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]<1])/n_rep
    min_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]<1])/n_rep
    
    se_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    se_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
    
    same.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]==1])/n_rep
    same.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]==1])/n_rep
  }
}

min_better_lasso3<- cbind(scenarios_full2,c(min_better.lasso))
min_better_ridge3<- cbind(scenarios_full2,c(min_better.ridge))

se_better_lasso3<- cbind(scenarios_full2,c(se_better.lasso))
se_better_ridge3<- cbind(scenarios_full2,c(se_better.ridge))

same_lasso3<- cbind(scenarios_full2,c(same.lasso))
same_ridge3<- cbind(scenarios_full2,c(same.ridge))


con_cscor<- list(min_better_lasso3=min_better_lasso3,min_better_ridge3=min_better_ridge3,
                 se_better_lasso3=se_better_lasso3,se_better_ridge3=se_better_ridge3,
                 same_lasso3=same_lasso3,same_ridge3=same_ridge3)
saveRDS(con_cscor,file = "br_con_cscor.RDS")


###############################################################
#with decaying coefficients with CS correlation

X.3<- function(Np,N,RHO){
  Omega <- matrix(RHO, ncol = Np, nrow = Np)
  diag(Omega)<- 1
  mvrnorm(N,rep(1,Np),Omega)
}#This is to generate train set with compound symmetry covariance
X.3.test<- function(Np,N,RHO){
  Omega <- matrix(RHO, ncol = Np, nrow = Np)
  diag(Omega)<- 1
  mvrnorm(0.3*N,rep(1,Np),Omega)
  
}

Y.6<- function(Q,X,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(N,1,pr)
}#This is to generate the outcome with decaying coefficients and compound symmetry correlation
Y.6.test<- function(Q,X,N,Beta){
  divided<- c(1:Q)
  Coef<- 1/divided
  z<- Beta*Coef*X[,1:Q]
  pr<- 1/(1+exp(-z))
  rbinom(0.3*N,1,pr)
}

min_better.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
min_better.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
se_better.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
se_better.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
same.lasso<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)
same.ridge<- matrix(0,nrow = nrow(scenarios_pre2),ncol=4)

br.lasso.min<- matrix(0,nrow=4,ncol=n_rep)
br.lasso.1se<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.min<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.1se<- matrix(0,nrow=4,ncol=n_rep)
ratio.lasso<- matrix(0,nrow=4,ncol=n_rep)
ratio.ridge<- matrix(0,nrow=4,ncol=n_rep)


a<- list()
fit<- list()


for (i in 1:nrow(scenarios_pre2)){
  q0<- scenarios_pre2[i,1]
  p0<- q0/scenarios_pre2[i,2]
  n0<- scenarios_pre2[i,3]
  beta0<- scenarios_pre2[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = c('ModelMetrics','glmnet','MixMatrix','MASS','doSNOW'))%dopar%{
    x0<-X.3(p0,n0)
    x0.test<- X.3.test(p0,n0)
    y0<- Y.6(q0,x0,n0,beta0)
    y0.test<- Y.6.test(q0,x0.test,n0,beta0)
    foreach(h = 1:4)%do%{
      
      fitting <- FALSE
      
      while(!fitting) {
        
        tmp <- tryCatch({
          
          fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
          fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
          fitting <- TRUE
          
        },
        error = function(e) {
          x0<-X.3(p0,n0)
          x0.test<- X.3.test(p0,n0)
          y0<- Y.6(q0,x0,n0,beta0)
          y0.test<- Y.6.test(q0,x0.test,n0,beta0)
          
        })
      }
      
      pre.lasso.min<- predict(fit.lasso,s = "lambda.min",newx = x0.test,type = "response")
      pre.lasso.1se<- predict(fit.lasso,s = "lambda.1se",newx = x0.test,type = "response")
      pre.ridge.min<- predict(fit.ridge,s = "lambda.min",newx = x0.test,type = "response")
      pre.ridge.1se<- predict(fit.ridge,s = "lambda.1se",newx = x0.test,type = "response")
      
      fit$lmin<- brier(y0.test,pre.lasso.min)
      fit$lse<- brier(y0.test,pre.lasso.1se)
      fit$rmin<- brier(y0.test,pre.ridge.min)
      fit$rse<- brier(y0.test,pre.ridge.1se)
      
      return(fit)
    }
  }
  for(j in 1:n_rep){
    for(h in 1:4){
      br.lasso.min[h,j]<-a[[j]][[h]]$lmin
      br.lasso.1se[h,j]<-a[[j]][[h]]$lse
      br.ridge.min[h,j]<-a[[j]][[h]]$rmin
      br.ridge.1se[h,j]<-a[[j]][[h]]$rse
      
      ratio.lasso[h,j]<- br.lasso.min[h,j]/br.lasso.1se[h,j]
      ratio.ridge[h,j]<- br.ridge.min[h,j]/br.ridge.1se[h,j]
    }
  }
  
  for (h in 1:4) {
    min_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]<1])/n_rep
    min_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]<1])/n_rep
    
    se_better.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    se_better.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
    
    same.lasso[i,h]<- length(ratio.lasso[h,][ratio.lasso[h,]==1])/n_rep
    same.ridge[i,h]<- length(ratio.ridge[h,][ratio.ridge[h,]==1])/n_rep
  }
}

min_better_lasso4<- cbind(scenarios_full2,c(min_better.lasso))
min_better_ridge4<- cbind(scenarios_full2,c(min_better.ridge))

se_better_lasso4<- cbind(scenarios_full2,c(se_better.lasso))
se_better_ridge4<- cbind(scenarios_full2,c(se_better.ridge))

same_lasso4<- cbind(scenarios_full2,c(same.lasso))
same_ridge4<- cbind(scenarios_full2,c(same.ridge))

dec_cscor<- list(min_better_lasso4=min_better_lasso4,min_better_ridge4=min_better_ridge4,
                 se_better_lasso4=se_better_lasso4,se_better_ridge4=se_better_ridge4,
                 same_lasso4=same_lasso4,same_ridge4=same_ridge4)
saveRDS(dec_cscor,file = "br_dec_cscor.RDS")
