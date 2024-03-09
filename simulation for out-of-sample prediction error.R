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



##prediction error estimation (MSE) Is it a good idea? maybe AUC? Due to Tibshirani's paper, the MSE estimation seems to make sense 
##even when it is high dimension. But what about the issue of logistic model? Or maybe log loss?  
##Search brier score!!!!!!

q.pre<- c(10,100,200) #the value of q for standard deviation test
qp.pre<- c(0.1,0.5,10/11) #the ratio between the number of non-zero coefficients to total number of coefficients

scenarios_pre<- expand.grid(q.pre=q.pre,qp.pre=qp.pre,n=n,beta.value=beta.value)

#work session
#with constant coefficients but without correlation

brier.lasso.min<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
brier.lasso.1se<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
brier.ridge.min<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
brier.ridge.1se<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
min_better.lasso<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
min_better.ridge<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
se_better.lasso<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
se_better.ridge<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
same.lasso<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)
same.ridge<- matrix(0,ncol = nrow(scenarios_pre),nrow=4)

br.lasso.min<- matrix(0,nrow=4,ncol=n_rep)
br.lasso.1se<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.min<- matrix(0,nrow=4,ncol=n_rep)
br.ridge.1se<- matrix(0,nrow=4,ncol=n_rep)
ratio.lasso<- matrix(0,nrow=4,ncol=n_rep)
ratio.ridge<- matrix(0,nrow=4,ncol=n_rep)



a<- list()
fit<- list()


for (i in 91:108){
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
  write.matrix(br.lasso.min,file=paste0("br.lasso.min",i,".csv"))
  write.matrix(br.lasso.1se,file=paste0("br.lasso.1se",i,".csv"))
  write.matrix(br.ridge.min,file=paste0("br.ridge.min",i,".csv"))
  write.matrix(br.ridge.1se,file=paste0("br.ridge.1se",i,".csv"))
  brier.lasso.min[,i]<- rowMeans(br.lasso.min)
  brier.lasso.1se[,i]<- rowMeans(br.lasso.1se)
  brier.ridge.min[,i]<- rowMeans(br.ridge.min)
  brier.ridge.1se[,i]<- rowMeans(br.ridge.1se)
  
  for (h in 1:4) {
    min_better.lasso[h,i]<- length(ratio.lasso[h,][ratio.lasso[h,]<1])/n_rep
    min_better.ridge[h,i]<- length(ratio.ridge[h,][ratio.ridge[h,]<1])/n_rep
    
    se_better.lasso[h,i]<- length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    se_better.ridge[h,i]<- length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
    
    same.lasso[h,i]<- length(ratio.lasso[h,][ratio.lasso[h,]==1])/n_rep
    same.ridge[h,i]<- length(ratio.ridge[h,][ratio.ridge[h,]==1])/n_rep
  }
}






system.time(for (i in 1:9){
  q0<- scenarios_pre[i,1]
  p0<- q0/scenarios_pre[i,2]
  n0<- scenarios_pre[i,3]
  beta0<- scenarios_pre[i,4]
  
  a<- foreach(j = 1:n_rep,.packages = c('ModelMetrics','glmnet','stats','doSNOW'))%dopar%{
    x0<-X.1(p0,n0)
    x0.test<- X.1.test(p0,n0)
    y0<- Y.1(q0,p0,n0,beta0)
    y0.test<- Y.1.test(q0,p0,n0,beta0)
    foreach(h = 1:4)%do%{
      
      fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
      fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)
      
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
  write.matrix(br.lasso.min,file=paste0("br.lasso.min",i,".csv"))
  write.matrix(br.lasso.1se,file=paste0("br.lasso.1se",i,".csv"))
  write.matrix(br.ridge.min,file=paste0("br.ridge.min",i,".csv"))
  write.matrix(br.ridge.1se,file=paste0("br.ridge.1se",i,".csv"))
  brier.lasso.min[,i]<- rowMeans(br.lasso.min)
  brier.lasso.1se[,i]<- rowMeans(br.lasso.1se)
  brier.ridge.min[,i]<- rowMeans(br.ridge.min)
  brier.ridge.1se[,i]<- rowMeans(br.ridge.1se)
  
  for (h in 1:4) {
    ra.lasso[h,i]<- 1-length(ratio.lasso[h,][ratio.lasso[h,]>1])/n_rep
    ra.ridge[h,i]<- 1-length(ratio.ridge[h,][ratio.ridge[h,]>1])/n_rep
  }
  
})
ra.lasso
ratio.lasso




for(i in 1:9){
  q0<- scenarios_sd[i+9,1]
  p0<- q0/scenarios_sd[i+9,2]
  n0<- scenarios_sd[i+9,3]
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
  
  cv.SE.lasso[,i]<- rowMeans(se.1se_lasso)
  cv.SE.ridge[,i]<- rowMeans(se.1se_ridge)
  cv.SD.lasso[,i]<- rowSds(mse.min_lasso)
  cv.SD.ridge[,i]<- rowSds(mse.min_ridge)
  SDtoSE.lasso[,i]<- cv.SD.lasso[,i]/cv.SE.lasso[,i]
  SDtoSE.ridge[,i]<- cv.SD.ridge[,i]/cv.SE.ridge[,i]
}

brier.lasso.min
brier.lasso.1se
brier.ridge.min
brier.ridge.1se


fit.lasso<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=1,parallel = TRUE)
fit.ridge<- cv.glmnet(x0,y0,nfolds = K[h],family="binomial",alpha=0,parallel = TRUE)

pre.lasso.min<- predict(fit.lasso,s = "lambda.min",newx = x0.test,type = "response")
pre.lasso.1se<- predict(fit.lasso,s = "lambda.1se",newx = x0.test,type = "response")
pre.ridge.min<- predict(fit.ridge,s = "lambda.min",newx = x0.test,type = "response")
pre.ridge.1se<- predict(fit.ridge,s = "lambda.1se",newx = x0.test,type = "response")
