setwd("C:/Users/joeyt/Desktop/research_project/HPC/data_to_HPC")
load(file = "sd_lasso.RData")
install.packages("Hmisc")
library(Hmisc)
library(dplyr)
install.packages("vctrs")
library(vctrs)
install.packages("BMA")
library(BMA)
library(rms)
library(parallel)
library(car)
library(FactoMineR)
install.packages("dummy")
library(dummy)
library(foreach)
library(MASS)
registerDoParallel(cores=6)
#define a exclusion function
combn_with_exclusion<- function(X, exclude){
  full <- X
  # remove any columns that have all elements of `exclude`
  full[, !apply(full, 2, function(y) all(exclude %in% y))]
}

inter<- function(VAR,n_VAR,exclude){
  comb_inter<-c(NA)
  for(i in 2:n_VAR){
    comb<-combn(VAR,i)
    for(h in 1:ncol(exclude)){
      comb<- combn_with_exclusion(comb, exclude[,h])
    }
    for(j in 1:ncol(comb)){
      comb_inter<- rbind(comb_inter,paste0(comb[1:nrow(comb),j],collapse = ":"))
      
    }
    comb_inter<- na.omit(comb_inter)
  }
  return(comb_inter)
}
#define a forward stepwise function using a score system combining R-squared, AIC and p-value
stepwise_forward<- function(FIT,VAR,n_VAR,exclude,y,P){
  for(i in 2:n_VAR){
    VAR<- sample(VAR,length(VAR),replace = FALSE)
    comb<-combn(VAR,i)
    for(h in 1:ncol(exclude)){
      comb<- combn_with_exclusion(comb, exclude[,h])
    }
    for(j in 1:ncol(comb)){
      comb_inter<- paste0(comb[1:nrow(comb),j],collapse = ":")
      fit<- update(FIT,paste("~ . +",comb_inter),train)
      n<- 0
      if(is.na(summary(fit)$coefficients[nrow(summary(fit)$coefficients),4])){
        FIT<-FIT
      }
      else if(summary(fit)$adj.r.squared>summary(FIT)$adj.r.squared){
        n<- n+1
        if(mean((exp(predict(fit))-y)^2)<mean((exp(predict(FIT))-y)^2)){
          n<- n+1
          if(summary(FIT)$coefficients[nrow(summary(FIT)$coefficients),4]<P){
            n<- n+0.5
          }
        }
      }
      if(n>1){
        FIT<-fit
      } else {
        FIT<- FIT
      }
    }
  }
  return(FIT)
}
b<- list()
a<-foreach(i = 1:5)%do%{
  b$m<-i
  b$n<-i^2
  return(b)
}
which.max(sapply(a,function(x) x$n))
folds <- cut(seq(1,nrow(DATA)),breaks=k,labels=FALSE)
cvstats<- list()
#cross validation for forward step
models<- foreach (i = 1:k,.packages = c('glmnet','rms'))%dopar% {
  # Set up the training set and validation set for this fold
  train <- DATA[folds != i, ]
  valid <- DATA[folds == i, ]
  comb_inter<- paste0(comb[1:nrow(comb),j],collapse = ":")
  cvstats$fit<- update(FIT,paste("~ . +",comb_inter),train)
  # Evaluate the model with interaction terms on the validation set
  pre <- predict(cvstats$fit, newdata=valid)
  cvstats$for_mse <- mean((exp(pre) - valid$true_se)^2)
  return(cvstats)
}
fit<- models[[which.min(sapply(cvstats, function(x) x$mse))]]$fit
#define a bidirectional stepwise function
set.seed(123)
bistep<- function(VAR,FIT,DATA){
  while (TRUE) {
    #divide the data into k folds

    #check which variables are remained to be added
    remaining_vars<- setdiff(VAR,attr(terms(FIT),"term.labels"))
    if(length(remaining_vars)==0){
      break
    }
    
    #forward step

      
      # Fit the model with one more term
      if(length(attr(terms(FIT),"term.labels"))==0){
        forward_candidates <- lapply(remaining_vars, 
                                     function(x) formula(paste("log(true_se)~", x)))
      } else {
        forward_candidates <- lapply(remaining_vars, 
                                     function(x) formula(paste( "log(true_se)~", 
                                                                paste(attr(terms(FIT),"term.labels"), 
                                                                      collapse = " + "), "+", x)))
      }
      forward_models<- lapply(forward_candidates, function(x) lm(x,data=DATA))
      forward_adj.r<- sapply(forward_models, function(x) summary(x)$adj.r.squared)
      forward_best<- forward_models[[which.max(forward_adj.r)]]
      

      
      # Fit the model with one less term
      if(length(attr(terms(FIT),"term.labels"))==0){
        if(summary(forward_best)$adj.r.squared > summary(FIT)$adj.r.squared){
          FIT<- forward_best
        }
      } else {
        backward_candidates<- lapply(attr(terms(FIT),"term.labels"), 
                                     function(x) {
                                       if(length(attr(terms(FIT),"term.labels"))==1){
                                         formula(paste("log(true_se)~ 1"))
                                       }else{
                                         formula(paste("log(true_se)~", 
                                                       paste(attr(terms(FIT),"term.labels")[-which(attr(terms(FIT),"term.labels") == x)], 
                                                             collapse = " + ")))
                                       }
                                     })
        backward_models<- lapply(backward_candidates, function(x) lm(x,data=DATA))
        backward_adj.r<- sapply(backward_models, function(x) summary(x)$adj.r.squared)
        backward_best<- backward_models[[which.max(backward_adj.r)]]
       
      # Compare the mse of the forward and backward steps
      if (mean((exp(predict(forward_best))-DATA$true_se)^2) < mean((exp(predict(FIT))-DATA$true_se)^2) || 
          mean((exp(predict(backward_best))-DATA$true_se)^2) < mean((exp(predict(FIT))-DATA$true_se)^2)) {
        if (mean((exp(predict(forward_best))-DATA$true_se)^2) < mean((exp(predict(backward_best))-DATA$true_se)^2)) {
          FIT <- forward_best
        } else {
          FIT <- backward_best
        }
      } else {
        break
      }
    }
  }
  return(FIT)
}

cv.mse<- function(FIT,K,DATA){
  k <- K
  
  # Split the dataset into k folds
  set.seed(123)  # Set seed for reproducibility
  fold_indices <- cut(seq(1, nrow(DATA)), breaks = k, labels = FALSE)
  
  # Create vectors to store performance metrics
  mse <- numeric(k)
  
  # Perform k-fold cross-validation
  for (i in 1:k) {
    # Create the training set and validation set
    train_set <- DATA[fold_indices != i, ]
    val_set <- DATA[fold_indices == i, ]
    
    # Fit the model on the training set
    fit <- lm(formula(FIT), data = train_set)
    
    # Make predictions on the validation set
    val_pred <- predict(fit, newdata = val_set)
    
    # Calculate the mean squared error (MSE)
    mse[i] <- mean((val_pred - log(val_set$true_se))^2)
  }
  
  # Calculate the average MSE
  avg_mse <- mean(mse)
  return(avg_mse)
}

cv.bistep<- function(VAR,FIT,DATA,K){
  mse_lowest<- cv.mse(FIT,K,DATA)
  while (TRUE) {
    #divide the data into k folds
    
    #check which variables are remained to be added
    remaining_vars<- setdiff(VAR,attr(terms(FIT),"term.labels"))
    if(length(remaining_vars)==0){
      break
    }
    
    #forward step
    
    
    # Fit the model with one more term
    if(length(attr(terms(FIT),"term.labels"))==0){
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste("log(true_se)~", x)))
    } else {
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste( "log(true_se)~", 
                                                              paste(attr(terms(FIT),"term.labels"), 
                                                                    collapse = " + "), "+", x)))
    }
    forward_models<- lapply(forward_candidates, function(x) lm(x,data=DATA))
    forward_adj.r<- sapply(forward_models, function(x) summary(x)$adj.r.squared)
    forward_best<- forward_models[[which.max(forward_adj.r)]]
    
    
    
    # Fit the model with one less term
    if(length(attr(terms(FIT),"term.labels"))==0){
      if(summary(forward_best)$adj.r.squared > summary(FIT)$adj.r.squared){
        FIT<- forward_best
      }
    } else {
      backward_candidates<- lapply(attr(terms(FIT),"term.labels"), 
                                   function(x) {
                                     if(length(attr(terms(FIT),"term.labels"))==1){
                                       formula(paste("log(true_se)~ 1"))
                                     }else{
                                       formula(paste("log(true_se)~", 
                                                     paste(attr(terms(FIT),"term.labels")[-which(attr(terms(FIT),"term.labels") == x)], 
                                                           collapse = " + ")))
                                     }
                                   })
      backward_models<- lapply(backward_candidates, function(x) lm(x,data=DATA))
      backward_adj.r<- sapply(backward_models, function(x) summary(x)$adj.r.squared)
      backward_best<- backward_models[[which.max(backward_adj.r)]]
      
      # Compare the mse of the forward and backward steps in a cv process
      mse_forward<- cv.mse(forward_best,K,DATA)
      mse_backward<- cv.mse(backward_best,K,DATA)
      if (mse_forward < mse_lowest || 
          mse_backward < mse_lowest) {
        if (mse_forward < mse_backward) {
          FIT <- forward_best
          mse_lowest<- mse_forward
        } else {
          FIT <- backward_best
          mse_lowest<- mse_backward
        }
      } else {
        break
      }
    }
  }
  return(FIT)
}


bistep2<- function(VAR,FIT,DATA){
  while (TRUE) {
    #divide the data into k folds
    
    #check which variables are remained to be added
    remaining_vars<- setdiff(VAR,attr(terms(FIT),"term.labels"))
    if(length(remaining_vars)==0){
      break
    }
    
    #forward step
    
    
    # Fit the model with one more term
    if(length(attr(terms(FIT),"term.labels"))==0){
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste("log(true_se)~", x)))
    } else {
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste( "log(true_se)~", 
                                                              paste(attr(terms(FIT),"term.labels"), 
                                                                    collapse = " + "), "+", x)))
    }
    forward_models<- lapply(forward_candidates, function(x) lm(x,data=DATA))
    forward_adj.r<- sapply(forward_models, function(x) summary(x)$adj.r.squared)
    forward_best<- forward_models[[which.max(forward_adj.r)]]
    
    
    
    # Fit the model with one less term
    if(length(attr(terms(FIT),"term.labels"))==0){
      if(summary(forward_best)$adj.r.squared > summary(FIT)$adj.r.squared){
        FIT<- forward_best
      }
    } else {
      backward_candidates<- lapply(attr(terms(FIT),"term.labels"), 
                                   function(x) {
                                     if(length(attr(terms(FIT),"term.labels"))==1){
                                       formula(paste("log(true_se)~ 1"))
                                     }else{
                                       formula(paste("log(true_se)~", 
                                                     paste(attr(terms(FIT),"term.labels")[-which(attr(terms(FIT),"term.labels") == x)], 
                                                           collapse = " + ")))
                                     }
                                   })
      backward_models<- lapply(backward_candidates, function(x) lm(x,data=DATA))
      backward_adj.r<- sapply(backward_models, function(x) summary(x)$adj.r.squared)
      backward_best<- backward_models[[which.max(backward_adj.r)]]
      
      # Compare the mse of the forward and backward steps
      if (summary(forward_best)$adj.r.squared > summary(FIT)$adj.r.squared || 
          summary(backward_best)$adj.r.squared > summary(FIT)$adj.r.squared) {
        if (summary(forward_best)$adj.r.squared > summary(backward_best)$adj.r.squared) {
          FIT <- forward_best
        } else {
          FIT <- backward_best
        }
      } else {
        break
      }
    }
  }
  return(FIT)
}
#Bayesian Model Averaging (BMA) based interaction selection function
bistep<- function(VAR,FIT,DATA){
  while (TRUE) {
    #divide the data into k folds
    
    #check which variables are remained to be added
    remaining_vars<- setdiff(VAR,attr(terms(FIT),"term.labels"))
    if(length(remaining_vars)==0){
      break
    }
    
    #forward step
    
    
    # Fit the model with one more term
    if(length(attr(terms(FIT),"term.labels"))==0){
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste("log(true_se)~", x)))
    } else {
      forward_candidates <- lapply(remaining_vars, 
                                   function(x) formula(paste( "log(true_se)~", 
                                                              paste(attr(terms(FIT),"term.labels"), 
                                                                    collapse = " + "), "+", x)))
    }
    forward_models<- lapply(forward_candidates, function(x) lm(x,DATA))
    forward_result<- sapply(forward_models, function(x) bma(x)$postprob)
    forward_best<- forward_models[[which.max(forward_result)]]
    
    
    
    # Fit the model with one less term
    if(length(attr(terms(FIT),"term.labels"))==0){
      if(summary(forward_best)$adj.r.squared > summary(FIT)$adj.r.squared){
        FIT<- forward_best
      }
    } else {
      backward_candidates<- lapply(attr(terms(FIT),"term.labels"), 
                                   function(x) {
                                     if(length(attr(terms(FIT),"term.labels"))==1){
                                       formula(paste("log(true_se)~ 1"))
                                     }else{
                                       formula(paste("log(true_se)~", 
                                                     paste(attr(terms(FIT),"term.labels")[-which(attr(terms(FIT),"term.labels") == x)], 
                                                           collapse = " + ")))
                                     }
                                   })
      backward_models<- lapply(backward_candidates, function(x) bic.glm(x,data=DATA))
      backward_result<- sapply(backward_models, function(x) bma(x)$postprob)
      backward_best<- backward_models[[which.max(backward_result)]]
      
      # Compare the mse of the forward and backward steps
      if (summary(forward_best)$adj.r > summary(FIT)$adj.r || 
          summary(backward_best)$adj.r > summary(FIT)$adj.r) {
        if (summary(forward_best)$adj.r > summary(backward_best)$adj.r) {
          FIT <- forward_best
        } else {
          FIT <- backward_best
        }
      } else {
        break
      }
    }
  }
  return(FIT)
}




#data management
sd_for_lowd<- rbind(data_00,data_01,data_10,data_11)
sd_for_lowd<- cbind(sd_for_lowd,sd)
save(sd_for_lowd,file = "sd_for_lowd.RData")
sd_for_lowd<- cbind(sd_for_lowd[,5:6],sd_for_lowd[,1:4],true_se=sd_for_lowd[,7])
sd.lasso<- rbind(sd.lasso.new,sd_for_lowd)
save(sd.lasso,file = "newdata_for_sdlasso.RData")
data1<- sd.lasso[sd.lasso$rho %in% 0,1:7]
data2<- sd.lasso[!sd.lasso$rho %in% 0,1:7]
data3<- sd.lasso[sd.lasso$coef %in% 0,1:7]
data4<- sd.lasso[sd.lasso$coef %in% 1,1:7]
data3.1<- data3[data3$rho %in% 0,1:7]
data3.2<- data3[!data3$rho %in% 0,1:7]
data4.1<- data4[data4$rho %in% 0,1:7]
data4.2<- data4[!data4$rho %in% 0,1:7]

##try to use n/q instead of q itself
data3.1_new<- mutate(data3.1,n.q=n/q) 
data3.2_new<- mutate(data3.2,n.q=n/q)
true_q_4.1<- unlist(lapply(data4.1$q, function(x) sum(1/c(1:x))))
true_q_4.2<- unlist(lapply(data4.2$q, function(x) sum(1/c(1:x))))
data4.1_new<- mutate(data4.1,n.q=n/q)
data4.2_new<- mutate(data4.2,n.q=n/q)

#using 1/ln(x+1) transformation
data4.1_new<- cbind(data4.1_new[,1:2],ln=1/(1+exp(-data4.1_new$n.q)),data4.1[,4:7])
data4.2_new<- cbind(data4.2_new[,1:2],ln=1/(1+exp(-data4.2_new$n.q)),data4.2[,4:7])
data4.1_new<- cbind(data4.1_new[,1:4],logn=log(data4.1_new$n),data4.1_new[,6:7])
data4.2_new<- cbind(data4.2_new[,1:4],logn=log(data4.2_new$n),data4.2_new[,6:7])

transformed_y<- boxcox(data3.1_new$true_se)

str_new1<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3), data = data3.1_new)
str_new2<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2), data = data3.2_new)
str_new3<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3), data = data4.1_new)
str_new4<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2), data = data4.2_new)

str_high1<- glm(log(true_se) ~ I(exp(n.q))+I(exp(n.q)^2)+I(exp(n.q)^3)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3),data = data3.1_new[data3.1_new$n.q<1,])
str_high2<- glm(log(true_se) ~ I(exp(n.q))+I(exp(n.q)^2)+I(exp(n.q)^3)+q.p+I(q.p^2)
               +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2),data = data3.2_new[data3.2_new$n.q<1,])
str_high3<- glm(log(true_se) ~ I(exp(n.q))+I(exp(n.q)^2)+I(exp(n.q)^3)+q.p+I(q.p^2)
               +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3),data = data4.1_new[data4.1_new$n.q<1,])
str_high4<- glm(log(true_se) ~ I(exp(n.q))+I(exp(n.q)^2)+I(exp(n.q)^3)+q.p+I(q.p^2)
               +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2),data = data4.2_new[data4.2_new$n.q<1,])

#interaction and variable matrix
var1 <- c("K","I(K^2)","I(K^3)","I(log(n.q))","I(log(n.q)^2)","I(log(n.q)^3)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)")
var2 <- c("K","I(K^2)","I(K^3)","I(log(n.q))","I(log(n.q)^2)","I(log(n.q)^3)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","rho","I(rho^2)")
var3 <- c("K","I(K^2)","I(K^3)","rcs(log(n.q),4)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","q","I(q^2)","I(q^3)")
var4 <- c("K","I(K^2)","I(K^3)","rcs(log(n.q),4)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","rho","I(rho^2)","q","I(q^2)","I(q^3)")
ex1<- cbind(combn(var1[1:3],2),combn(var1[4:6],2),combn(var1[7:8],2),combn(var1[9:11],2))
ex2<- cbind(combn(var2[1:3],2),combn(var2[4:6],2),combn(var2[7:8],2),combn(var2[9:11],2),combn(var2[12:13],2))
ex1
ex2
inter_var1<- inter(var1,3,ex1)
inter_var1
all_var1<- rbind(matrix(var1),inter_var1)
all_var1
inter_var2<- inter(var2,3,ex2)
inter_var2
all_var2<- rbind(matrix(var2),inter_var2)
all_var2
exclusion3<- cbind(combn(var3[1:3],2),combn(var3[5:6],2),combn(var3[7:9],2),combn(var3[10:12],2))
exclusion3
exclusion4<- cbind(combn(var4[1:3],2),combn(var4[5:6],2),combn(var4[7:9],2),combn(var4[10:11],2),combn(var4[12:14],2))
exclusion4
inter_var3<- inter(var3,2,exclusion3)
inter_var3
all_var3<- rbind(matrix(var3),inter_var3)
all_var3
inter_var4<- inter(var4,2,exclusion4)
all_var4<- rbind(matrix(var4),inter_var4)
all_var4
inter_var5<- inter(var3,3,exclusion3)
all_var5<- rbind(matrix(var3),inter_var5)
all_var5
inter_var6<- inter(var4,3,exclusion4)
inter_var6
all_var6<- rbind(matrix(var4),inter_var6)
all_var6
##stepwise bidirection fit
##high dimension
step_3.1high<- bistep(all_var1,str_high1,data3.1_new[data3.1_new$n.q<1,])
step_3.2high<- bistep(all_var2,str_high2,data3.2_new[data3.2_new$n.q<1,])
step_4.1high<- bistep(all_var1,str_high3,data4.1_new[data4.1_new$n.q<1,])
step_4.2high<- bistep(all_var2,str_high4,data4.2_new[data4.2_new$n.q<1,])
step_new3.1<- cv.bistep(all_var3,str_new1,data3.1_new,5)
step_new3.2<- cv.bistep(all_var4,str_new2,data3.2_new,5)
step_new4.1<- cv.bistep(all_var3,str_new3,data4.1_new,5)
step_new4.2<- cv.bistep(all_var4,str_new4,data4.2_new,5)
forward1<- stepwise_forward(str_new1,var3,3,exclusion3,data3.1_new$true_se,0.05)
forward2<- stepwise_forward(str_new2,var4,3,exclusion4,data3.2_new$true_se,0.05)
forward3<- stepwise_forward(str_new3,var3,3,exclusion3,data4.1_new$true_se,0.05)
forward4<- stepwise_forward(str_new4,var4,3,exclusion4,data4.2_new$true_se,0.05)

summary(step_3.1high)
summary(step_3.2high)
summary(step_4.1high)
summary(step_4.2high)
summary(step_new3.1)
summary(step_new3.2)
summary(step_new4.1)
summary(step_new4.2)
summary(forward1)
summary(forward2)
summary(forward3)
summary(forward4)

coefficients1 <- coef(step_new3.1)
data1<- data.frame(coefficients1,p_value=summary(step_new3.1)$coefficients[, "Pr(>|t|)"])
data1<- mutate(parameter=rownames(data1),data1)
data1<- cbind(parameter=data1$parameter,coefficient=data1$coefficients1,p_value=data1$p_value)
coefficients2 <- coef(step_new3.2)
data2<- data.frame(coefficients2,p_value=summary(step_new3.2)$coefficients[, "Pr(>|t|)"])
data2<- mutate(parameter=rownames(data2),data2)
data2<- cbind(parameter=data2$parameter,coefficient=data2$coefficients2,p_value=data2$p_value)
coefficients3 <- coef(step_new4.1)
data3<- data.frame(coefficients3,p_value=summary(step_new4.1)$coefficients[, "Pr(>|t|)"])
data3<- mutate(parameter=rownames(data3),data3)
data3<- cbind(parameter=data3$parameter,coefficient=data3$coefficients3,p_value=data3$p_value)
coefficients4 <- coef(step_new4.2)
data4<- data.frame(coefficients4,p_value=summary(step_new4.2)$coefficients[, "Pr(>|t|)"])
data4<- mutate(parameter=rownames(data4),data4)
data4<- cbind(parameter=data4$parameter,coefficient=data4$coefficients4,p_value=data4$p_value)

write.csv(data1,file = "con_nocor_model.csv")
write.csv(data2,file = "con_cor_model.csv")
write.csv(data3,file = "dec_nocor_model.csv")
write.csv(data4,file = "dec_cor_model.csv")

cex_size <- 1.5  # Adjust the scaling factor as per your preference
par(cex.lab = cex_size,cex.main=cex_size)
plot(predict(step_3.1high),residuals(step_3.1high))
plot(predict(step_3.2high),residuals(step_3.2high))
plot(predict(step_4.1high),residuals(step_4.1high))
plot(predict(step_4.2high),residuals(step_4.2high))
plot(predict(step_new3.1),residuals(step_new3.1),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="constant coefficient with independent predictors")
plot(predict(step_new3.2),residuals(step_new3.2),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="constant coefficient with correlated predictors")
plot(predict(step_new4.1),residuals(step_new4.1),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="decaying coefficient with independent predictors")
plot(predict(step_new4.2),residuals(step_new4.2),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="decaying coefficient with correlated predictors")
plot(predict(forward1),residuals(forward1))
plot(predict(forward2),residuals(forward2))
plot(predict(forward3),residuals(forward3))
plot(predict(forward4),residuals(forward4))

pre_high1<- predict(step_3.1high)
pre_high2<- predict(step_3.2high)
pre_high3<- predict(step_4.1high)
pre_high4<- predict(step_4.2high)
prenew1<- predict(step_new3.1)
prenew2<- predict(step_new3.2)
prenew3<- predict(step_new4.1)
prenew4<- predict(step_new4.2)
preforward1<- predict(forward1)
preforward2<- predict(forward2)
preforward3<- predict(forward3)
preforward4<- predict(forward4)

outcome3.1high<- data3.1_new[data3.1_new$n.q<1,7]
outcome3.2high<- data3.2_new[data3.2_new$n.q<1,7]
outcome4.1high<- data4.1_new[data4.1_new$n.q<1,7]
outcome4.2high<- data4.2_new[data4.2_new$n.q<1,7]
outcome3.1<- log(data3.1[,7])
outcome3.2<- log(data3.2[,7])
outcome4.1<- log(data4.1[,7])
outcome4.2<- log(data4.2[,7])

biashigh1<- (exp(pre_high1)-outcome3.1high)/outcome3.1high
biashigh2<- (exp(pre_high2)-outcome3.2high)/outcome3.2high
biashigh3<- (exp(pre_high3)-outcome4.1high)/outcome4.1high
biashigh4<- (exp(pre_high4)-outcome4.2high)/outcome4.2high
bias_new1<- (prenew1-outcome3.1)/outcome3.1
bias_new2<- (prenew2-outcome3.2)/outcome3.2
bias_new3<- (prenew3-outcome4.1)/outcome4.1
bias_new4<- (prenew4-outcome4.2)/outcome4.2
biasforward1<- (preforward1-outcome3.1)/outcome3.1
biasforward2<- (preforward2-outcome3.2)/outcome3.2
biasforward3<- (preforward3-outcome4.1)/outcome4.1
biasforward4<- (preforward4-outcome4.2)/outcome4.2


mse3<- mean((exp(pre3)-outcome4.1)^2)
msenew1<- mean((exp(prenew1)-outcome3.1)^2)
msenew2<- mean((exp(prenew2)-outcome3.2)^2)
msenew3<- mean((exp(prenew3)-outcome4.1)^2)
msenew4<- mean((exp(prenew4)-outcome4.2)^2)

plot(outcome3.1)
plot(outcome3.2)
plot(outcome4.1)
plot(outcome4.2)

mse3
msenew1
msenew2
msenew3
msenew4

par(mfrow=c(1,1))
par(mfrow=c(2,2))
plot(biashigh1,data3.1_new[data3.1_new$n.q<1,]$n.q)
plot(biashigh2,data3.2_new[data3.2_new$n.q<1,]$n.q)
plot(biashigh3,data4.1_new[data4.1_new$n.q<1,]$n.q)
plot(biashigh4,data4.2_new[data4.2_new$n.q<1,]$n.q)
plot(data3.1_new$n.q,bias_new1,ylab="in-sample bias1",xlab="n/q",ylim = c(-0.08,0.08),main="constant coefficient with independent predictors")
plot(data3.2_new$n.q,bias_new2,ylab="in-sample bias2",xlab="n/q",ylim = c(-0.08,0.08),main="constant coefficient with correlated predictors")
plot(data4.1_new$n.q,bias_new3,ylab="in-sample bias3",xlab="n/q",ylim = c(-0.08,0.08),main="decaying coefficient with independent predictors")
plot(data4.2_new$n.q,bias_new4,ylab="in-sample bias4",xlab="n/q",ylim = c(-0.08,0.08),main="decaying coefficient with correlated predictors")
data4.2_new[which(abs(bias_new4)>0.05),]
plot(biasforward1,data3.1_new$n.q)
plot(biasforward2,data3.2_new$n.q)
plot(biasforward3,data4.1_new$n.q)
plot(biasforward4,data4.2_new$n.q)

test_con_nocor_sd<- mutate(testdataset_con_nocor_sd,n.q=n/q)
test_con_cor_sd<- mutate(testdataset_con_cor_sd,n.q=n/q)
test_dec_nocor_sd<- mutate(testdataset_dec_nocor_sd,n.q=n/q)
test_dec_cor_sd<- mutate(testdataset_dec_cor_sd,n.q=n/q)

test_pre1<- predict(step_new3.1,newdata=test_con_nocor_sd)
test_pre2<- predict(step_new3.2,newdata=test_con_cor_sd)
test_pre3<- predict(step_new4.1,newdata=test_dec_nocor_sd)
test_pre4<- predict(step_new4.2,newdata=test_dec_cor_sd)

bias_test1<- (test_pre1-log(test_con_nocor_sd$sd))/log(test_con_nocor_sd$sd)
bias_test2<- (test_pre2-log(test_con_cor_sd$sd))/log(test_con_cor_sd$sd)
bias_test3<- (test_pre3-log(test_dec_nocor_sd$sd))/log(test_dec_nocor_sd$sd)
bias_test4<- (test_pre4-log(test_dec_cor_sd$sd))/log(test_dec_cor_sd$sd)

plot(test_con_nocor_sd$n.q,bias_test1,ylab="out-of-sample bias1",xlab="n/q",ylim = c(-0.06,0.06),main="constant coefficient with independent predictors")
plot(test_con_cor_sd$n.q,bias_test2,ylab="out-of-sample bias2",xlab="n/q",ylim = c(-0.06,0.06),main="constant coefficient with correlated predictors")
plot(test_dec_nocor_sd$n.q,bias_test3,ylab="out-of-sample bias3",xlab="n/q",ylim = c(-0.06,0.06),main="decaying coefficient with independent predictors")
plot(test_dec_cor_sd$n.q,bias_test4,ylab="out-of-sample bias4",xlab="n/q",ylim = c(-0.06,0.06),main="decaying coefficient with correlated predictors")

plot(log(test_con_nocor_sd$sd))
plot(log(test_con_cor_sd$sd))
plot(log(test_dec_nocor_sd$sd))
plot(log(test_dec_cor_sd$sd))

test_r1<- cor(test_pre1,log(test_con_nocor_sd$sd))^2
test_r2<- cor(test_pre2,log(test_con_cor_sd$sd))^2
test_r3<- cor(test_pre3,log(test_dec_nocor_sd$sd))^2
test_r4<- cor(test_pre4,log(test_dec_cor_sd$sd))^2

mean(abs(bias_test1))
mean(abs(bias_test2))
mean(abs(bias_test3))
mean(abs(bias_test4))

qqnorm(residuals(step_3.1high))
qqline(residuals(step_3.1high))
qqnorm(residuals(step_3.2high))
qqline(residuals(step_3.2high))
qqnorm(residuals(step_4.1high))
qqline(residuals(step_4.1high))
qqnorm(residuals(step_4.2high))
qqline(residuals(step_4.2high))

qqnorm(residuals(step_new3.1),xlab = "",main = "constant coefficient with independent predictors")
qqline(residuals(step_new3.1))
qqnorm(residuals(step_new3.2),xlab = "",main = "constant coefficient with correlated predictors")
qqline(residuals(step_new3.2))
qqnorm(residuals(step_new4.1),xlab = "",main = "decaying coefficient with independent predictors")
qqline(residuals(step_new4.1))
qqnorm(residuals(step_new4.2),xlab = "",main = "decaying coefficient with correlated predictors")
qqline(residuals(step_new4.2))

install.packages("ppcor")
library(ppcor)
pcor.test(data3.1_new$q)

stats:::logLik

##random forest
#split data
library(randomForest)
library(caret)

set.seed(123)
train<- sample(nrow(sd.lasso), 0.7*nrow(sd.lasso), replace = FALSE)
TrainSet<- sd.lasso[train,]
ValidSet<- sd.lasso[-train,]
train3.1 <- sample(nrow(data3.1_new), 0.7*nrow(data3.1_new), replace = FALSE)
TrainSet3.1 <- data3.1_new[train3.1,c(2,3,4,5,7)]
ValidSet3.1 <- data3.1_new[-train3.1,c(2,3,4,5,7)]
train3.2 <- sample(nrow(data3.2_new), 0.7*nrow(data3.2_new), replace = FALSE)
TrainSet3.2 <- data3.2_new[train3.2,c(2,3,4,5,6,7)]
ValidSet3.2 <- data3.2_new[-train3.2,c(2,3,4,5,6,7)]
train4.1 <- sample(nrow(data4.1_new), 0.7*nrow(data4.1_new), replace = FALSE)
TrainSet4.1 <- data4.1_new[train4.1,c(2,3,4,5,7)]
ValidSet4.1 <- data4.1_new[-train4.1,c(2,3,4,5,7)]
train4.2 <- sample(nrow(data4.2_new), 0.7*nrow(data4.2_new), replace = FALSE)
TrainSet4.2 <- data4.2_new[train4.2,c(2,3,4,5,6,7)]
ValidSet4.2 <- data4.2_new[-train4.2,c(2,3,4,5,6,7)]
#svm
library(e1071)
svm_model4<- svm(TrainSet4.2,kernel = "poly")
summary(svm_model4)
predTrain4 <- predict(svm_model4, TrainSet4.2)
biassvm4<- (predTrain4-TrainSet4.2$true_se)/TrainSet4.2$true_se
plot(biassvm4)
predValid3.1 <- predict(rf3.1, ValidSet3.1)
biasrf_valid3.1<- (predValid3.1-ValidSet3.1$true_se)/ValidSet3.1$true_se
plot(biasrf_valid3.1,ValidSet3.1$n.q)
importance(rf3.1)

#fit default model
rf<- randomForest(true_se ~.,data = TrainSet,importance= TRUE)
rf

rf1<- randomForest(true_se ~.,data = TrainSet,ntree=1000,mtry=6,importance=TRUE)
rf1

predTrain <- predict(rf1, TrainSet)
biasrf<- (predTrain-TrainSet$true_se)/TrainSet$true_se
plot(biasrf)
predValid <- predict(rf1, ValidSet)
biasrf_valid<- (predValid-ValidSet$true_se)/ValidSet$true_se
plot(biasrf_valid)
importance(rf1)
#fit stratified model
for (i in 1:1000) {
  nt<- c(1:1000)
}
rf3.1<- randomForest(true_se ~.,data = TrainSet3.1,ntree=1000,mtry=4,max_features=5
                     ,min_sample_split=70,maxdepth=19,importance=TRUE)
rf3.1
predTrain3.1 <- predict(rf3.1, TrainSet3.1)
biasrf3.1<- (predTrain3.1-TrainSet3.1$true_se)/TrainSet3.1$true_se
mape3.1<- mean(abs(biasrf3.1))
mape3.1
plot(biasrf3.1)
predValid3.1 <- predict(rf3.1, ValidSet3.1)
biasrf_valid3.1<- (predValid3.1-ValidSet3.1$true_se)/ValidSet3.1$true_se
plot(biasrf_valid3.1,ValidSet3.1$n.q)
importance(rf3.1)


rf3.2<- randomForest(true_se ~.,data = TrainSet3.2,ntree=1000,mtry=6,importance=TRUE)
rf3.2
predTrain3.2 <- predict(rf3.2, TrainSet3.2)
biasrf3.2<- (predTrain3.2-TrainSet3.2$true_se)/TrainSet3.2$true_se
plot(biasrf3.2)
predValid3.2 <- predict(rf3.2, ValidSet3.2)
biasrf_valid3.2<- (predValid3.2-ValidSet3.2$true_se)/ValidSet3.2$true_se
plot(biasrf_valid3.2)
importance(rf3.2)



rf4.1<- randomForest(true_se ~.,data = TrainSet4.1,ntree=1000,mtry=4,max_features=5
                     ,min_sample_split=70,maxdepth=20,importance=TRUE)
rf4.1
predTrain4.1 <- predict(rf4.1, TrainSet4.1)
biasrf4.1<- (predTrain4.1-TrainSet4.1$true_se)/TrainSet4.1$true_se
mape4.1<- mean(abs(biasrf4.1))
mape4.1
plot(biasrf4.1)
predValid4.1 <- predict(rf4.1, ValidSet4.1)
biasrf_valid4.1<- (predValid4.1-ValidSet4.1$true_se)/ValidSet4.1$true_se
plot(biasrf_valid4.1,ValidSet4.1$n.q)
importance(rf4.1)

rf3.2<- randomForest(true_se ~.,data = TrainSet3.2,ntree=1000,mtry=6,importance=TRUE)
rf3.2
predTrain3.2 <- predict(rf3.2, TrainSet3.2)
biasrf3.2<- (predTrain3.2-TrainSet3.2$true_se)/TrainSet3.2$true_se
plot(biasrf3.2)
predValid3.2 <- predict(rf3.2, ValidSet3.2)
biasrf_valid3.2<- (predValid3.2-ValidSet3.2$true_se)/ValidSet3.2$true_se
plot(biasrf_valid3.2)
importance(rf3.2)


##take log of the interaction of n and q
varnew1 <- c("K","I(K^2)","I(K^3)","q","I(q^2)","I(q.p)","I(1/q.p^2)"
             ,"I(log(n))","I(log(n)^2)","I(log(n)^3)")
varnew2 <- c("coef","K","I(K^2)","I(K^3)","q","I(q^2)","I(1/q.p)","I(1/q.p^2)"
             ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","rho","I(rho^2)")
varnew3<- c("coef","I(log(n/q))","I(log(n/q)^2)","I(log(n/K))","I(log(n/K)^2)","I(1/q.p)","I(1/q.p^2)")
varnew4<- c("coef","I(log(n/q))","I(log(n/q)^2)","I(log(n/K))","I(log(n/K)^2)","I(1/q.p)","I(1/q.p^2)","rho","I(rho^2)")
exnew3<- cbind(combn(varnew1[2:4],2),combn(varnew1[5:6],2),combn(varnew1[7:8],2),combn(varnew1[9:11],2),
               t(expand.grid(varnew1[2:4],varnew1[9:11])),t(expand.grid(varnew1[5:6],varnew1[9:11])))
exnew3
exnew4<- cbind(combn(varnew2[2:4],2),combn(varnew2[5:6],2),combn(varnew2[7:8],2),combn(varnew2[9:11],2),
               combn(varnew2[12:13],2),t(expand.grid(varnew2[2:4],varnew2[9:11])),t(expand.grid(varnew2[5:6],varnew2[9:11])))
exnew4
exnew5<- cbind(combn(varnew3[2:3],2),combn(varnew3[4:5],2),combn(varnew3[6:7],2))
exnew6<- cbind(combn(varnew4[2:3],2),combn(varnew4[4:5],2),combn(varnew4[6:7],2),combn(varnew4[8:9],2))
internew1<- inter(varnew1,3,exnew3)
internew1
allnew1<- rbind(matrix(varnew1),internew1)
allnew1<- rbind(allnew1,inter(varnew3,3,exnew5))
allnew1<- rbind(matrix(c("I(log(n/q))","I(log(n/q)^2)","I(log(n/K))","I(log(n/K)^2)")),allnew1)
allnew1
internew2<- inter(varnew2,3,exnew4)
allnew2<- rbind(matrix(varnew2),internew2)
allnew2<- rbind(allnew2,inter(varnew4,3,exnew6))
allnew2<- rbind(matrix(c("I(log(n/q))","I(log(n/q)^2)","I(log(n/K))","I(log(n/K)^2)")),allnew2)

bi_new1<- bistep_rss(allnew1,str_fit1,data1)
bi_new2<- bistep(allnew2,str_fit2,data2)
summary(bi_new1)
summary(bi_new2)
plot(predict(bi_new1),residuals(bi_new1))
plot(predict(bi_new2),residuals(bi_new2))
prenew1<- predict(bi_new1,data1[,1:6])
prenew2<- predict(bi_new2,data2[,1:6])
bias_new1<- (prenew1-outcome1)/outcome1
bias_new2<- (prenew2-outcome2)/outcome2
plot(prenew1,bias_new1)
plot(prenew2,bias_new2)
msenew1<- mean((prenew1-outcome1)^2)
msenew2<- mean((prenew2-outcome2)^2)
msenew1
msenew2
qqnorm(residuals(bi_new1))
qqline(residuals(bi_new1))
qqnorm(residuals(bi_new2))
qqline(residuals(bi_new2))

plot(data1$coef,bias_new1)
plot(data1$K,bias_new1)

##stepwise forward fit
str_for<- stepwise_forward(str_fit,var4,3,exclusion4,0.05)
str_for1<- stepwise_forward(str_fit1,var3,3,exclusion3,0.05)
str_for2<- stepwise_forward(str_fit2,var4,3,exclusion4,0.05)
summary(str_for)
summary(str_for1)
summary(str_for2)
extractAIC(str_for1)
extractAIC(str_for2)
plot(predict(str_for),residuals(str_for))
plot(predict(str_for1),residuals(str_for1))
plot(predict(str_for2),residuals(str_for2))
outcome1<- data1[,7]
outcome2<- data2[,7]
pre1<- predict(str_for1,data1[,1:6])
pre2<- predict(str_for2,data2[,1:6])
bias_str1<- (pre1-outcome1)/outcome1
bias_str2<- (pre2-outcome2)/outcome2
mse1<- mean((pre1-outcome1)^2)
mse2<- mean((pre2-outcome2)^2)
mse1
mse2
plot(bias_str1)
plot(bias_str2)
qqnorm(residuals(str_for1))
qqline(residuals(str_for1))
qqnorm(residuals(str_for2))
qqline(residuals(str_for2))
plot(bias_str1,data1$coef)


