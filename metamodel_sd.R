setwd("")
#load train dataset
load(file = "newdata_for_sdlasso.RData")
library(dplyr)
library(rms)
library(MASS)

set.seed(123)

#define a function excluding factors that are not going to contribute to interaction
combn_with_exclusion<- function(X, exclude){
  full <- X
  # remove any columns that have all elements of `exclude`
  full[, !apply(full, 2, function(y) all(exclude %in% y))]
}

#function for interaction building
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


#function for calculate cv mse
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

#define a bidirectional stepwise function
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




###############################################
#data management,data will be divided into four parts in order to fit four piecewise metamodels

data_coef0<- sd.lasso[sd.lasso$coef %in% 0,1:7]
data_coef1<- sd.lasso[sd.lasso$coef %in% 1,1:7]
data_coef0_rho0<- data_coef0[data_coef0$rho %in% 0,1:7]
data_coef0_rho1<- data_coef0[!data_coef0$rho %in% 0,1:7]
data_coef1_rho0<- data_coef1[data_coef1$rho %in% 0,1:7]
data_coef1_rho1<- data_coef1[!data_coef1$rho %in% 0,1:7]

##try to fit with n/q instead of q itself
data_coef0_rho0<- mutate(data_coef0_rho0,n.q=n/q) 
data_coef0_rho1<- mutate(data_coef0_rho1,n.q=n/q)
data_coef1_rho0<- mutate(data_coef1_rho0,n.q=n/q)
data_coef1_rho1<- mutate(data_coef1_rho1,n.q=n/q)


#fit an original model for each piece
str_new1<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3), data = data_coef0_rho0)
str_new2<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2), data = data_coef0_rho1)
str_new3<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3), data = data_coef1_rho0)
str_new4<- lm(log(true_se) ~ rcs(log(n.q),4)+q.p+I(q.p^2)
              +I(log(n))+I(log(n)^2)+I(log(n)^3)+K+I(K^2)+I(K^3)+rho+I(rho^2), data = data_coef1_rho1)

#interaction and variable matrix
var_rho0 <- c("K","I(K^2)","I(K^3)","rcs(log(n.q),4)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","q","I(q^2)","I(q^3)")
var_rho1 <- c("K","I(K^2)","I(K^3)","rcs(log(n.q),4)","q.p","I(q.p^2)"
          ,"I(log(n))","I(log(n)^2)","I(log(n)^3)","rho","I(rho^2)","q","I(q^2)","I(q^3)")

exclusion_rho0<- cbind(combn(var_rho0[1:3],2),combn(var_rho0[5:6],2),combn(var_rho0[7:9],2),combn(var_rho0[10:12],2))
exclusion_rho0
exclusion_rho1<- cbind(combn(var_rho1[1:3],2),combn(var_rho1[5:6],2),combn(var_rho1[7:9],2),combn(var_rho1[10:11],2),combn(var_rho1[12:14],2))
exclusion_rho1
inter_var_rho0<- inter(var_rho0,2,exclusion_rho0)
inter_var_rho0
all_var_rho0<- rbind(matrix(var_rho0),inter_var_rho0)
all_var_rho0
inter_var_rho1<- inter(var_rho1,2,exclusion_rho1)
all_var_rho1<- rbind(matrix(var_rho1),inter_var_rho1)
all_var_rho1


##stepwise bidirection fit
step_new1<- cv.bistep(inter_var_rho0,str_new1,data_coef0_rho0,5)
step_new2<- cv.bistep(all_var_rho1,str_new2,data_coef0_rho1,5)
step_new3<- cv.bistep(inter_var_rho0,str_new3,data_coef1_rho0,5)
step_new4<- cv.bistep(all_var_rho1,str_new4,data_coef1_rho1,5)


summary(step_new1)
summary(step_new2)
summary(step_new3)
summary(step_new4)


coefficients1 <- coef(step_new1)
data1<- data.frame(coefficients1,p_value=summary(step_new1)$coefficients[, "Pr(>|t|)"])
data1<- mutate(parameter=rownames(data1),data1)
data1<- cbind(parameter=data1$parameter,coefficient=data1$coefficients1,p_value=data1$p_value)
coefficients2 <- coef(step_new2)
data2<- data.frame(coefficients2,p_value=summary(step_new2)$coefficients[, "Pr(>|t|)"])
data2<- mutate(parameter=rownames(data2),data2)
data2<- cbind(parameter=data2$parameter,coefficient=data2$coefficients2,p_value=data2$p_value)
coefficients3 <- coef(step_new3)
data3<- data.frame(coefficients3,p_value=summary(step_new3)$coefficients[, "Pr(>|t|)"])
data3<- mutate(parameter=rownames(data3),data3)
data3<- cbind(parameter=data3$parameter,coefficient=data3$coefficients3,p_value=data3$p_value)
coefficients4 <- coef(step_new4)
data4<- data.frame(coefficients4,p_value=summary(step_new4)$coefficients[, "Pr(>|t|)"])
data4<- mutate(parameter=rownames(data4),data4)
data4<- cbind(parameter=data4$parameter,coefficient=data4$coefficients4,p_value=data4$p_value)

write.csv(data1,file = "con_nocor_model.csv")
write.csv(data2,file = "con_cor_model.csv")
write.csv(data3,file = "dec_nocor_model.csv")
write.csv(data4,file = "dec_cor_model.csv")

cex_size <- 1.5  # Adjust the scaling factor as per your preference
par(cex.lab = cex_size,cex.main=cex_size)
plot(predict(step_new1),residuals(step_new1),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="constant coefficient with independent predictors")
plot(predict(step_new2),residuals(step_new2),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="constant coefficient with correlated predictors")
plot(predict(step_new3),residuals(step_new3),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="decaying coefficient with independent predictors")
plot(predict(step_new4),residuals(step_new4),xlab="predicted value",ylab="residuals",ylim=c(-0.21,0.21),main="decaying coefficient with correlated predictors")



prenew1<- predict(step_new1)
prenew2<- predict(step_new2)
prenew3<- predict(step_new3)
prenew4<- predict(step_new4)



outcome1<- log(data_coef0_rho0[,7])
outcome2<- log(data_coef0_rho1[,7])
outcome3<- log(data_coef1_rho0[,7])
outcome4<- log(data_coef1_rho1[,7])


bias_new1<- (prenew1-outcome1)/outcome1
bias_new2<- (prenew2-outcome2)/outcome2
bias_new3<- (prenew3-outcome3)/outcome3
bias_new4<- (prenew4-outcome4)/outcome4



par(mfrow=c(2,2))
plot(data_coef0_rho0$n.q,bias_new1,ylab="in-sample bias1",xlab="n/q",ylim = c(-0.08,0.08),main="constant coefficient with independent predictors")
plot(data_coef0_rho1$n.q,bias_new2,ylab="in-sample bias2",xlab="n/q",ylim = c(-0.08,0.08),main="constant coefficient with correlated predictors")
plot(data_coef1_rho0$n.q,bias_new3,ylab="in-sample bias3",xlab="n/q",ylim = c(-0.08,0.08),main="decaying coefficient with independent predictors")
plot(data_coef1_rho1$n.q,bias_new4,ylab="in-sample bias4",xlab="n/q",ylim = c(-0.08,0.08),main="decaying coefficient with correlated predictors")
par(mfrow=c(1,1))

#read test dataset
testdataset_con_nocor_sd<- read.csv(file = "testdataset_con_nocor_sd.csv")
testdataset_con_cor_sd<- read.csv(file = "testdataset_con_cor_sd.csv")
testdataset_dec_nocor_sd<- read.csv(file = "testdataset_dec_nocor_sd.csv")
testdataset_dec_cor_sd<- read.csv(file = "testdataset_dec_cor_sd.csv")

test_con_nocor_sd<- mutate(testdataset_con_nocor_sd,n.q=n/q)
test_con_cor_sd<- mutate(testdataset_con_cor_sd,n.q=n/q)
test_dec_nocor_sd<- mutate(testdataset_dec_nocor_sd,n.q=n/q)
test_dec_cor_sd<- mutate(testdataset_dec_cor_sd,n.q=n/q)

test_pre1<- predict(step_new1,newdata=test_con_nocor_sd)
test_pre2<- predict(step_new2,newdata=test_con_cor_sd)
test_pre3<- predict(step_new3,newdata=test_dec_nocor_sd)
test_pre4<- predict(step_new4,newdata=test_dec_cor_sd)

bias_test1<- (test_pre1-log(test_con_nocor_sd$sd))/log(test_con_nocor_sd$sd)
bias_test2<- (test_pre2-log(test_con_cor_sd$sd))/log(test_con_cor_sd$sd)
bias_test3<- (test_pre3-log(test_dec_nocor_sd$sd))/log(test_dec_nocor_sd$sd)
bias_test4<- (test_pre4-log(test_dec_cor_sd$sd))/log(test_dec_cor_sd$sd)

plot(test_con_nocor_sd$n.q,bias_test1,ylab="out-of-sample bias1",xlab="n/q",ylim = c(-0.06,0.06),main="constant coefficient with independent predictors")
plot(test_con_cor_sd$n.q,bias_test2,ylab="out-of-sample bias2",xlab="n/q",ylim = c(-0.06,0.06),main="constant coefficient with correlated predictors")
plot(test_dec_nocor_sd$n.q,bias_test3,ylab="out-of-sample bias3",xlab="n/q",ylim = c(-0.06,0.06),main="decaying coefficient with independent predictors")
plot(test_dec_cor_sd$n.q,bias_test4,ylab="out-of-sample bias4",xlab="n/q",ylim = c(-0.06,0.06),main="decaying coefficient with correlated predictors")


#MAPE calculation
mean(abs(bias_test1))
mean(abs(bias_test2))
mean(abs(bias_test3))
mean(abs(bias_test4))


qqnorm(residuals(step_new1),xlab = "",main = "constant coefficient with independent predictors")
qqline(residuals(step_new1))
qqnorm(residuals(step_new2),xlab = "",main = "constant coefficient with correlated predictors")
qqline(residuals(step_new2))
qqnorm(residuals(step_new3),xlab = "",main = "decaying coefficient with independent predictors")
qqline(residuals(step_new3))
qqnorm(residuals(step_new4),xlab = "",main = "decaying coefficient with correlated predictors")
qqline(residuals(step_new4))


