library(ggplot2)
library(gridExtra)
setwd("C:/Users/fang035/OneDrive - Wageningen University & Research/Bureaublad/publication/1 se rule/cleaned R file/prediction error estimation")
K<- c(3,5,10,20)
n<-c(100,1000) #n is the sample size
rho<- c(0.01,0.1,0.5,0.9) #correlation between
beta.value<- c(0.2,0.4,0.6,0.8,1.0,1.5,2.0,3.0) #constant beta value
q.pre<- c(10,100,200) #the value of q for standard deviation test
qp.pre<- c(0.1,0.5,10/11) #the ratio between the number of non-zero coefficients to total number of coefficients

scenarios_pre<- expand.grid(q=q.pre,qp=qp.pre,n=n,beta.value=beta.value,K=K)
scenarios_pre2<- expand.grid(q=q.pre,qp=qp.pre,n=n,beta.value=beta.value,rho=rho,K=K)

##constant coefficient with independent predictors
br_con_nocor<- readRDS(file = "br_con_nocor.RDS")
con_nocor_lasso<- br_con_nocor$se_better_lasso$`c(se_better.lasso_1)`/
                  br_con_nocor$min_better_lasso$min_better_lasso
con_nocor_lasso<- cbind(scenarios_pre,"P(se)/P(min)"=con_nocor_lasso)
con_nocor_ridge<- br_con_nocor$se_better_ridge$`c(se_better.ridge_1)`/
                  br_con_nocor$min_better_ridge$`c(min_better.ridge_1)`
con_nocor_ridge<- cbind(scenarios_pre,"P(se)/P(min)"=con_nocor_ridge)

ggplot(data=con_nocor_lasso[con_nocor_lasso$n==100 & con_nocor_lasso$q==10,], 
       mapping = aes(x = beta.value, y = `P(se)/P(min)`,
        group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "a.Lasso prediction error for constant coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), labeller = "label_both")

ggplot(data=con_nocor_ridge[con_nocor_ridge$n==100 & con_nocor_ridge$q==10,], 
       mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "b.Ridge prediction error for constant coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), labeller = "label_both")

ggplot(data=con_nocor_lasso, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
      group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Lasso prediction error for constant coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,2.7)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

ggplot(data=con_nocor_ridge, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
      group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Ridge prediction error for constant coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,2.7)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

##decaying coefficient with independent predictors
br_dec_nocor<- readRDS(file = "br_dec_nocor.RDS")
dec_nocor_lasso<- br_dec_nocor$se_better_lasso2$`c(se_better.lasso_2)`/
                  br_dec_nocor$min_better_lasso2$`c(min_better.lasso_2)`
dec_nocor_lasso<- cbind(scenarios_pre,"P(se)/P(min)"=dec_nocor_lasso)
dec_nocor_ridge<- br_dec_nocor$se_better_ridge2$`c(se_better.ridge_2)`/
                  br_dec_nocor$min_better_ridge2$`c(min_better.ridge_2)`
dec_nocor_ridge<- cbind(scenarios_pre,"P(se)/P(min)"=dec_nocor_ridge)

ggplot(data=dec_nocor_lasso[dec_nocor_lasso$q %in% c(10,200)&dec_nocor_lasso$n==100,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
      group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "a.Lasso prediction error for decaying coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

ggplot(data=dec_nocor_ridge[dec_nocor_ridge$q %in% c(10,200)&dec_nocor_ridge$n==100,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
      group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "b.Ridge prediction error for decaying coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

ggplot(data=dec_nocor_lasso, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Lasso prediction error for decaying coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

ggplot(data=dec_nocor_ridge, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Ridge prediction error for decaying coefficient with independent predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")


##constant coefficient with correlated predictors
br_con_cscor<- readRDS(file = "br_con_cscor.RDS")
con_cscor_lasso<- br_con_cscor$se_better_lasso3$`c(se_better.lasso3)`/
                  br_con_cscor$min_better_lasso3$`c(min_better.lasso3)`
con_cscor_lasso<- cbind(scenarios_pre2,"P(se)/P(min)"=con_cscor_lasso)
con_cscor_ridge<- br_con_cscor$se_better_ridge3$`c(se_better.ridge3)`/
                  br_con_cscor$min_better_ridge3$`c(min_better.ridge3)`
con_cscor_ridge<- cbind(scenarios_pre2,"P(se)/P(min)"=con_cscor_ridge)

ggplot(data=con_cscor_lasso[con_cscor_lasso$q %in% c(10,200)&con_cscor_lasso$rho==0.5,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
                                           group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "a.Lasso prediction error for constant coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,labeller = "label_both")

ggplot(data=con_cscor_ridge[con_cscor_ridge$q %in% c(10,200)&con_cscor_ridge$rho==0.5,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
                                           group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "b.Ridge prediction error for constant coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,labeller = "label_both")


ggplot(data=con_cscor_lasso, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Lasso prediction error for constant coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,labeller = "label_both")

ggplot(data=con_cscor_ridge, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Ridge prediction error for constant coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,labeller = "label_both")



##decaying coefficient with correlated predictors
br_dec_cscor<- readRDS(file = "br_dec_cscor.RDS")
dec_cscor_lasso<- br_dec_cscor$se_better_lasso4$`c(se_better.lasso_4)`/
                  br_dec_cscor$min_better_lasso4$`c(min_better.lasso_4)`
dec_cscor_lasso<- cbind(scenarios_pre2,"P(se)/P(min)"=dec_cscor_lasso)
dec_cscor_ridge<- br_dec_cscor$se_better_ridge4$`c(se_better.ridge_4)`/
                  br_dec_cscor$min_better_ridge4$`c(min_better.ridge_4)`
dec_cscor_ridge<- cbind(scenarios_pre2,"P(se)/P(min)"=dec_cscor_ridge)

ggplot(data=dec_cscor_lasso[dec_cscor_lasso$q %in% c(10,200)&dec_cscor_lasso$rho==0.5,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
                                           group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "a.Lasso prediction error for decaying coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,scale="free",labeller = "label_both")

ggplot(data=dec_cscor_ridge[dec_cscor_ridge$q %in% c(10,200)&dec_cscor_ridge$rho==0.5,], mapping = aes(x = beta.value, y = `P(se)/P(min)`,
                                           group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "b.Ridge prediction error for decaying coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0,3)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,scale="free",labeller = "label_both")


ggplot(data=dec_cscor_lasso, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Lasso prediction error for decaying coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,scale="free",labeller = "label_both")

ggplot(data=dec_cscor_ridge, mapping = aes(x = beta.value, y = `P(se)/P(min)`,
       group=interaction(K,qp,n,q,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  labs(title = "Ridge prediction error for decaying coefficient with correlated predictors",
       x = "beta value",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(vars(n,q,rho), ncol = 4,scale="free",labeller = "label_both")


