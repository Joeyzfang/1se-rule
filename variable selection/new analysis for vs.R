library(ggplot2)
library(gridExtra)
library(matrixStats)
#unzip vs_F1.zip and load data from it
setwd(".../vs_F1")
load(file = "scenarios_vs_ind.RData")
load(file = "scenarios_vs_dep.RData")
load(file = "F1_dec_cor_min.RData")
load(file = "F1_dec_cor_1se.RData")
load(file = "F1_dec_nocor_min.RData")
load(file = "F1_dec_nocor_1se.RData")
load(file = "F1_con_cor_min.RData")
load(file = "F1_con_cor_1se.RData")
load(file = "F1_con_nocor_min.RData")
load(file = "F1_con_nocor_1se.RData")


par(mfrow = c(2, 2))
par(mfrow = c(1, 1))
##F1 score analysis
#con_nocor
mean_F1_con_nocor_min<- rowMeans(F1_con_nocor_min,na.rm = TRUE)
rowSds(F1_con_nocor_min,na.rm = TRUE)
mean_F1_con_nocor_1se<- rowMeans(F1_con_nocor_1se,na.rm = TRUE)
ratio_con_nocor<- mean_F1_con_nocor_min/mean_F1_con_nocor_1se
ratio_con_nocor<- cbind(scenarios_vs_ind,ratio_con_nocor)
plot(ratio_con_nocor$ratio_con_nocor,xlab="",ylab = "F1(min)/F1(se)",main="constant_independent")
ggplot(data=ratio_con_nocor, mapping = aes(x = q, y = `ratio_con_nocor`,
      group=interaction(K,qp,n),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  labs(x = "true predictor",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, scales = "fixed",labeller = "label_both")

#con_cor
mean_F1_con_cor_min<- rowMeans(F1_con_cor_min,na.rm = TRUE)
mean_F1_con_cor_1se<- rowMeans(F1_con_cor_1se,na.rm = TRUE)
ratio_con_cor<- mean_F1_con_cor_min/mean_F1_con_cor_1se
ratio_con_cor<- cbind(scenarios_vs_dep,ratio_con_cor)
plot(ratio_con_cor$ratio_con_cor,xlab = "",ylab = "F1(min)/F1(se)",main="constant_correlated")
ggplot(data=ratio_con_cor, mapping = aes(x = q, y = `ratio_con_cor`,
        group=interaction(K,qp,n,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  labs(x = "true predictor",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n+rho, scales = "fixed",labeller = "label_both")

#dec_nocor
mean_F1_dec_nocor_min<- rowMeans(F1_dec_nocor_min,na.rm = TRUE)
mean_F1_dec_nocor_1se<- rowMeans(F1_dec_nocor_1se,na.rm = TRUE)
ratio_dec_nocor<- mean_F1_dec_nocor_min/mean_F1_dec_nocor_1se
ratio_dec_nocor<- cbind(scenarios_vs_ind,ratio_dec_nocor)
plot(ratio_dec_nocor$ratio_dec_nocor,xlab = "",ylab = "F1(min)/F1(se)",main="decaying_independent")
ggplot(data=ratio_dec_nocor, mapping = aes(x = q, y = `ratio_dec_nocor`,
                                           group=interaction(K,qp,n),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  labs(x = "true predictor",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, scales = "fixed",labeller = "label_both")

#dec_cor
mean_F1_dec_cor_min<- rowMeans(F1_dec_cor_min,na.rm = TRUE)
mean_F1_dec_cor_1se<- rowMeans(F1_dec_cor_1se,na.rm = TRUE)
ratio_dec_cor<- mean_F1_dec_cor_min/mean_F1_dec_cor_1se
ratio_dec_cor<- cbind(scenarios_vs_dep,ratio_dec_cor)
plot(ratio_dec_cor$ratio_dec_cor,xlab = "",ylab = "F1(min)/F1(se)",main="decaying_correlated")
ggplot(data=ratio_dec_cor[ratio_dec_cor$n==100,], mapping = aes(x = q, y = `ratio_dec_cor`,
                                         group=interaction(K,qp,n,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  labs(title = "a.F1 ratio for data with decaying coefficient and correlated predictors",
    y="F1(min)/F1(se)",
    x = "true predictor",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.9,2.14)+
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))+
  theme(legend.position = "bottom")+
  facet_wrap(~n+rho, scales = "fixed",labeller = "label_both")



##constant coefficient with independent predictors
moretrue_con_nocor<- readRDS(file = "con_nocor_moretruepre_new.RDS")
cex_size <- 1.5  # Adjust the scaling factor as per your preference
par(cex.lab = cex_size,cex.main=cex_size)
plot(moretrue_dec_nocor[which(moretrue_dec_nocor$beta.value==1),6],xlab = "",ylab = "P（se)-P(min)",ylim=c(-1,0),main="decaying_independent")
lessfake_con_nocor<- readRDS(file = "con_nocor_lessfakepre_new.RDS")
ggplot(data=moretrue_con_nocor[moretrue_con_nocor$n==100,], mapping = aes(x = q, y = `P_se-P_min`,
                group=interaction(K,qp,beta.value),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~beta.value, scales = "free")
ggplot(data=moretrue_con_nocor[moretrue_con_nocor$n==1000,], mapping = aes(x = q, y = `P_se-P_min`,
                                                                          group=interaction(K,qp,beta.value),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~beta.value, scales = "free")
ggplot(data=moretrue_con_nocor[moretrue_con_nocor$beta.value==1,], mapping = aes(x = q, y = `P_se-P_min`,
      group=interaction(K,qp,n),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, scales = "free")



moretrue_dec_nocor<- readRDS(file = "dec_nocor_moretruepre_new.RDS")
plot(moretrue_dec_nocor[which(moretrue_dec_nocor$beta.value==1),6],xlab = "",ylab = "P(se)-P(min)",main="decaying_independent")
ggplot(data=moretrue_dec_nocor[moretrue_dec_nocor$beta.value==1,], mapping = aes(x = q, y = `P_se-P_min`,
                               group=interaction(K,qp,n),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, scales = "free")

moretrue_con_cscor<- readRDS(file = "con_cscor_moretruepre_new.RDS")
moretrue_dec_cscor<- readRDS(file = "dec_cscor_moretruepre_new.RDS")

plot(moretrue_con_cscor[which(moretrue_con_cscor$beta.value==1&moretrue_con_cscor$rho %in% c(0.1,0.9)),7],xlab = "",ylab = "P(se)-P(min)",ylim=c(-1,0),main="constant_correlated")
plot(moretrue_dec_cscor[which(moretrue_dec_cscor$beta.value==1&moretrue_dec_cscor$rho %in% c(0.1,0.9)),7],xlab = "",ylab = "P(se)-P(min)",ylim=c(-1,0),main="decaying_correlated")


s1<- rep("constant",288)
s2<- rep("decaying",288)
s3<- rep("constant",864)
s4<- rep("decaying",864)
s<- rbind(as.matrix(s1),as.matrix(s2),as.matrix(s3),as.matrix(s4))
moretrue_con_nocor<- cbind(moretrue_con_nocor[,1:4],rho=rep(0,288),
                           moretrue_con_nocor[,5:6])
moretrue_dec_nocor<- cbind(moretrue_dec_nocor[,1:4],rho=rep(0,288),
                           moretrue_dec_nocor[,5:6])
all_vs<- rbind(moretrue_con_nocor,moretrue_dec_nocor,
               moretrue_con_cscor,moretrue_dec_cscor)
all_vs<- cbind(s,all_vs)
saveRDS(all_vs,file = "all_scenarios_vs.RDS")

ggplot(data=all_vs[all_vs$beta.value==1 & all_vs$rho==0,], mapping = aes(x = q, y = `P_se-P_min`,
      group=interaction(K,qp,n),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "variable selesction for independent predictors",
       x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(s~n, ncol = 2,scales = "free")

ggplot(data=all_vs[all_vs$beta.value==1 & all_vs$rho %in% c(0.1,0.9),], mapping = aes(x = q, y = `P_se-P_min`,
       group=interaction(K,qp,n,rho),color=factor(K),shape=factor(qp)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "variable selesction for correlated predictors",
       x = "q (true predictor)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(-1,0)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(s~n+rho, ncol = 2,scales = "free")
