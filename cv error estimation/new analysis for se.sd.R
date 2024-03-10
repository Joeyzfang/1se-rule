##cal_se/true_se simulation results analysis
#load data
library(ggplot2)
library(gridExtra)
install.packages("cowplot")
library(cowplot)
setwd("C:/Users/joeyt/Desktop/research_project/HPC/data_to_HPC")
load(file = "setosd_data.RData")
all<- setosd_data[1:3920,]
#split dataset into 4 groups respectively
con<- all[all$coef==0,]
dec<- all[all$coef==1,]
con_noco<- con[con$rho==0,]
con_csco<- con[!con$rho==0,]
dec_noco<- dec[dec$rho==0,]
dec_csco<- dec[!dec$rho==0,]
con_noco<- con_noco[con_noco$n %in% c(100,1000),]
con_noco<- con_noco[con_noco$q %in% c(10,100,200),]
dec_noco<- dec_noco[dec_noco$n %in% c(100,1000),]
dec_noco<- dec_noco[dec_noco$q %in% c(10,100,200),]
con_csco<- con_csco[con_csco$n %in% c(100,1000),]
con_csco<- con_csco[con_csco$q %in% c(10,100,200),]
dec_csco<- dec_csco[dec_csco$n %in% c(100,1000),]
dec_csco<- dec_csco[dec_csco$q %in% c(10,100,200),]

all_se_sd<- rbind(con_csco,con_noco,dec_noco,dec_csco)
plot(all_se_sd[,7], ylab = "SE(CV)/SD(CV)",xlab = "")


#constant coefficient structure with no correlation
noco<- rbind(con_noco,dec_noco)
s<- rep("constant",144)
s[73:144]<- rep("decaying",72)
s
noco<- cbind(noco,s=s)
ggplot(data=noco[noco$n ==100,], mapping = aes(x = q, y = se.sd,
                                                       group=interaction(K,q.p,s),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal()+
  theme(legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))+
  theme(legend.position = "bottom")+
  facet_wrap(~s+n, ncol = 2,labeller = "label_both")


ggplot(data=con_noco[con_noco$n ==100,], mapping = aes(x = q, y = se.sd,
                                                           group=interaction(K,q.p,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal()+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")


a<- ggplot(data=con_noco[con_noco$n ==100,], mapping = aes(x = q, y = se.sd,
          group=interaction(K,q.p,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal()+
  facet_wrap(vars(n,q), ncol = 2,labeller = "label_both")

b<- ggplot(data=dec_noco[dec_noco$n ==100,], mapping = aes(x = q, y = se.sd,
          group=interaction(K,q.p,rho),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal()

grid.arrange(a, b,ncol = 2)

ggplot(data=dec_csco[dec_csco$n ==100,], mapping = aes(x = rho, y = se.sd,
      group=interaction(K,q.p,q),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "correlation coefficient",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))+
  theme(legend.position = "bottom")+
  facet_wrap(~q, ncol = 3,labeller = "label_both")



ggplot(data=con_noco, mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, ncol = 2,labeller = "label_both")

ggplot(data=dec_noco, mapping = aes(x = q, y = se.sd,
                                    group=interaction(K,q.p,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")+
  facet_wrap(~n, ncol = 2,labeller = "label_both")

ggplot(data=dec_csco[dec_csco$rho %in% c(0.1,0.5,0.9),], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p,rho,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))+
  theme(legend.position = "bottom")+
  facet_wrap(vars(rho,n), ncol = 2,labeller = "label_both")

ggplot(data=con_csco[con_csco$rho %in% c(0.1,0.5,0.9),], mapping = aes(x = q, y = se.sd,
        group=interaction(K,q.p,rho,n),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "",
       x = "true predictor",
       y = "SE(CV)/SD(CV)",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1.05)+
  theme_minimal() +
  theme(legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14))+
  theme(legend.position = "bottom")+
  facet_wrap(vars(rho,n), ncol = 2,labeller = "label_both")


p1<- ggplot(data=con_noco[con_noco$n==100,], mapping = aes(x = q, y = se.sd,
                                                      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100",
       x = "true predictor",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p2<- ggplot(data=con_noco[con_noco$n==1000,], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000",
       x = "true predictor",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p1, p2, ncol = 2)


#decaying coefficient structure with no correlation
p3<- ggplot(data=dec_noco[dec_noco$n==100,], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p4<- ggplot(data=dec_noco[dec_noco$n==1000,], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.35,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p3, p4, ncol = 2)


#constant coefficient structure with CS correlation
p5<- ggplot(data=con_csco[con_csco$n==100 & con_csco$rho==0.01,], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100, rho=0.01",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p6<- ggplot(data=con_csco[con_csco$n==100 & con_csco$rho==0.1,], mapping = aes(x = q, y = se.sd,
            group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100, rho=0.1",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p7<- ggplot(data=con_csco[con_csco$n==100 & con_csco$rho==0.5,], mapping = aes(x = q, y = se.sd,
          group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q ,n=100, rho=0.5",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p8<- ggplot(data=con_csco[con_csco$n==100 & con_csco$rho==0.9,], mapping = aes(x = q, y = se.sd,
            group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100, rho=0.9",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p5, p6, ncol = 2)
grid.arrange(p7, p8, ncol = 2)

p9<- ggplot(data=con_csco[con_csco$n==1000 & con_csco$rho==0.01,], mapping = aes(x = q, y = se.sd,
                                                                                group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000, rho=0.01",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p10<- ggplot(data=con_csco[con_csco$n==1000 & con_csco$rho==0.1,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000, rho=0.1",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p11<- ggplot(data=con_csco[con_csco$n==1000 & con_csco$rho==0.5,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q ,n=1000, rho=0.5",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p12<- ggplot(data=con_csco[con_csco$n==1000 & con_csco$rho==0.9,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000, rho=0.9",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.6,1.1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p9, p10, ncol = 2)
grid.arrange(p11, p12, ncol = 2)


#decaying coefficient structure with CS correlation
p13<-ggplot(data=dec_csco[dec_csco$n==100 & dec_csco$rho==0.01,], mapping = aes(x = q, y = se.sd,
      group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100,rho=0.01",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p14<-ggplot(data=dec_csco[dec_csco$n==100 & dec_csco$rho==0.1,], mapping = aes(x = q, y = se.sd,
                                                                                group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100,rho=0.1",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p15<-ggplot(data=dec_csco[dec_csco$n==100 & dec_csco$rho==0.5,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100,rho=0.5",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p16<-ggplot(data=dec_csco[dec_csco$n==100 & dec_csco$rho==0.9,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=100,rho=0.9",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p13, p14, ncol = 2)
grid.arrange(p15, p16, ncol = 2)

p17<-ggplot(data=dec_csco[dec_csco$n==1000 & dec_csco$rho==0.01,], mapping = aes(x = q, y = se.sd,
                                                                                group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000,rho=0.01",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p18<-ggplot(data=dec_csco[dec_csco$n==1000 & dec_csco$rho==0.1,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000,rho=0.1",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p19<-ggplot(data=dec_csco[dec_csco$n==1000 & dec_csco$rho==0.5,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000,rho=0.5",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
p20<-ggplot(data=dec_csco[dec_csco$n==1000 & dec_csco$rho==0.9,], mapping = aes(x = q, y = se.sd,
                                                                               group=interaction(K,q.p),color=factor(K),shape=factor(q.p)))+
  geom_point(size=3)+
  geom_line(size=1)+
  labs(title = "se/sd vs. q,n=1000,rho=0.9",
       x = "q (true predictor)",
       y = "se/sd",
       color = "K",shape = "q/p") +
  scale_color_manual(values = c("blue", "red", "yellow", "green")) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond")) +
  ylim(0.5,1)+
  theme_minimal() +
  theme(legend.position = "bottom")
grid.arrange(p17, p18, ncol = 2)
grid.arrange(p19, p20, ncol = 2)
