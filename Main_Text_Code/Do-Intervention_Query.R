##in the discrete case, do-intervention query
#set random seed
set.seed(65432)

#load R packages required
library(gRbase)
library(graph)
library(BiocGenerics)
library(grid)
library(Rgraphviz)
suppressPackageStartupMessages(library(bnlearn))
suppressPackageStartupMessages(library(gRain))
suppressPackageStartupMessages(library(openxlsx2))
suppressPackageStartupMessages(library(ggplot2))
library(ggforce)

#some preparing assignments
T_cpt<-vector()
C_cpt1<-vector()
C_cpt2<-vector()
C_cpt3<-vector()
C_cpt4<-vector()

T_cpt1_d<-vector()

IC_cpt<-vector()
c_a1<-vector()
c_a2<-vector()
c_a3<-vector()
c_a4<-vector()

c_3<-vector()
c_03<-vector()
t<-100


#given the Bayesian network (BN) structure
asianet <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung][xray|either][dysp|bronc:either]")
asianet4 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung]")
asianet2a <- model2network("[asia][smoke|tub][tub|asia][lung|smoke:tub][bronc|smoke:either][either|tub:lung:smoke]")
asianet2b <- model2network("[asia][smoke|asia][tub|asia:smoke][lung|smoke:tub][bronc|either:smoke][either|lung:smoke:tub]")
asianet2c <- model2network("[asia][smoke|tub][tub|asia][lung|bronc:smoke:tub][bronc|smoke:tub][either|bronc:lung:tub]")
asianet2d <- model2network("[asia][smoke|asia][tub|asia:bronc:lung][lung|asia:bronc:smoke][bronc|asia:smoke][either|bronc:lung:tub]")

yn <- c("yes", "no")

#generating discrete distributions compatitive with the BN structure
for (i in 1:t) {
  a<-runif(18,0.000001,0.999999)
  #a<-runif(18,0.01,0.99)
  cptA <- matrix(c(a[1], 1-a[1]), ncol = 2, dimnames = list(NULL, yn))
  cptS <- matrix(c(a[2], 1-a[2]), ncol = 2, dimnames = list(NULL, yn))
  cptT <- matrix(c(a[3], 1-a[3], a[4], 1-a[4]), ncol = 2, dimnames = list(tub = yn, 
                                                                          asia = yn))
  cptL <- matrix(c(a[5],1-a[5],a[6], 1-a[6]), ncol = 2, dimnames = list(lung = yn, 
                                                                        smoke = yn))
  cptB <- matrix(c(a[7], 1-a[7], a[8], 1-a[8]), ncol = 2, dimnames = list(bronc = yn, 
                                                                          smoke = yn))
  cptE <- c(a[9], 1-a[9], a[10], 1-a[10], a[11], 1-a[11], a[12], 1-a[12])
  dim(cptE) <- c(2, 2, 2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn)
  cptX <- matrix(c(a[13], 1-a[13], a[14], 1-a[14]), ncol = 2, dimnames = list(xray = yn, 
                                                                              either = yn))
  cptD <- c(a[15], 1-a[15], a[16], 1-a[16], a[17], 1-a[17], a[18], 1-a[18])
  dim(cptD) <- c(2, 2, 2)
  dimnames(cptD) <- list(dysp = yn, either = yn, bronc = yn)
  asianet.fit <- custom.fit(asianet, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                 lung = cptL, bronc = cptB, either = cptE, xray = cptX, dysp = cptD))
  
  fitted_mutilated <- mutilated(asianet.fit, evidence = list(tub = "no"))
  asia.jtree <- compile(as.grain(asianet.fit))
  asia.jtree1 <- compile(as.grain(fitted_mutilated))
  asia.ev1<- setEvidence(asia.jtree1, nodes = c("xray", "dysp"), states = c("yes", "yes"))                                                                         
  
  cp2<-as.data.frame(querygrain(asia.jtree1,nodes = "asia", type ='marginal'))
  cp2a<-as.data.frame(querygrain(asia.ev1,nodes = "asia", type ='marginal'))
  cp3<-as.data.frame(querygrain(asia.jtree1,nodes = "smoke", type ='marginal'))
  cp3a<-as.data.frame(querygrain(asia.ev1,nodes = "smoke", type ='marginal'))
  cp3b<-as.data.frame(querygrain(asia.ev1,nodes ="smoke", type ='marginal'))
  
  cp3c<-as.data.frame(querygrain(asia.ev1,nodes = "smoke", type ='marginal'))
  
  cp5<-as.data.frame(querygrain(asia.jtree1,nodes = c("lung","smoke"), type ='conditional'))
  cp5a<-as.data.frame(querygrain(asia.ev1,nodes = c("lung","smoke"), type ='conditional'))#l|s,t
  cp5b<-as.data.frame(querygrain(asia.ev1,nodes = c("lung","smoke"), type ='conditional'))
  

  cp5d<-as.data.frame(querygrain(asia.ev1,nodes = c("lung","bronc","smoke"), type ='conditional'))#l|b,t,s
  cp5e<-as.data.frame(querygrain(asia.ev1,nodes = c("lung","bronc","smoke"), type ='conditional'))#l|a,s,b
  
  
  cp6<-as.data.frame(querygrain(asia.jtree1,nodes = c("either","lung"), type ='conditional'))#e|t,l
  cp6a<-as.data.frame(querygrain(asia.ev1,nodes = c("either","lung","smoke"), type ='conditional'))#e|l,s
  cp6b<-as.data.frame(querygrain(asia.ev1,nodes = c("either","lung"), type ='conditional'))#e|l
  
  cp6d<-as.data.frame(querygrain(asia.ev1,nodes = c("either","bronc","lung"), type ='conditional'))#e|b,t,l
  
  cp7<-as.data.frame(querygrain(asia.jtree1,nodes = c("bronc","smoke"), type ='conditional'))
  cp7a<-as.data.frame(querygrain(asia.ev1,nodes = c("bronc","smoke","either"), type ='conditional'))#b|s,e
  cp7b<-as.data.frame(querygrain(asia.ev1,nodes = c("bronc","smoke"), type ='conditional'))
  
  cp7d<-as.data.frame(querygrain(asia.ev1,nodes = c("bronc","smoke"), type ='conditional'))#b|s,t
  cp7e<-as.data.frame(querygrain(asia.ev1,nodes = c("bronc","smoke"), type ='conditional'))#b|s,a
  
  cp8<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","either"), type ='conditional'))
  cp9<-as.data.frame(querygrain(asia.jtree1,nodes = c("dysp","either","bronc"), type ='conditional'))#d|e,b
  
  #true values calculated according to the factorization formula of the joint distributions
  T_cp<-as.data.frame(querygrain(asia.ev1,nodes = c("asia","smoke","lung","either","bronc"), type ='joint'))#s,b,a,e,l
  cp12<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","dysp"), type ='joint'))#s,b,a,e,l
  T_cpt1_d[i]<-T_cp$no.yes.yes.no[2]
  T_cpt[i]<-cp2$asia[2]*1*cp3$smoke[2]*cp5$no[1]*cp6$yes[1]*cp7$no[2]*cp8$yes[1]*cp9$yes.no[1]/cp12$yes[1]#(cp10$yes.no[1]*cp11$no[1])
  #values calculated according to the factorization formula of the conditional distributions (by our theory results)
  C_cpt1[i]<-cp2a$asia[2]*1*cp3a$smoke[2]*cp5a$no[1]*cp6a$yes.no[1]*cp7a$no.yes[2]
  C_cpt2[i]<-cp2a$asia[2]*1*cp3c$smoke[2]*cp5a$no[1]*cp6a$yes.no[1]*cp7a$no.yes[2]
  C_cpt3[i]<-cp2a$asia[2]*1*cp3a$smoke[2]*cp5d$no.no[1]*cp6d$no.yes[1]*cp7d$no[2]
  C_cpt4[i]<-cp2a$asia[2]*1*cp3c$smoke[2]*cp5e$no.no[1]*cp6d$no.yes[1]*cp7e$no[2]
  #values calculated according to the factorization formula of the conditional distributions, ignring the changes of the structure of the conditional distributions
  IC_cpt[i]<-cp2a$asia[2]*1*cp3b$smoke[2]*cp5b$no[1]*cp6b$yes[1]*cp7b$no[2]
}

for (j in 1:t) {
  c_a1[j]<-(C_cpt1[j]-T_cpt[j])/T_cpt[j]#relative bias of C_Method_1
  c_a2[j]<-(C_cpt2[j]-T_cpt[j])/T_cpt[j]#relative bias of C_Method_2
  c_a3[j]<-(C_cpt3[j]-T_cpt[j])/T_cpt[j]#relative bias of C_Method_3
  c_a4[j]<-(C_cpt4[j]-T_cpt[j])/T_cpt[j]#relative bias of C_Method_4
  
  c_3[j]<-(IC_cpt[j]-T_cpt[j])/T_cpt[j]#relative bias of IC_Method

  
  c_03[j]<-IC_cpt[j]-T_cpt[j]#bias of IC_Method
}

#save the true values and the values obtained by C_Method and IC_Method
s1<-cbind(T_cpt,C_cpt1,C_cpt2,C_cpt3,IC_cpt)

#save results of relative biases and biases of C_Method and IC_Method
s2<-cbind(c_03,c_a1,c_a2,c_a3,c_a4,c_3)


#the scatter plot for the relative biases of C_Method and IC_Method
sample<-1:t
type <- rep(c('C-Method_1','C-Method_2','C-Method_3','C-Method_4','IC-Method'),each = t)
value <- c(c_a1,c_a2,c_a3,c_a4,c_3)
df <- data.frame(sample =sample, type = type, value = value)
Method<-factor(type)

p1<-ggplot(df,mapping = aes(x=sample,y=value,shape=Method,color=Method))+
  geom_point(size=6,alpha=1)+
  theme_bw()+ 
  xlab("Indexes of 100 Different Discrete Probability Distributions")+
  ylab("Relative Bias")+
  scale_shape_manual(values = c(3,4,5,6,20))+
  facet_zoom(ylim=c(-1,1))+ 
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(family = 'serif'),
        legend.position = 'top',
        legend.direction = 'horizontal')
p1
