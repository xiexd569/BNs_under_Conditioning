##in the discrete case, do probability reasoning after related local conditional probabilities are estimated by simulation data
##including the runtime of J_Method_cd, J_Method_cr,J_Method_fr

#set random seed
set.seed(65432)

#load R packages required
library(ggalt)
library(gRbase)
library(BiocGenerics)
library(graph)
library(grid)
library(Rgraphviz)
suppressPackageStartupMessages(library(bnlearn))
suppressPackageStartupMessages(library(gRain))
suppressPackageStartupMessages(library(openxlsx2))
suppressPackageStartupMessages(library(ggplot2))
library(ggforce)
library(igraph)



#some preparing assignments
T_cpt<-vector()


t<-1
t1<-100
s<-7

c_02<-vector()
c_04<-vector()
c_06<-vector()

c_0a<-matrix(rep(0,s*t),nrow=t)
c_0b<-matrix(rep(0,s*t),nrow=t)
c_04a<-matrix(rep(0,s*t),nrow=t)


c_12<-vector()
c_14<-vector()
c_16<-vector()

c_1a<-matrix(rep(0,s*t),nrow=t)
c_1b<-matrix(rep(0,s*t),nrow=t)
c_14a<-matrix(rep(0,s*t),nrow=t)


c_0ae<-vector()
c_0be<-vector()
c_02e<-vector()

c_1ae<-vector()
c_1be<-vector()
c_12e<-vector()



elapse_J<-vector()
elapse_J_p<-vector()
elapse_J_pt<-vector()
elapse_J_t<-vector()
elapse_J_a<-vector()


elapse_J_cliq<-vector()
elapse_J_cliq_a<-vector()

elapse_J_t1<-matrix(rep(0,s*t),nrow=t)
elapse_J_pt1<-matrix(rep(0,s*t),nrow=t)
elapse_J_cliq1<-matrix(rep(0,s*t),nrow=t)

elapse_C1<-vector()
elapse_C1_t<-vector()
elapse_C1_t1<-matrix(rep(0,s*t),nrow=t)
elapse_C1_pt1<-matrix(rep(0,s*t),nrow=t)



elapse_clique_infer1<-vector()
elapse_clique_infer2<-matrix(rep(0,s*t),nrow=t)
elapse_clique_infer<-vector()

asianet.fit<-list()

J_cpt<-vector()
J_cpt2<-vector()
J_cpt3<-vector()

S<-c(asia="no",tub="yes",either="yes",lung="yes",bronc="no",smoke="no",xray="yes",dysp="yes")
Z<-c(100,500,1000,2500, 5000,7500,10000)


#given the Bayesian network (BN) structures
asianet1 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung][xray|either][dysp|bronc:either]")
asianet4 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung]")
asianet2a <- model2network("[asia][smoke|tub][tub|asia][lung|smoke:tub][bronc|smoke:either][either|tub:lung:smoke]")
asianet2b <- model2network("[asia][smoke|asia][tub|asia:smoke][lung|smoke:tub][bronc|either:smoke][either|lung:smoke:tub]")
asianet2c <- model2network("[asia][smoke|tub][tub|asia][lung|bronc:smoke:tub][bronc|smoke:tub][either|bronc:lung:tub]")
asianet2d <- model2network("[asia][smoke|asia][tub|asia:bronc:lung][lung|asia:bronc:smoke][bronc|asia:smoke][either|bronc:lung:tub]")
yn <- c("yes", "no")

#generating discrete distributions compatitive with the BN structure 1, from which data are generated
for (i in 1:t) {
  a<-runif(18,0.000001,0.999999)
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
  asianet.fit[[i]]<- custom.fit(asianet1, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                      lung = cptL, bronc = cptB, either = cptE, xray = cptX, dysp = cptD))
  
  asia.jtree <- compile(as.grain(asianet.fit[[i]]))
  asia.ev1<- setEvidence(asia.jtree, nodes = c("xray", "dysp"), states = c("yes", "yes")) 
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia"), type ='conditional')
  cptT <- matrix(c(a2[1], 1-a2[1], a2[3], 1-a2[3]), ncol = 2, dimnames = list(tub = yn, 
                                                                              asia = yn))
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung"), type ='conditional')
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7])
  dim(cptE) <- c(2, 2, 2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn)
  
  a4<-querygrain(asia.ev1,nodes = "smoke", type ='marginal')
  cptS<-matrix(c(a4$smoke[1], 1-a4$smoke[1]), ncol = 2, dimnames = list(NULL, yn))
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke"), type ='conditional')
  cptL <- matrix(c(a5[1],1-a5[1],a5[3], 1-a5[3]), ncol = 2, dimnames = list(lung = yn, 
                                                                            smoke = yn))
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke"), type ='conditional') 
  cptB <- matrix(c(a6[1], 1-a6[1], a6[3], 1-a6[3]), ncol = 2, dimnames = list(bronc = yn,smoke = yn))
}
  

#the true value
for (i in 1:t) {
  asia.jtree <- compile(as.grain(asianet.fit[[i]]))
  #asia.ev1<- setFinding(asia.jtree, nodes = c("xray", "dysp"), states = c("yes", "yes")) 
  
  cp2<-as.data.frame(querygrain(asia.jtree,nodes = "asia", type ='marginal'))
  cp3<-as.data.frame(querygrain(asia.jtree,nodes = "smoke", type ='marginal'))
  cp4<-as.data.frame(querygrain(asia.jtree,nodes = c("tub","asia"), type ='conditional'))
  cp5<-as.data.frame(querygrain(asia.jtree,nodes = c("lung","smoke"), type ='conditional'))
  cp6<-as.data.frame(querygrain(asia.jtree,nodes = c("either","tub","lung"), type ='conditional'))
  cp7<-as.data.frame(querygrain(asia.jtree,nodes = c("bronc","smoke"), type ='conditional'))
  cp8<-as.data.frame(querygrain(asia.jtree,nodes = c("xray","either"), type ='conditional'))
  cp9<-as.data.frame(querygrain(asia.jtree,nodes = c("dysp","either","bronc"), type ='conditional'))
  cp10<-as.data.frame(querygrain(asia.jtree,nodes = c("xray","dysp"), type ='joint'))
  
  #true values calculated according to the factorization formula of the joint distribution
  T_cpt[i]<-cp2$asia[2]*cp3$smoke[2]*cp4$no[1]*cp5$no[1]*cp6$yes.yes[1]*cp7$no[2]*cp8$yes[1]*cp9$yes.no[1]/cp10$yes[1]
}


#runtime of three kinds of J_Methods
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    for (k in 1:t1) {
      data1<-rbn(asianet.fit[[i]],n)
      
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes1 <- bn.fit(asianet1, data = data1, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_J_p[k] <- difftime(end1,start1,units="secs")
      
      start11 <- Sys.time()
      asia.jtree1<- compile(as.grain(bn.bayes1))
      end11 <- Sys.time()
      elapse_J_cliq[k] <- difftime(end11,start11,units="secs")
      
      
      start<-Sys.time()
      b22<-as.data.frame(querygrain(asia.jtree1,nodes = c("asia","tub","smoke","lung","either","bronc","xray","dysp"), type ='joint'))
      b10<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","dysp"), type ='joint'))
      J_cpt2[k]<-b22$no.yes.yes.yes.yes.no.yes[2]/b10$yes[1]
      end<-Sys.time()
      elapse_clique_infer1[k] <- difftime(end,start,units="secs")
      
      start <- Sys.time()
      #values calculated according to the factorization formula of the local CPTs obtained by of the data from the joint distributions
      asia.ev1<- setEvidence(asia.jtree1, nodes = c("xray", "dysp"), states = c("yes", "yes")) 
      z11<-as.data.frame(querygrain(asia.ev1,nodes = c("asia","tub","smoke","lung","either","bronc"), type ='joint'))#先最内，再外排
      J_cpt[k]<-z11$no.yes.yes.no.no[2]
      end <- Sys.time()
      elapse_J[k] <- difftime(end,start,units="secs")
      
      start<-Sys.time()
      b1<-as.vector(bn.bayes1$asia$prob)
      b2<-as.data.frame(bn.bayes1$bronc$prob)
      b3<-as.data.frame(bn.bayes1$dysp$prob)
      b4<-as.data.frame(bn.bayes1$either$prob)
      b5<-as.data.frame(bn.bayes1$lung$prob)
      b6<-as.data.frame(bn.bayes1$smoke$prob)
      b7<-as.data.frame(bn.bayes1$tub$prob)
      b8<-as.data.frame(bn.bayes1$xray$prob)
      b10<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","dysp"), type ='joint'))
      J_cpt3[k]<-b1[2]*b2$Freq[4]*b3$Freq[3]*b4$Freq[1]*b5$Freq[3]*b6$Freq[2]*b7$Freq[3]*b8$Freq[1]/b10$yes[1]
      end<-Sys.time()
      
      elapse_C1[k]<-difftime(end,start,units="secs")
      
      c_02[k]<-J_cpt[k]-T_cpt[i]
      c_12[k]<-(J_cpt[k]-T_cpt[i])^2
      
      c_04[k]<-J_cpt2[k]-T_cpt[i]
      c_14[k]<-(J_cpt2[k]-T_cpt[i])^2
      
      c_06[k]<-J_cpt3[k]-T_cpt[i]
      c_16[k]<-(J_cpt3[k]-T_cpt[i])^2
      
    }
    c_0a[i,m]<-mean(c_02)#Bias
    c_1a[i,m]<-sqrt(mean(c_12))#RMSE
    
    c_0b[i,m]<-mean(c_04)#Bias
    c_1b[i,m]<-sqrt(mean(c_14))#RMSE
    
    c_04a[i,m]<-mean(c_06)#Bias
    c_14a[i,m]<-sqrt(mean(c_16))#RMSE
    
    elapse_J_cliq1[i,m]<-sum(elapse_J_cliq)
    elapse_J_t1[i,m] <- sum(elapse_J)
    elapse_J_pt1[i,m] <- sum(elapse_J_p)
    elapse_clique_infer2[i,m]<-sum(elapse_clique_infer1)
    elapse_C1_t1[i,m]<-sum(elapse_C1)
  }
  end2<-Sys.time()
  elapse_J_a[m] <- difftime(end2,start2,units="secs")#total time
  elapse_J_cliq_a[m] <- sum(elapse_J_cliq1[,m])#the time of constructing clique tree
  elapse_J_t[m] <- sum(elapse_J_t1[,m])#the time of J-Method_cd
  elapse_J_pt[m] <- sum(elapse_J_pt1[,m])#the time of parameter estimation
  elapse_clique_infer[m]<-sum(elapse_clique_infer2[,m])#the time of J-Method_cr
  elapse_C1_t[m]<-sum(elapse_C1_t1[,m])#the time of J-Method_fr
  m<-m+1
}



#related running time
s11<-rbind(t(elapse_J_a),t(elapse_J_pt),t(elapse_J_cliq_a),t(elapse_J_t),t(elapse_clique_infer),t(elapse_C1_t))



s44<-round(s11,3)
row.names(s44)<-c("J_a","J_pt","J_cliq","J_clique_d","J_clique_r","J_factor_r")

s55<-t(s44)
s55


c_0ae<-apply(c_0a,2,mean)
c_1ae<-apply(c_1a,2,mean)

c_0be<-apply(c_0b,2,mean)
c_1be<-apply(c_1b,2,mean)

c_02e<-apply(c_04a,2,mean)
c_12e<-apply(c_14a,2,mean)

sample<- Z
type <- rep(c("J_Method_cd","J_Method_cr","J_Method_fr"),each = s)
value <- c(c_0ae,c_0be,c_02e)
df <- data.frame(sample =sample, type = type, value = value)
Method<-factor(type)
p2<-ggplot(data = df,  mapping = aes(x = sample, y = value, geom="smooth", linetype=Method, shape=Method, color = Method)) +
  geom_point(size=6)+
  theme_bw()+
  geom_xspline(aes(colour = Method), linetype =5, size=1.5)+
  scale_size(range=c(5, 6))+
  xlab('Sample Size')+
  ylab('Bias')+
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(family = 'serif'),
        legend.position = 'top',
        legend.direction = 'horizontal')

p2

sample1<- Z
type1 <- rep(c("J_Method_cd","J_Method_cr","J_Method_fr"),each = s)
value1 <- c(c_1ae,c_1be,c_12e)
df1 <- data.frame(sample =sample1, type = type1, value = value1)
Method<-factor(type1)
p4<-ggplot(data = df1,  mapping = aes(x = sample, y = value, geom="smooth", linetype=Method,shape=Method,color =Method)) +
  geom_point(size=6)+
  theme_bw()+
  geom_xspline(aes(colour = Method), linetype =5, size=1.5)+
  scale_size(range=c(5, 6))+
  xlab('Sample Size')+
  ylab('RMSE')+
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(family = 'serif'),
        legend.position = 'top',
        legend.direction = 'horizontal')

p4

