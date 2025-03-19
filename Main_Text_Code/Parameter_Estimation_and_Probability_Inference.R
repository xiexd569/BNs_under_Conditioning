##in the discrete case, do probability reasoning after related local conditional probabilities are estimated by simulation data
##including the runtime of J_Method_fr, IC_Method and C_Method
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
J_cpt<-vector()
IC_cpt<-vector()
C_cpt<-vector()


t<-100
t1<-100
s<-7

c_02<-vector()
c_04<-vector()
c_06<-vector()
c_0a<-matrix(rep(0,s*t),nrow=t)
c_0b<-matrix(rep(0,s*t),nrow=t)

c_04a<-matrix(rep(0,s*t),nrow=t)
c_04b<-matrix(rep(0,s*t),nrow=t)
c_04c<-matrix(rep(0,s*t),nrow=t)
c_04d<-matrix(rep(0,s*t),nrow=t)

c_12<-vector()
c_14<-vector()
c_16<-vector()

c_1a<-matrix(rep(0,s*t),nrow=t)
c_1b<-matrix(rep(0,s*t),nrow=t)

c_14a<-matrix(rep(0,s*t),nrow=t)
c_14b<-matrix(rep(0,s*t),nrow=t)
c_14c<-matrix(rep(0,s*t),nrow=t)
c_14d<-matrix(rep(0,s*t),nrow=t)

c_0ae<-vector()
c_0be<-vector()
c_02e<-vector()
c_04e<-vector()
c_06e<-vector()
c_08e<-vector()

c_1ae<-vector()
c_1be<-vector()
c_12e<-vector()
c_14e<-vector()
c_16e<-vector()
c_18e<-vector()


elapse_J<-vector()
elapse_IC<-vector()
elapse_C1<-vector()
elapse_C2<-vector()
elapse_C3<-vector()
elapse_C4<-vector()

elapse_J_p<-vector()
elapse_IC_p<-vector()
elapse_C1_p<-vector()
elapse_C2_p<-vector()
elapse_C3_p<-vector()
elapse_C4_p<-vector()

elapse_J_t<-vector()
elapse_IC_t<-vector()
elapse_C1_t<-vector()
elapse_C2_t<-vector()
elapse_C3_t<-vector()
elapse_C4_t<-vector()

elapse_J_pt<-vector()
elapse_IC_pt<-vector()
elapse_C1_pt<-vector()
elapse_C2_pt<-vector()
elapse_C3_pt<-vector()
elapse_C4_pt<-vector()

elapse_J_a<-vector()
elapse_IC_a<-vector()
elapse_C1_a<-vector()
elapse_C2_a<-vector()
elapse_C3_a<-vector()
elapse_C4_a<-vector()

elapse_J_cliq<-vector()
elapse_J_cliq_a<-vector()

elapse_J_t1<-matrix(rep(0,s*t),nrow=t)
elapse_J_pt1<-matrix(rep(0,s*t),nrow=t)
elapse_J_cliq1<-matrix(rep(0,s*t),nrow=t)
elapse_IC_t1<-matrix(rep(0,s*t),nrow=t)
elapse_IC_pt1<-matrix(rep(0,s*t),nrow=t)
elapse_C1_t1<-matrix(rep(0,s*t),nrow=t)
elapse_C1_pt1<-matrix(rep(0,s*t),nrow=t)

elapse_C2_t1<-matrix(rep(0,s*t),nrow=t)
elapse_C2_pt1<-matrix(rep(0,s*t),nrow=t)
elapse_C3_t1<-matrix(rep(0,s*t),nrow=t)
elapse_C3_pt1<-matrix(rep(0,s*t),nrow=t)
elapse_C4_t1<-matrix(rep(0,s*t),nrow=t)
elapse_C4_pt1<-matrix(rep(0,s*t),nrow=t)

elapse_clique_infer1<-vector()
elapse_clique_infer2<-matrix(rep(0,s*t),nrow=t)
elapse_clique_infer<-vector()

asianet.fit<-list()
asianet.fit_IC<-list()
asianet.fit_1<-list()
asianet.fit_2<-list()
asianet.fit_3<-list()
asianet.fit_4<-list()

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
  asia.ev1<- setFinding(asia.jtree, nodes = c("xray", "dysp"), states = c("yes", "yes")) 
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
  cptB <- matrix(c(a6[1], 1-a6[1], a6[3], 1-a6[3]), ncol = 2, dimnames = list(bronc = yn, 
                                                                              smoke = yn))
  
  ##the distribution of IC-Method
  asianet.fit_IC[[i]]<-custom.fit(asianet4, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  
  
  
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia"), type ='conditional')
  cptT <- matrix(c(a2[1], 1-a2[1], a2[3], 1-a2[3]), ncol = 2, dimnames = list(tub = yn, 
                                                                              asia = yn))
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung","smoke"), type ='conditional')
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7],
            a3[9], 1-a3[9], a3[11], 1-a3[11], a3[13], 1-a3[13], a3[15], 1-a3[15])
  dim(cptE) <- c(2, 2, 2,2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn, smoke=yn)
  
  a4<-querygrain(asia.ev1,nodes = c("smoke","tub"), type ='conditional')
  cptS<-matrix(c(a4[1], 1-a4[1], a4[3], 1-a4[3]), ncol = 2, dimnames = list(smoke = yn, 
                                                                            tub = yn))
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke","tub"), type ='conditional')
  cptL <- c(a5[1], 1-a5[1], a5[3], 1-a5[3], a5[5], 1-a5[5], a5[7], 1-a5[7])
  dim(cptL) <- c(2, 2, 2)
  dimnames(cptL) <- list(lung = yn, smoke = yn, tub = yn)
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke","either"), type ='conditional') 
  cptB <- c(a6[1], 1-a6[1], a6[3], 1-a6[3], a6[5], 1-a6[5], a6[7], 1-a6[7])
  dim(cptB) <- c(2, 2, 2)
  dimnames(cptB) <- list(bronc = yn, smoke = yn, either = yn)
  
  
  
  
  ##the distribution of C-Method_1
  asianet.fit_1[[i]]<-custom.fit(asianet2a, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia","smoke"), type ='conditional')
  cptT <- c(a2[1], 1-a2[1], a2[3], 1-a2[3], a2[5], 1-a2[5], a2[7], 1-a2[7])
  dim(cptT) <- c(2, 2, 2)
  dimnames(cptT) <- list(tub = yn, asia = yn, smoke = yn)
  
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung","smoke"), type ='conditional')
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7],
            a3[9], 1-a3[9], a3[11], 1-a3[11], a3[13], 1-a3[13], a3[15], 1-a3[15])
  dim(cptE) <- c(2, 2, 2,2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn, smoke=yn)
  
  a4<-querygrain(asia.ev1,nodes = c("smoke","asia"), type ='conditional')
  cptS<-matrix(c(a4[1], 1-a4[1], a4[3], 1-a4[3]), ncol = 2, dimnames = list(smoke = yn, 
                                                                            asia = yn))
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke","tub"), type ='conditional')
  cptL <- c(a5[1], 1-a5[1], a5[3], 1-a5[3], a5[5], 1-a5[5], a5[7], 1-a5[7])
  dim(cptL) <- c(2, 2, 2)
  dimnames(cptL) <- list(lung = yn, smoke = yn, tub = yn)
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke","either"), type ='conditional') 
  cptB <- c(a6[1], 1-a6[1], a6[3], 1-a6[3], a6[5], 1-a6[5], a6[7], 1-a6[7])
  dim(cptB) <- c(2, 2, 2)
  dimnames(cptB) <- list(bronc = yn, smoke = yn, either = yn)
  
  
  
  ##the distribution of C-Method_2
  asianet.fit_2[[i]]<-custom.fit(asianet2b, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia"), type ='conditional')
  cptT <- c(a2[1], 1-a2[1], a2[3], 1-a2[3])
  dim(cptT) <- c(2, 2)
  dimnames(cptT) <- list(tub = yn, asia = yn)
  
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung","bronc"), type ='conditional')
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7],
            a3[9], 1-a3[9], a3[11], 1-a3[11], a3[13], 1-a3[13], a3[15], 1-a3[15])
  dim(cptE) <- c(2, 2, 2,2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn, bronc=yn)
  
  a4<-querygrain(asia.ev1,nodes = c("smoke","tub"), type ='conditional')
  cptS<-matrix(c(a4[1], 1-a4[1], a4[3], 1-a4[3]), ncol = 2, dimnames = list(smoke = yn, 
                                                                            tub = yn))
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke","tub","bronc"), type ='conditional')
  cptL <- c(a5[1], 1-a5[1], a5[3], 1-a5[3], a5[5], 1-a5[5], a5[7], 1-a5[7],
            a5[9], 1-a5[9], a5[11], 1-a5[11], a5[13], 1-a5[13], a5[15], 1-a5[15])
  dim(cptL) <- c(2, 2, 2,2)
  dimnames(cptL) <- list(lung = yn, smoke = yn, tub = yn, bronc=yn)
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke","tub"), type ='conditional') 
  cptB <- c(a6[1], 1-a6[1], a6[3], 1-a6[3], a6[5], 1-a6[5], a6[7], 1-a6[7])
  dim(cptB) <- c(2, 2, 2)
  dimnames(cptB) <- list(bronc = yn, smoke = yn, tub = yn)
  
  
  ##the Method of C-Method_3
  asianet.fit_3[[i]]<-custom.fit(asianet2c, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia","bronc","lung"), type ='conditional')
  cptT <- c(a2[1], 1-a2[1], a2[3], 1-a2[3],a2[5], 1-a2[5], a2[7], 1-a2[7],
            a2[9], 1-a2[9], a2[11], 1-a2[11],a2[13], 1-a2[13], a2[15], 1-a2[15])
  dim(cptT) <- c(2, 2,2,2)
  dimnames(cptT) <- list(tub = yn, asia = yn,bronc=yn,lung=yn)
  
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung","bronc"), type ='conditional')
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7],
            a3[9], 1-a3[9], a3[11], 1-a3[11], a3[13], 1-a3[13], a3[15], 1-a3[15])
  dim(cptE) <- c(2, 2, 2,2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn, bronc=yn)
  
  a4<-querygrain(asia.ev1,nodes = c("smoke","asia"), type ='conditional')
  cptS<-matrix(c(a4[1], 1-a4[1], a4[3], 1-a4[3]), ncol = 2, dimnames = list(smoke = yn, 
                                                                            asia = yn))
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke","asia","bronc"), type ='conditional')
  cptL <- c(a5[1], 1-a5[1], a5[3], 1-a5[3], a5[5], 1-a5[5], a5[7], 1-a5[7],
            a5[9], 1-a5[9], a5[11], 1-a5[11], a5[13], 1-a5[13], a5[15], 1-a5[15])
  dim(cptL) <- c(2, 2, 2,2)
  dimnames(cptL) <- list(lung = yn, smoke = yn, asia = yn, bronc=yn)
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke","asia"), type ='conditional') 
  cptB <- c(a6[1], 1-a6[1], a6[3], 1-a6[3], a6[5], 1-a6[5], a6[7], 1-a6[7])
  dim(cptB) <- c(2, 2, 2)
  dimnames(cptB) <- list(bronc = yn, smoke = yn, asia = yn)
  
  ##the distribution of C-Method_4
  asianet.fit_4[[i]]<-custom.fit(asianet2d, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  
}

#the true value
for (i in 1:t) {
  asia.jtree <- compile(as.grain(asianet.fit[[i]]))
  asia.ev1<- setFinding(asia.jtree, nodes = c("xray", "dysp"), states = c("yes", "yes")) 
  
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

#J_Method
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
      b10<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","dysp"), type ='joint'))
      end<-Sys.time()
      elapse_clique_infer1[k] <- difftime(end,start,units="secs")
      
      start <- Sys.time()
      b1<-as.vector(bn.bayes1$asia$prob)
      b2<-as.data.frame(bn.bayes1$bronc$prob)
      b3<-as.data.frame(bn.bayes1$dysp$prob)
      b4<-as.data.frame(bn.bayes1$either$prob)
      b5<-as.data.frame(bn.bayes1$lung$prob)
      b6<-as.data.frame(bn.bayes1$smoke$prob)
      b7<-as.data.frame(bn.bayes1$tub$prob)
      b8<-as.data.frame(bn.bayes1$xray$prob) 
      #values calculated according to the factorization formula of the local CPTs obtained by of the data from the joint distributions
      J_cpt[k]<-b1[2]*b2$Freq[4]*b3$Freq[3]*b4$Freq[1]*b5$Freq[3]*b6$Freq[2]*b7$Freq[3]*b8$Freq[1]/b10$yes[1]
      end <- Sys.time()
      elapse_J[k] <- difftime(end,start,units="secs")
      
      c_02[k]<-J_cpt[k]-T_cpt[i]
      c_12[k]<-(J_cpt[k]-T_cpt[i])^2
      
    }
    c_0a[i,m]<-mean(c_02)#Bias
    c_1a[i,m]<-sqrt(mean(c_12))#RMSE
    elapse_J_cliq1[i,m]<-sum(elapse_J_cliq)
    elapse_J_t1[i,m] <- sum(elapse_J)
    elapse_J_pt1[i,m] <- sum(elapse_J_p)
    elapse_clique_infer2[i,m]<-sum(elapse_clique_infer1)
  }
  end2<-Sys.time()
  elapse_J_a[m] <- difftime(end2,start2,units="secs")#total time
  elapse_J_cliq_a[m] <- sum(elapse_J_cliq1[,m])#the time for clique tree
  elapse_J_t[m] <- sum(elapse_J_t1[,m])#the time for inference
  elapse_J_pt[m] <- sum(elapse_J_pt1[,m])#the time for parameter estimation
  elapse_clique_infer[m]<-sum(elapse_clique_infer2[,m])#the time for calculating normalized constant
  m<-m+1
}

#IC_Method
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    
    for (k in 1:t1) {
      data2<-rbn(asianet.fit_IC[[i]],n)
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes4 <- bn.fit(asianet4, data = data2, method = "bayes", iss =10)
      end1 <- Sys.time()
      elapse_IC_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      e1<-as.vector(bn.bayes4$asia$prob)
      e2<-as.data.frame(bn.bayes4$bronc$prob)
      e3<-as.data.frame(bn.bayes4$either$prob)
      e4<-as.data.frame(bn.bayes4$lung$prob)
      e5<-as.data.frame(bn.bayes4$smoke$prob)
      e6<-as.data.frame(bn.bayes4$tub$prob)
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions, ignring the changes of the structure of the conditional distributions
      IC_cpt[k]<-e1[2]*e2$Freq[4]*e3$Freq[1]*e4$Freq[3]*e5$Freq[2]*e6$Freq[3]
      end <- Sys.time()
      elapse_IC[k] <- difftime(end,start,units="secs")
      
      c_06[k]<-IC_cpt[k]-T_cpt[i]
      c_16[k]<-(IC_cpt[k]-T_cpt[i])^2
    }
    c_0b[i,m]<-mean(c_06)#Bias
    c_1b[i,m]<-sqrt(mean(c_16))#RMSE
    elapse_IC_t1[i,m] <- sum(elapse_IC)
    elapse_IC_pt1[i,m] <- sum(elapse_IC_p)
  }
  end2<-Sys.time()
  elapse_IC_a[m] <- difftime(end2,start2,units="secs")#total time
  elapse_IC_t[m] <- sum(elapse_IC_t1[,m])#the time for inference
  elapse_IC_pt[m] <- sum(elapse_IC_pt1[,m])#the time for parameter estimation
  m<-m+1
}



#C_Method-order 1
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    for (k in 1:t1) {
      data2<-rbn(asianet.fit_1[[i]],n)
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes2 <- bn.fit(asianet2a, data = data2, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_C1_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      c1<-as.vector(bn.bayes2$asia$prob)
      c2<-as.data.frame(bn.bayes2$bronc$prob)
      c3<-as.data.frame(bn.bayes2$either$prob)
      c4<-as.data.frame(bn.bayes2$lung$prob)
      c5<-as.data.frame(bn.bayes2$smoke$prob)
      c6<-as.data.frame(bn.bayes2$tub$prob)
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      C_cpt[k]<-c1[2]*c2$Freq[6]*c3$Freq[5]*c4$Freq[3]*c5$Freq[2]*c6$Freq[3]
      end <- Sys.time()
      elapse_C1[k] <- difftime(end,start,units="secs")
      
      c_04[k]<-C_cpt[k]-T_cpt[i]
      c_14[k]<-(C_cpt[k]-T_cpt[i])^2
    }
    
    c_04a[i,m]<-mean(c_04)#Bias
    c_14a[i,m]<-sqrt(mean(c_14))#RMSE 
    elapse_C1_t1[i,m] <- sum(elapse_C1)
    elapse_C1_pt1[i,m] <- sum(elapse_C1_p)
  }
  end2<-Sys.time()
  elapse_C1_a[m]<-difftime(end2,start2,units="secs")#total time
  elapse_C1_t[m] <- sum(elapse_C1_t1[,m])#the time for inference
  elapse_C1_pt[m] <- sum(elapse_C1_pt1[,m])# the time for parameter estimation 
  m<-m+1
}


#C_Method-order 2
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    for (k in 1:t1) {
      data2<-rbn(asianet.fit_2[[i]],n)
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes2 <- bn.fit(asianet2b, data = data2, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_C2_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      c1<-as.vector(bn.bayes2$asia$prob)
      c2<-as.data.frame(bn.bayes2$bronc$prob)
      c3<-as.data.frame(bn.bayes2$either$prob)
      c4<-as.data.frame(bn.bayes2$lung$prob)
      c5<-as.data.frame(bn.bayes2$smoke$prob)
      c6<-as.data.frame(bn.bayes2$tub$prob)
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      C_cpt[k]<-c1[2]*c2$Freq[6]*c3$Freq[5]*c4$Freq[3]*c5$Freq[4]*c6$Freq[7]
      end <- Sys.time()
      elapse_C2[k] <- difftime(end,start,units="secs")
      
      c_04[k]<-C_cpt[k]-T_cpt[i]
      c_14[k]<-(C_cpt[k]-T_cpt[i])^2
    }
    
    c_04b[i,m]<-mean(c_04)#Bias
    c_14b[i,m]<-sqrt(mean(c_14))#RMSE 
    elapse_C2_t1[i,m] <- sum(elapse_C2)
    elapse_C2_pt1[i,m] <- sum(elapse_C2_p)
  }
  end2<-Sys.time()
  elapse_C2_a[m]<-difftime(end2,start2,units="secs")#total time
  elapse_C2_t[m] <- sum(elapse_C2_t1[,m])#the time for inference
  elapse_C2_pt[m] <- sum(elapse_C2_pt1[,m])#the time for parameter estimation
  m<-m+1
}




#C_Method-order 3
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    for (k in 1:t1) {
      data2<-rbn(asianet.fit_3[[i]],n)
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes2 <- bn.fit(asianet2c, data = data2, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_C3_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      c1<-as.vector(bn.bayes2$asia$prob)
      c2<-as.data.frame(bn.bayes2$bronc$prob)
      c3<-as.data.frame(bn.bayes2$either$prob)
      c4<-as.data.frame(bn.bayes2$lung$prob)
      c5<-as.data.frame(bn.bayes2$smoke$prob)
      c6<-as.data.frame(bn.bayes2$tub$prob)
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      C_cpt[k]<-c1[2]*c2$Freq[4]*c3$Freq[3]*c4$Freq[7]*c5$Freq[2]*c6$Freq[3]
      end <- Sys.time()
      elapse_C3[k] <- difftime(end,start,units="secs")
      
      c_04[k]<-C_cpt[k]-T_cpt[i]
      c_14[k]<-(C_cpt[k]-T_cpt[i])^2
    }
    
    c_04c[i,m]<-mean(c_04)#Bias
    c_14c[i,m]<-sqrt(mean(c_14))#RMSE 
    elapse_C3_t1[i,m] <- sum(elapse_C3)
    elapse_C3_pt1[i,m] <- sum(elapse_C3_p)
  }
  end2<-Sys.time()
  elapse_C3_a[m]<-difftime(end2,start2,units="secs")#total time 
  elapse_C3_t[m] <- sum(elapse_C3_t1[,m])#the time for inference
  elapse_C3_pt[m] <- sum(elapse_C3_pt1[,m])#the tine for parameter estimation
  m<-m+1
}


#C_Method-order 4
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {
    for (k in 1:t1) {
      data2<-rbn(asianet.fit_4[[i]],n)
      
      
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes2 <- bn.fit(asianet2d, data = data2, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_C4_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      c1<-as.vector(bn.bayes2$asia$prob)
      c2<-as.data.frame(bn.bayes2$bronc$prob)
      c3<-as.data.frame(bn.bayes2$either$prob)
      c4<-as.data.frame(bn.bayes2$lung$prob)
      c5<-as.data.frame(bn.bayes2$smoke$prob)
      c6<-as.data.frame(bn.bayes2$tub$prob)
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      C_cpt[k]<-c1[2]*c2$Freq[8]*c3$Freq[3]*c4$Freq[15]*c5$Freq[4]*c6$Freq[7]
      end <- Sys.time()
      elapse_C4[k] <- difftime(end,start,units="secs")
      
      c_04[k]<-C_cpt[k]-T_cpt[i]
      c_14[k]<-(C_cpt[k]-T_cpt[i])^2
    }
    
    c_04d[i,m]<-mean(c_04)#Bias
    c_14d[i,m]<-sqrt(mean(c_14))#RMSE 
    elapse_C4_t1[i,m] <- sum(elapse_C4)
    elapse_C4_pt1[i,m] <- sum(elapse_C4_p)
  }
  end2<-Sys.time()
  elapse_C4_a[m]<-difftime(end2,start2,units="secs")#total time
  elapse_C4_t[m] <- sum(elapse_C4_t1[,m])#the time for inference
  elapse_C4_pt[m] <- sum(elapse_C4_pt1[,m])#the time for parameter estimation
  m<-m+1
}

#load data

load("D:/Users/xie_x/Downloads/Asia.rda")

#set network structure and conditioning set
network_string <- modelstring(bn)
asianet<-model2network(network_string)
C<-c("xray","dysp")


##Finding ancestors and themselves
ancestors_with_self <- function(G, node) {
  # Retrieve the ancestors of the given node
  ancestors_set <- bnlearn::ancestors(G, node)
  
  # Include the node itself in the ancestor set
  return(union(ancestors_set, node))
}

#Finding minimal filling edge set by the ROS algorithm
minimal_fill_edge_set <- function(G, alpha, C) {
  
  topological_order <- alpha
  nodes <- C 
  E_0 <- vector("list", length = 0)  # Initialize an empty list
  
  # Update graph by removing outgoing edges from nodes in C
  for (v in C) {
    ch <- bnlearn::children(G, v)
    for (u in ch) {
      G <- drop.arc(G, v, u)
    }
  }
  
  # Compute the set of ancestors for the nodes in C
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest <- ancestor_union_all 
  
  # Induce a subgraph containing only the relevant nodes
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)  
  G_00 <- as.bn(G_0)
  
  # Order vertices according to the given topological order
  k <- length(nodes(G_00))
  vertices <- intersect(rev(topological_order), nodes_of_interest)
  v_structures <- list()
  
  # Apply minimal fill-in algorithm if there are at least 3 vertices
  if (k >= 3) {
    for (v in vertices[1:(k-2)]) {
      parents <- bnlearn::parents(G_00, v)
      if (length(parents) >= 2) {
        comb <- combn(parents, 2)
        for (i in 1:ncol(comb)) {
          x <- comb[1, i]
          y <- comb[2, i]
          
          # Ensure that x and y are not already connected
          if (all(!(x %in% bnlearn::parents(G_00, y)) & !(y %in% bnlearn::parents(G_00, x)))) {
            if (which(topological_order == x) < which(topological_order == y)) {
              G_00 <- set.arc(G_00, x, y, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(x, y, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1  # Add edge to list
            } else {
              G_00 <- set.arc(G_00, y, x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(y, x, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1  # Add edge to list
            }
          }
        }
      }
    }
  }  
  
  # Convert edge list to vector if not empty
  if (length(E_0) != 0) {
    E_0 <- t(unlist(E_0))
  }
  
  # Remove nodes in C from the final graph
  for (v in C) {
    G_00 <- remove.node(G_00, v)
  }
  
  # Return the modified graph and the added edges
  list(E_alpha_R = E_0, G_alpha_R = G_00)
}

#four topological orders
alpha1<-c("asia","tub","smoke","lung","either","bronc","xray","dysp")
alpha2<-c("asia","smoke","tub","lung","either","bronc","xray","dysp")
alpha3<-c("asia","tub","smoke","bronc","lung","either","xray","dysp")
alpha4<-c("asia","smoke","bronc","lung","tub","either","xray","dysp")

start <- Sys.time()
for(i in 1:10000){
  z1<-minimal_fill_edge_set(asianet,alpha1,C)
}
end <- Sys.time()
t11<- difftime(end,start,units="secs")


start <- Sys.time()
for(i in 1:10000){
  z2<-minimal_fill_edge_set(asianet,alpha2,C)
}
end <- Sys.time()
t22<- difftime(end,start,units="secs")


start <- Sys.time()
for(i in 1:10000){
  z3<-minimal_fill_edge_set(asianet,alpha3,C)
}
end <- Sys.time()
t33<- difftime(end,start,units="secs")


start <- Sys.time()
for(i in 1:10000){
  z4<-minimal_fill_edge_set(asianet,alpha4,C)
}
end <- Sys.time()
t44<- difftime(end,start,units="secs")



#related running time
s1<-rbind(t(elapse_J_a),t(elapse_J_cliq_a),t(elapse_clique_infer),t(elapse_J_pt),t(elapse_J_t),t(elapse_IC_a),t(elapse_IC_pt),t(elapse_IC_t),
          t(elapse_C1_a),t(elapse_C1_pt),t(elapse_C1_t), t(elapse_C2_a),t(elapse_C2_pt),t(elapse_C2_t), t(elapse_C3_a),t(elapse_C3_pt),t(elapse_C3_t),
          t(elapse_C4_a), t(elapse_C4_pt),t(elapse_C4_t))



s4<-round(s1,3)
row.names(s4)<-c("J_a","J_cliq","J_cliq_infer","J_pt","J_t","IC_a","IC_pt","IC_t","C1_a","C1_pt","C1_t","C2_a",
                 "C2_pt","C2_t","C3_a","C3_pt","C3_t","C4_a","C4_pt","C4_t")


s5<-rbind(t(elapse_J_cliq_a),t(elapse_clique_infer),t(elapse_J_pt),t(elapse_J_t),t(elapse_IC_pt),t(elapse_IC_t),
          t(elapse_C1_pt),t(elapse_C1_t), t(elapse_C2_pt),t(elapse_C2_t), t(elapse_C3_pt),t(elapse_C3_t),
          t(elapse_C4_pt),t(elapse_C4_t))
row.names(s5)<-c("J_cliq","J_cliq_infer","J_pt","J_t","IC_pt","IC_t","C1_pt","C1_t",
                 "C2_pt","C2_t","C3_pt","C3_t","C4_pt","C4_t")
s6<-round(t(s5),3)

t1a<-c(t11,t22,t33,t44)
t1b<-round(t1a,3)


c_0ae<-apply(c_0a,2,mean)
c_1ae<-apply(c_1a,2,mean)
c_0be<-apply(c_0b,2,mean)
c_1be<-apply(c_1b,2,mean)

c_02e<-apply(c_04a,2,mean)
c_12e<-apply(c_14a,2,mean)
c_04e<-apply(c_04b,2,mean)
c_14e<-apply(c_14b,2,mean)
c_06e<-apply(c_04c,2,mean)
c_16e<-apply(c_14c,2,mean)
c_08e<-apply(c_04d,2,mean)
c_18e<-apply(c_14d,2,mean)


#the boxplots of biases for the C_Method with 4 different orders
f0a<-t(cbind(t(c_0a[,1]),t(c_0a[,2]),t(c_0a[,3]),t(c_0a[,4]),t(c_0a[,5]),t(c_0a[,6]),t(c_0a[,7])))
typea1 <- rep(Z,each = t)
typea1 <-as.factor(typea1)
f0a<-data.frame(f0a)
f0a<-cbind(f0a,typea1)

f1a<-t(cbind(t(c_0b[,1]),t(c_0b[,2]),t(c_0b[,3]),t(c_0b[,4]),t(c_0b[,5]),t(c_0b[,6]),t(c_0b[,7])))
typea1 <- rep(Z,each = t)
typea1 <-as.factor(typea1)
f1a<-data.frame(f1a)
f1a<-cbind(f1a,typea1)



f2a<-t(cbind(t(c_04a[,1]),t(c_04a[,2]),t(c_04a[,3]),t(c_04a[,4]),t(c_04a[,5]),t(c_04a[,6]),t(c_04a[,7])))
typea1 <- rep(Z,each = t)
typea1 <-as.factor(typea1)
f2a<-data.frame(f2a)
f2a<-cbind(f2a,typea1)

f2b<-t(cbind(t(c_04b[,1]),t(c_04b[,2]),t(c_04b[,3]),t(c_04b[,4]),t(c_04b[,5]),t(c_04b[,6]),t(c_04b[,7])))
typea2 <- rep(Z,each = t)
typea2<-as.factor(typea2)
f2b<-data.frame(f2b)
f2b<-cbind(f2b,typea2)

f2c<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
f2c<-t(cbind(t(c_04c[,1]),t(c_04c[,2]),t(c_04c[,3]),t(c_04c[,4]),t(c_04c[,5]),t(c_04c[,6]),t(c_04c[,7])))
typea3 <- rep(Z,each =t)
typea3<-as.factor(typea3)
f2c<-data.frame(f2c)
f2c<-cbind(f2c,typea3)

f2d<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
f2d<-t(cbind(t(c_04d[,1]),t(c_04d[,2]),t(c_04d[,3]),t(c_04d[,4]),t(c_04d[,5]),t(c_04d[,6]),t(c_04d[,7])))
typea4 <- rep(Z,each = t)
typea4<-as.factor(typea4)
f2d<-data.frame(f2d)
f2d<-cbind(f2d,typea4)


g<-cbind(f0a[,1],f1a[,1],f2a[,1],f2b[,1],f2c[,1], f2d)
names(g)<-c("J_Method","IC_Method","C_Method_1","C_Method_2","C_Method_3","C_Method_4")
typeb<-rep(c("J_Method","IC_Method","C_Method_1","C_Method_2","C_Method_3","C_Method_4"),each=s*t)
typeb<-as.factor(typeb)


g1<-as.matrix(g[,1:6])
g1<-as.numeric(g1)
typeb1 <- rep(Z,each = t)
typeb1 <- t(cbind(t(typeb1),t(typeb1),t(typeb1),t(typeb1),t(typeb1),t(typeb1)))
typeb2<-rep(c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4"),each=s*t)
g1<-data.frame(g1)
g1<-cbind(g1,typeb1,typeb2)
names(g1)<-c("Bias","Sample", "Method")
g1$Method<-as.factor(g1$Method)
g1$Sample<-as.factor(g1$Sample)




p1<-ggplot(g1,mapping = aes(x=Sample, y=Bias, fill=Method))+
  theme_bw()+
  xlab("Sample Size")+
  ylab("Bias")+
  stat_boxplot(geom="errorbar",width=0.3,linewidth=1,position=position_dodge(0.6))+
  geom_boxplot(aes(fill=Method),
               position=position_dodge(0.6),
               size=0.5,
               width=0.5,
               #color="blue",
               outlier.color = "black",
               outlier.fill = "white",
               outlier.shape = 1,
               outlier.size = 4,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(family = 'serif'),
        legend.position = 'top',
        legend.direction = 'horizontal')

p1


#the plots of means of biases for the C_Method with 4 different orders
sample<- Z
type <- rep(c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4"),each = s)
value <- c(c_0ae,c_0be,c_02e,c_04e,c_06e,c_08e)
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


#the boxplots of RMSEs for the C_Method with 4 different orders
r0a<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r0a<-t(cbind(t(c_1a[,1]),t(c_1a[,2]),t(c_1a[,3]),t(c_1a[,4]),t(c_1a[,5]),t(c_1a[,6]),t(c_1a[,7])))
r0a<-data.frame(r0a)


r1a<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r1a<-t(cbind(t(c_1b[,1]),t(c_1b[,2]),t(c_1b[,3]),t(c_1b[,4]),t(c_1b[,5]),t(c_1b[,6]),t(c_1b[,7])))
r1a<-data.frame(r1a)

r2a<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r2a<-t(cbind(t(c_14a[,1]),t(c_14a[,2]),t(c_14a[,3]),t(c_14a[,4]),t(c_14a[,5]),t(c_14a[,6]),t(c_14a[,7])))
r2a<-data.frame(r2a)


r2b<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r2b<-t(cbind(t(c_14b[,1]),t(c_14b[,2]),t(c_14b[,3]),t(c_14b[,4]),t(c_14b[,5]),t(c_14b[,6]),t(c_14b[,7])))
r2b<-data.frame(r2b)

r2c<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r2c<-t(cbind(t(c_14c[,1]),t(c_14c[,2]),t(c_14c[,3]),t(c_14c[,4]),t(c_14c[,5]),t(c_14c[,6]),t(c_14c[,7])))
typea3 <- rep(Z,each = t)
typea3<-as.factor(typea1)
r2c<-data.frame(r2c)

r2d<-matrix(rep(0,s*t),nrow=s*t,ncol=1)
r2d<-t(cbind(t(c_14d[,1]),t(c_14d[,2]),t(c_14d[,3]),t(c_14d[,4]),t(c_14d[,5]),t(c_14d[,6]),t(c_14d[,7])))
typea4 <- rep(Z,each = t)
typea4<-as.factor(typea4)
r2d<-data.frame(r2d)
r2d<-cbind(r2d,typea4)


t2<-cbind(r0a,r1a,r2a,r2b,r2c, r2d)
names(t2)<-c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4","Sample Size")
typeb1<-rep(c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4"),each=s*t)
typeb1<-as.factor(typeb1)


t3<-as.matrix(t2[,1:6])
t3<-as.numeric(t3)
typea1 <- rep(Z,each = t)
typea1 <- t(cbind(t(typea1),t(typea1),t(typea1),t(typea1),t(typea1),t(typea1)))
typeb1<-rep(c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4"),each=s*t)
t3<-data.frame(t3)
t3<-cbind(t3,typea1,typeb1)
names(t3)<-c("RMSE","Sample", "Method")
t3$Method<-as.factor(t3$Method)
t3$Sample<-as.factor(t3$Sample)


p3<-ggplot(t3,mapping = aes(x=Sample, y=RMSE, fill=Method))+
  theme_bw()+
  xlab("Sample Size")+
  ylab("RMSE")+
  stat_boxplot(geom="errorbar",width=0.3,size=1,position=position_dodge(0.6))+
  geom_boxplot(aes(fill=Method),
               position=position_dodge(0.6),
               size=0.5,
               width=0.5,
               #color="blue",
               outlier.color = "black",
               outlier.fill = "white",
               outlier.shape = 1,
               outlier.size = 4,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(family = 'serif'),
        legend.position = 'top',
        legend.direction = 'horizontal')
p3



#the plots of means of RMSEs for the C_Method with 4 different orders
sample1<- Z
type1 <- rep(c("J-Method","IC-Method","C-Method_1","C-Method_2","C-Method_3","C-Method_4"),each = s)
value1 <- c(c_1ae,c_1be,c_12e,c_14e,c_16e,c_18e)
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


