##in the discrete case, do probability reasoning after related local conditional probabilities are estimated by simulation data

#set random seed
# set.seed(19260817)
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


# t<-100
# t1<-100
# t<-10
# t1<-10
t<-1
t1<-10000
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
#Z<-c(100,500,1000,2500, 5000,7500,10000)
Z<-c(10000)
yn<-c("yes","no")

#given the Bayesian network (BN) structures
asianet1 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung][xray|either][dysp|bronc:either]")
asianet4 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung]")
asianet2a <- model2network("[asia][smoke|tub][tub|asia][lung|smoke:tub][bronc|either:smoke][either|tub:lung:smoke]")


#generating discrete distributions compatitive with the BN structure 1, from which data are generated
for (i in 1:t) {
  # t<-1
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
  asia.ev1<- setFinding( asia.jtree, nodes = c("xray", "dysp"), states = c("yes", "yes"))
  
  a1<-querygrain(asia.ev1,nodes = "asia", type ='marginal')
  cptA<-matrix(c(a1$asia[1], 1-a1$asia[1]), ncol = 2, dimnames = list(NULL, yn))
  a2<-querygrain(asia.ev1,nodes =c("tub","asia"), type ='conditional')
  cptT <- matrix(c(a2[1], 1-a2[1], a2[3], 1-a2[3]), ncol = 2, dimnames = list(tub = yn, 
                                                                              asia = yn))
  a3<-querygrain(asia.ev1,nodes = c("either","tub","lung","smoke"), type ='conditional')#e|t,l,s
  cptE <- c(a3[1], 1-a3[1], a3[3], 1-a3[3], a3[5], 1-a3[5], a3[7], 1-a3[7], 
            a3[9],1-a3[9],a3[11],1-a3[11],a3[13],1-a3[13],a3[15],1-a3[15])
  dim(cptE) <- c(2, 2, 2,2)
  dimnames(cptE) <- list(either = yn, lung = yn, tub = yn, smoke=yn)
  
  a4<-querygrain(asia.ev1,nodes = c("smoke","tub"), type ='conditional')
  cptS<-matrix(c(a4[1], 1-a4[1], a4[3], 1-a4[3]), ncol = 2, dimnames = list(smoke = yn, 
                                                                            tub = yn))
  
  
  a5<-querygrain(asia.ev1,nodes = c("lung","smoke","tub"), type ='conditional')#l|s,t
  cptL <- c(a5[1], 1-a5[1], a5[3], 1-a5[3], a5[5], 1-a5[5], a5[7], 1-a5[7])
  dim(cptL) <- c(2, 2, 2)
  dimnames(cptL) <- list(lung = yn, smoke = yn, tub = yn)
  
  
  a6<-querygrain(asia.ev1,nodes = c("bronc","smoke","either"), type ='conditional')#b|s,e
  cptB <- c(a6[1], 1-a6[1], a6[3], 1-a6[3],a6[5],1-a6[5],a6[7],1-a6[7])
  dim(cptB) <- c(2, 2,2)
  dimnames(cptB) <- list(bronc = yn, smoke = yn,either=yn)
  
  # a7<-querygrain(asia.ev1,nodes = c("dysp","either","bronc"), type ='conditional')#d|e,b
  # cptD <- c(a7[1], 1-a7[1], a7[3], 1-a7[3], a7[5], 1-a7[5], a7[7], 1-a7[7])
  # dim(cptD) <- c(2, 2, 2)
  # dimnames(cptD) <- list(dysp = yn, either = yn, bronc = yn)
  
  
  ##C1-方法的分布
  asianet.fit_1[[i]]<-custom.fit(asianet2a, dist = list(asia = cptA, smoke = cptS, tub = cptT, 
                                                        lung = cptL, bronc = cptB, either = cptE))
  
}

set.seed(65432)
m<-1
#Z<-100
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {     
    data1<-rbn(asianet.fit[[i]],n)
    for (k in 1:t1) {
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
      # b22<-as.data.frame(querygrain(asia.jtree1,nodes = c("asia","tub","smoke","lung","either","bronc","xray","dysp"), type ='joint'))
      # b23<-as.matrix(b22)
      b10<-as.data.frame(querygrain(asia.jtree1,nodes = c("xray","dysp"), type ='joint'))
      end<-Sys.time()
      elapse_clique_infer1[k] <- difftime(end,start,units="secs")
      
      start <- Sys.time()
      #values calculated according to the factorization formula of the local CPTs obtained by of the data from the joint distributions
      b1<-as.vector(bn.bayes1$asia$prob)
      b2<-as.data.frame(bn.bayes1$bronc$prob)
      b3<-as.data.frame(bn.bayes1$dysp$prob)
      b4<-as.data.frame(bn.bayes1$either$prob)
      b5<-as.data.frame(bn.bayes1$lung$prob)
      b6<-as.data.frame(bn.bayes1$smoke$prob)
      b7<-as.data.frame(bn.bayes1$tub$prob)
      b8<-as.data.frame(bn.bayes1$xray$prob)
      J_cpt[k]<- b1[2]*b2$Freq[4]*b3$Freq[3]*b4$Freq[1]*b5$Freq[3]*b6$Freq[2]*b7$Freq[3]*b8$Freq[1]/b10$yes[1]
      end <- Sys.time()
      elapse_J[k] <- difftime(end,start,units="secs")
      
      # c_02[k]<-unlist(J_cpt[k])-T_cpt[i]
      # c_12[k]<-(unlist(J_cpt[k])-T_cpt[i])^2
      
    }
    # # aa<-Sys.time()
    # c_0a[i,m]<-mean(c_02)#Bias
    # c_1a[i,m]<-sqrt(mean(c_12))#RMSE
    elapse_J_cliq1[i,m]<-sum(elapse_J_cliq)
    elapse_J_t1[i,m] <- sum(elapse_J)
    elapse_J_pt1[i,m] <- sum(elapse_J_p)
    elapse_clique_infer2[i,m]<-sum(elapse_clique_infer1)
  }
  end2<-Sys.time()
  elapse_J_a[m] <- difftime(end2,start2,units="secs")#总时间
  elapse_J_cliq_a[m] <- sum(elapse_J_cliq1[,m])#团树
  elapse_J_t[m] <- sum(elapse_J_t1[,m])#推断
  elapse_J_pt[m] <- sum(elapse_J_pt1[,m])#估计
  elapse_clique_infer[m]<-sum(elapse_clique_infer2[,m])#用团树估计归一化常数
  m<-m+1
}


#C_Method-order 1
# set.seed(19260817+k)
set.seed(65432)
m<-1
for(n in Z){
  start2<-Sys.time()
  for (i in 1:t) {   
    data2<-rbn(asianet.fit_1[[i]],n)
    for (k in 1:t1) {
      start1 <- Sys.time()
      #Bayesian estimators
      bn.bayes2 <- bn.fit(asianet2a, data = data2, method = "bayes", iss = 10)
      end1 <- Sys.time()
      elapse_C1_p[k] <- difftime(end1,start1,units="secs")
      
      start <- Sys.time()
      c1<-as.vector(bn.bayes2$asia$prob)#a
      c2<-as.data.frame(bn.bayes2$bronc$prob)#b|e,s
      c3<-as.data.frame(bn.bayes2$either$prob)#e|l,t
      c4<-as.data.frame(bn.bayes2$lung$prob)#l|s,t
      c5<-as.data.frame(bn.bayes2$smoke$prob)#s|t
      c6<-as.data.frame(bn.bayes2$tub$prob)#t|a
      #c7<-as.data.frame(bn.bayes2$dysp$prob)#d|b,e
      #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      C_cpt[k]<-c1[2]*c2$Freq[6]*c3$Freq[5]*c4$Freq[3]*c5$Freq[2]*c6$Freq[3]
      # 
      # c1<-as.vector(bn.bayes2$asia$prob)
      # c2<-as.data.frame(bn.bayes2$bronc$prob)
      # c3<-as.data.frame(bn.bayes2$either$prob)
      # c4<-as.data.frame(bn.bayes2$lung$prob)
      # c5<-as.data.frame(bn.bayes2$smoke$prob)
      # c6<-as.data.frame(bn.bayes2$tub$prob)
      # c7<-as.data.frame(bn.bayes2$dysp$prob)
      # #values calculated according to the factorization formula of local CPTs obtained by the data from the conditional distributions (by our theory results)
      # C_cpt[k]<-c1[2]*c2$Freq[6]*c3$Freq[5]*c4$Freq[3]*c5$Freq[2]*c6$Freq[3]
      end <- Sys.time()
      elapse_C1[k] <- difftime(end,start,units="secs")
      
      # c_04[k]<-C_cpt[k]-T_cpt[i]
      # c_14[k]<-(C_cpt[k]-T_cpt[i])^2
    }
    
    # c_04a[i,m]<-mean(c_04)#Bias
    # c_14a[i,m]<-sqrt(mean(c_14))#RMSE 
    elapse_C1_t1[i,m] <- sum(elapse_C1)
    elapse_C1_pt1[i,m] <- sum(elapse_C1_p)
  }
  end2<-Sys.time()
  elapse_C1_a[m]<-difftime(end2,start2,units="secs")
  elapse_C1_t[m] <- sum(elapse_C1_t1[,m])#推断
  elapse_C1_pt[m] <- sum(elapse_C1_pt1[,m])#估计
  m<-m+1
}
# elapse_C1_pt


#ROS算法找DAG极小独立图
ancestors_with_self <- function(G, node) {
  # 获取该节点的祖先
  ancestors_set <- bnlearn::ancestors(G, node)
  # 将节点本身加入祖先集合
  return(union(ancestors_set, node))
}

# 输入：DAG G, 拓扑排序 alpha, 和子集 C
# 获取 DAG 的拓扑排序


minimal_fill_edge_set <- function(G, alpha, C) {
  # G<-asianet1
  # C<-c("xray","dysp")
  # 初始化
  topological_order<-alpha
  nodes<-C 
  E_0 <- vector("list", length = 0)  # 初始化空列表
  #更新G
  for(v in C){
    ch<- bnlearn::children(G, v)
    for(u in ch){
      G<-drop.arc(G,v,u)
    }
  }
  
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest<-ancestor_union_all 
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)  # 只考虑子集C的诱导子图
  G_00<-as.bn(G_0)
  
  # 按照 alpha 排序顶点
  k <- length(nodes(G_00))
  vertices <- intersect(rev(topological_order),nodes_of_interest)
  v_structures <- list()
  # 如果顶点数 >= 3，执行填充算法
  if (k >= 3) {
    for (v in vertices[1:(k-2)]) {
      # v<-"dysp"
      parents <- bnlearn::parents(G_00, v)
      if (length(parents) >= 2) {
        comb <- combn(parents, 2)
        for (i in 1:ncol(comb)) {
          x <- comb[1, i]
          y <- comb[2, i]
          if (all(!(x %in% bnlearn::parents(G_00,y)) & !(y %in% bnlearn::parents(G_00,x)))) {
            if(which(topological_order == x)<which(topological_order == y)){
              G_00<- set.arc(G_00, x, y, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1<- paste(x,y, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1    # 添加元素
            }
            else
            {
              G_00<- set.arc(G_00, y, x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1<- paste(y,x, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1    # 添加元素
            }
            
          }
        }
      }
    }
  }  
  if(length(E_0)!=0){
    E_0<-t(unlist(E_0))                 # 最后转换为向量 
  }
  
  for(v in C){
    G_00<-remove.node(G_00,v)
    
  }
  #   # 对每个元素中的单词按字母排序，并重新拼接
  #   sorted_edges <- sapply(E_0, function(edge) {
  #      paste(sort(unlist(strsplit(edge, " "))), collapse = " ")
  #      })  
  # 返回结果
  list(E_alpha_R = E_0, G_alpha_R = G_00)
}

load("D:/Users/xie_x/Downloads/Asia.rda")

# asianet <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung][xray|either][dysp|bronc:either]")
# 
# asianet1<-as.igraph(asianet)

network_string <- modelstring(bn)
asianet<-model2network(network_string)
#asianet1<-as.igraph(asianet)
C<-"xray"
# C<-c("xray","dysp")
# C<-c("xray","dysp","asia")
# C<-c("xray","dysp","asia","smoke")
# C<-c("xray","dysp","asia","smoke","tub")
# C<-c("xray","dysp","asia","smoke","tub","lung")
# C<-c("xray","dysp","asia","smoke","tub","lung","bronc")

alpha1<-c("asia","tub","smoke","lung","either","bronc","xray","dysp")



start <- Sys.time()
for(i in 1:10000){
  z1<-minimal_fill_edge_set(asianet,alpha1,C)
}
end <- Sys.time()
t11<- difftime(end,start,units="secs")




#related running time
s1<-rbind(t(elapse_J_a),t(elapse_J_cliq_a),t(elapse_clique_infer),t(elapse_J_pt),t(elapse_J_t),
          t(elapse_C1_a),t(elapse_C1_pt),t(elapse_C1_t))



s4<-round(s1,3)
row.names(s4)<-c("J_a","J_cliq","J_cliq_infer","J_pt","J_t","C1_a","C1_pt","C1_t")


s5<-rbind(t(elapse_J_cliq_a),t(elapse_clique_infer),t(elapse_J_pt),t(elapse_J_t),
          t(elapse_C1_pt),t(elapse_C1_t))
row.names(s5)<-c("J_cliq","J_cliq_infer","J_pt","J_t","C1_pt","C1_t")
s6<-round(t(s5),3)

#s4<-t(s4)
t1a<-t11
t1b<-round(t1a,3)





