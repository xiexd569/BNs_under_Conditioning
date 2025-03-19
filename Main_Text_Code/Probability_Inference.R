##in the discrete case, do probability reasoning assuming that related local conditional probabilities are known
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
T_cpt <- vector()
C_cpt1 <- vector()
C_cpt2 <- vector()
C_cpt3 <- vector()
C_cpt4 <- vector()

IC_cpt <- vector()
c_a1 <- vector()
c_a2 <- vector()
c_a3 <- vector()
c_a4 <- vector()

c_3 <- vector()
c_a01 <- vector()
c_a02 <- vector()
c_a03 <- vector()
c_a04 <- vector()

c_03 <- vector()
t <- 100


#given the Bayesian network (BN) structure
asianet <- model2network(
  "[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung][xray|either][dysp|bronc:either]"
)
asianet4 <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][either|tub:lung]")
asianet2a <- model2network(
  "[asia][smoke|tub][tub|asia][lung|smoke:tub][bronc|smoke:either][either|tub:lung:smoke]"
)
asianet2b <- model2network(
  "[asia][smoke|asia][tub|asia:smoke][lung|smoke:tub][bronc|either:smoke][either|lung:smoke:tub]"
)
asianet2c <- model2network(
  "[asia][smoke|tub][tub|asia][lung|bronc:smoke:tub][bronc|smoke:tub][either|bronc:lung:tub]"
)
asianet2d <- model2network(
  "[asia][smoke|asia][tub|asia:bronc:lung][lung|asia:bronc:smoke][bronc|asia:smoke][either|bronc:lung:tub]"
)

yn <- c("yes", "no")

#generating discrete distributions compatitive with the BN structure
for (i in 1:t) {
  a <- runif(18, 0.000001, 0.999999)
  cptA <- matrix(c(a[1], 1 - a[1]), ncol = 2, dimnames = list(NULL, yn))
  cptS <- matrix(c(a[2], 1 - a[2]), ncol = 2, dimnames = list(NULL, yn))
  cptT <- matrix(c(a[3], 1 - a[3], a[4], 1 - a[4]),
                 ncol = 2,
                 dimnames = list(tub = yn, asia = yn))
  cptL <- matrix(c(a[5], 1 - a[5], a[6], 1 - a[6]),
                 ncol = 2,
                 dimnames = list(lung = yn, smoke = yn))
  cptB <- matrix(c(a[7], 1 - a[7], a[8], 1 - a[8]),
                 ncol = 2,
                 dimnames = list(bronc = yn, smoke = yn))
  cptE <- c(a[9], 1 - a[9], a[10], 1 - a[10], a[11], 1 - a[11], a[12], 1 -
              a[12])
  dim(cptE) <- c(2, 2, 2)
  dimnames(cptE) <- list(either = yn,
                         lung = yn,
                         tub = yn)
  cptX <- matrix(c(a[13], 1 - a[13], a[14], 1 - a[14]),
                 ncol = 2,
                 dimnames = list(xray = yn, either = yn))
  cptD <- c(a[15], 1 - a[15], a[16], 1 - a[16], a[17], 1 - a[17], a[18], 1 -
              a[18])
  dim(cptD) <- c(2, 2, 2)
  dimnames(cptD) <- list(dysp = yn,
                         either = yn,
                         bronc = yn)
  asianet.fit <- custom.fit(
    asianet,
    dist = list(
      asia = cptA,
      smoke = cptS,
      tub = cptT,
      lung = cptL,
      bronc = cptB,
      either = cptE,
      xray = cptX,
      dysp = cptD
    )
  )
  
  asia.jtree <- compile(as.grain(asianet.fit))
  asia.ev1 <- setFinding(asia.jtree,
                         nodes = c("xray", "dysp"),
                         states = c("yes", "yes"))
  
  cp2 <- as.data.frame(querygrain(asia.jtree, nodes = "asia", type = 'marginal'))
  cp2a <- as.data.frame(querygrain(asia.ev1, nodes = "asia", type = 'marginal'))
  cp2b <- as.data.frame(querygrain(asia.ev1, nodes = "asia", type = 'marginal'))
  
  cp2c <- as.data.frame(querygrain(asia.ev1, nodes = "asia", type = 'marginal'))
  cp2d <- as.data.frame(querygrain(asia.ev1, nodes = "asia", type = 'marginal'))
  cp2e <- as.data.frame(querygrain(asia.ev1, nodes = "asia", type = 'marginal'))
  
  
  
  cp3 <- as.data.frame(querygrain(asia.jtree, nodes = "smoke", type = 'marginal'))
  cp3a <- as.data.frame(querygrain(asia.ev1, nodes = c("smoke", "tub"), type =
                                     'conditional'))
  cp3b <- as.data.frame(querygrain(asia.ev1, nodes = "smoke", type = 'marginal'))
  
  cp3c <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("smoke", "asia"),
    type = 'conditional'
  ))
  cp3d <- as.data.frame(querygrain(asia.ev1, nodes = c("smoke", "tub"), type =
                                     'conditional'))
  cp3e <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("smoke", "asia"),
    type = 'conditional'
  ))
  
  cp4 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("tub", "asia"),
    type = 'conditional'
  ))
  cp4_1 <- as.data.frame(querygrain(asia.jtree, nodes = "tub", type = 'marginal'))
  cp4a <- as.data.frame(querygrain(asia.ev1, nodes = c("tub", "asia"), type =
                                     'conditional'))
  cp4b <- as.data.frame(querygrain(asia.ev1, nodes = c("tub", "asia"), type =
                                     'conditional'))
  
  cp4c <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("tub", "asia", "smoke"),
    type = 'conditional'
  ))
  cp4d <- as.data.frame(querygrain(asia.ev1, nodes = c("tub", "asia"), type =
                                     'conditional'))
  cp4e <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("tub", "asia", "bronc", "lung"),
    type = 'conditional'
  ))
  
  cp5 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("lung", "smoke"),
    type = 'conditional'
  ))
  cp5a <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("lung", "smoke", "tub"),
    type = 'conditional'
  ))
  cp5b <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("lung", "smoke"),
    type = 'conditional'
  ))
  
  cp5c <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("lung", "smoke", "tub"),
    type = 'conditional'
  ))
  cp5d <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("lung", "bronc", "smoke", "tub"),
    type = 'conditional'
  ))
  cp5e <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("lung", "asia", "bronc", "smoke"),
    type = 'conditional'
  ))
  
  
  cp6 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("either", "tub", "lung"),
    type = 'conditional'
  ))
  cp6a <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("either", "tub", "lung", "smoke"),
    type = 'conditional'
  ))
  cp6b <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("either", "tub", "lung"),
    type = 'conditional'
  ))
  
  cp6c <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("either", "lung", "smoke", "tub"),
    type = 'conditional'
  ))
  cp6d <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("either", "bronc", "lung", "tub"),
    type = 'conditional'
  ))
  cp6e <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("either", "bronc", "lung", "tub"),
    type = 'conditional'
  ))
  
  cp7 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("bronc", "smoke"),
    type = 'conditional'
  ))
  cp7a <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("bronc", "smoke", "either"),
    type = 'conditional'
  ))
  cp7b <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("bronc", "smoke"),
    type = 'conditional'
  ))
  
  cp7c <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("bronc", "smoke", "either"),
    type = 'conditional'
  ))
  cp7d <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("bronc", "smoke", "tub"),
    type = 'conditional'
  ))
  cp7e <- as.data.frame(querygrain(
    asia.ev1,
    nodes = c("bronc", "smoke", "asia"),
    type = 'conditional'
  ))
  
  cp8 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("xray", "either"),
    type = 'conditional'
  ))
  cp9 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("dysp", "either", "bronc"),
    type = 'conditional'
  ))
  cp12 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("xray", "dysp"),
    type = 'joint'
  ))
  cp10 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("xray", "dysp", "smoke"),
    type = 'conditional'
  ))
  cp11 <- as.data.frame(querygrain(
    asia.jtree,
    nodes = c("dysp", "smoke"),
    type = 'conditional'
  ))
  
  #true values calculated according to the factorization formula of the joint distributions
  T_cpt[i] <- cp2$asia[2] * cp3$smoke[2] * cp4$no[1] * cp5$no[1] * cp6$yes.yes[1] *
    cp7$no[2] * cp8$yes[1] * cp9$yes.no[1] / cp12$yes[1]#(cp10$yes.no[1]*cp11$no[1])
  #values calculated according to the factorization formula of the conditional distributions (by our theory results)
  C_cpt1[i] <- cp2a$asia[2] * cp3a$yes[2] * cp4a$no[1] * cp5a$no.yes[1] *
    cp6a$yes.yes.no[1] * cp7a$no.yes[2]
  
  C_cpt2[i] <- cp2c$asia[2] * cp3c$no[2] * cp4c$no.no[1] * cp5c$no.yes[1] *
    cp6c$yes.no.yes[1] * cp7c$no.yes[2]
  C_cpt3[i] <- cp2d$asia[2] * cp3d$yes[2] * cp4d$no[1] * cp5d$no.no.yes[1] *
    cp6d$no.yes.yes[1] * cp7d$no.yes[2]
  C_cpt4[i] <- cp2e$asia[2] * cp3e$no[2] * cp4e$no.no.yes[1] * cp5e$no.no.no[1] *
    cp6e$no.yes.yes[1] * cp7e$no.no[2]
  #values calculated according to the factorization formula of the conditional distributions, ignring the changes of the structure of the conditional distributions
  IC_cpt[i] <- cp2b$asia[2] * cp3b$smoke[2] * cp4b$no[1] * cp5b$no[1] *
    cp6b$yes.yes[1] * cp7b$no[2]
}

for (j in 1:t) {
  c_a1[j] <- (C_cpt1[j] - T_cpt[j]) / T_cpt[j]#relative bias of C_Method_1
  c_a2[j] <- (C_cpt2[j] - T_cpt[j]) / T_cpt[j]#relative bias of C_Method_2
  c_a3[j] <- (C_cpt3[j] - T_cpt[j]) / T_cpt[j]#relative bias of C_Method_3
  c_a4[j] <- (C_cpt4[j] - T_cpt[j]) / T_cpt[j]#relative bias of C_Method_4
  
  c_3[j] <- (IC_cpt[j] - T_cpt[j]) / T_cpt[j]#relative bias of IC_Method
  
  
  c_a01[j] <- C_cpt1[j] - T_cpt[j]#bias of C_Method_1
  c_a02[j] <- C_cpt2[j] - T_cpt[j]#bias of C_Method_2
  c_a03[j] <- C_cpt3[j] - T_cpt[j]#bias of C_Method_3
  c_a04[j] <- C_cpt4[j] - T_cpt[j]#bias of C_Method_4
  
  c_03[j] <- IC_cpt[j] - T_cpt[j]#bias of IC_Method
}

#save the true values and the values obtained by C_Method and IC_Method
s1 <- cbind(T_cpt, C_cpt1, C_cpt2, C_cpt3, IC_cpt)

#save results of relative biases and biases of C_Method and IC_Method
s2 <- cbind(c_a01, c_a02, c_a03, c_a04, c_03, c_a1, c_a2, c_a3, c_a4, c_3)




#the boxplot for the biases of C_Method and IC_Method
sample <- 1:t
type <- rep(c(
  'C-Method_1',
  'C-Method_2',
  'C-Method_3',
  'C-Method_4',
  'IC-Method'
),
each = t)
value <- c(c_a01, c_a02, c_a03, c_a04, c_03)
df <- data.frame(sample = sample,
                 type = type,
                 value = value)
Method <- factor(type)
p1 <- ggplot(
  data = df,
  mapping = aes(
    x = sample,
    y = value,
    geom = "smooth",
    shape = Method,
    color = Method
  )
) +
  xlab('Indexes of 100 Different Discrete Probability Distributions') +
  ylab('Bias') +
  geom_point(size = 6, alpha = 1) +
  theme_bw() +
  scale_shape_manual(values = c(3, 4, 5, 6, 20)) +
  facet_zoom(ylim = c(-0.005, 0.005)) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(family = 'serif'),
    legend.position = 'top',
    legend.direction = 'horizontal'
  )
p1




#the boxplot for the relative biases of C_Method and IC_Method
sample <- 1:t
type <- rep(c(
  'C-Method_1',
  'C-Method_2',
  'C-Method_3',
  'C-Method_4',
  'IC-Method'
),
each = t)
value <- c(c_a1, c_a2, c_a3, c_a4, c_3)
df <- data.frame(sample = sample,
                 type = type,
                 value = value)
Method <- factor(type)

p2 <- ggplot(df, mapping = aes(
  x = sample,
  y = value,
  shape = Method,
  color = Method
)) +
  geom_point(size = 6, alpha = 1) +
  theme_bw() +
  xlab("Indexes of 100 Different Discrete Probability Distributions") +
  ylab("Relative Bias") +
  scale_shape_manual(values = c(3, 4, 5, 6, 20)) +
  facet_zoom(ylim = c(-1, 1)) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(family = 'serif'),
    legend.position = 'top',
    legend.direction = 'horizontal'
  )
p2


