rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(biomaRt)
library(stringr)
library(coloc)
library(mvtnorm)

#order is MB, PB, malawi, mali
#lacc1.1
beta.pbmb.l1 <- c(0.360398,0.285898)
se.pbmb.l1 <- c(0.150027,0.119981)
cor.l1 <- 0.00533748
beta.overall.l1 <- 0.304906
se.overall.l1 <- 0.105869
beta.country.l1 <- c(0.240405, 0.406681)
se.country.l1 <- c(0.135321, 0.169981)
beta.china.l1 <- log(1.97)
se.china.l1 <- (log(2.62)-log(1.97))/1.96

#slc29a2
beta.pbmb.slc <- c(0.316495, 0.218808)
se.pbmb.slc <- c(0.177449, 0.127115)
cor.slc <- 0.0062168
beta.overall.slc <- 0.243745	
se.overall.slc <- 0.114872
beta.country.slc <- c(0.134301, 0.483083)	
se.country.slc <- c(0.138671, 0.205067)
beta.china.slc <- log(1.14)
se.china.slc <- (log(1.19)-log(1.14))/1.96

#assign at top
beta <- beta.pbmb.l1
se <- se.pbmb.l1
cor <- cor.l1

n <- length(beta)

#rho = 1, i.e. fixed effects
rho <- 1

#expected alternative effect size
sigma <- 0.2

#study variability matrix
study.mat <- matrix(cor,nrow=n,ncol=n)
diag(study.mat) <- se^2

#prior effects matrix
prior.mat <- matrix(rho,nrow=n,ncol=n)
diag(prior.mat) <- 1
prior.mat <- sigma^2 * prior.mat

#log10 BF of alt model/null
BF <- (dmvnorm(beta,sigma=(study.mat+prior.mat),log=TRUE) - dmvnorm(beta,sigma=(study.mat),log=TRUE))/log(10)

integer.base.b <-function(x, b=2){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base=2))+1)
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

#table of all possible models of correlation
subset <- integer.base.b(1:(2^(n)))[,-1]
subset <- subset[nrow(subset):1,]

related <- array(rho,c(n,n))
diag(related) <- 1

prior.mat <- array(NA,c(nrow(subset),n,n))

study.mat <- matrix(0,nrow=n,ncol=n)
diag(study.mat) <- se^2
for(i in 1:(nrow(subset))){
  prior.mat[i,,] <- related * (sigma * subset[i,]) * rep(sigma * subset[i,], each = n)
}

BF <- rep(1,nrow(subset))

#Calculates the Bayes Factor for each prior matrix
for(i in 1:nrow(subset)){
  BF[i] <- (dmvnorm(beta,sigma=(study.mat+prior.mat[i,,]),log=TRUE) - dmvnorm(beta,sigma=(study.mat),log=TRUE))/log(10)
}

prior.prob <- rep(1,nrow(subset))

prior.prob[1] <- 1/(length(prior.prob))

#each model has equal weight
prior.prob[2:length(prior.prob)] <- (1-prior.prob[1])/(length(prior.prob)-1)

#converts BF to posterior prob
prob <- BF + log10(prior.prob/prior.prob[1])
prob <- 10^(prob-(max(prob)+log10(sum(10^(prob-max(prob))))))

prob.pbmb.l1 <- prob




#assign at top
beta <- beta.pbmb.slc
se <- se.pbmb.slc
cor <- cor.slc

n <- length(beta)

#rho = 1, i.e. fixed effects
rho <- 1

#expected alternative effect size
sigma <- 0.2

#study variability matrix
study.mat <- matrix(cor,nrow=n,ncol=n)
diag(study.mat) <- se^2

#prior effects matrix
prior.mat <- matrix(rho,nrow=n,ncol=n)
diag(prior.mat) <- 1
prior.mat <- sigma^2 * prior.mat

#log10 BF of alt model/null
BF <- (dmvnorm(beta,sigma=(study.mat+prior.mat),log=TRUE) - dmvnorm(beta,sigma=(study.mat),log=TRUE))/log(10)

integer.base.b <-function(x, b=2){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base=2))+1)
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

#table of all possible models of correlation
subset <- integer.base.b(1:(2^(n)))[,-1]
subset <- subset[nrow(subset):1,]

related <- array(rho,c(n,n))
diag(related) <- 1

prior.mat <- array(NA,c(nrow(subset),n,n))

study.mat <- matrix(0,nrow=n,ncol=n)
diag(study.mat) <- se^2
for(i in 1:(nrow(subset))){
  prior.mat[i,,] <- related * (sigma * subset[i,]) * rep(sigma * subset[i,], each = n)
}

BF <- rep(1,nrow(subset))

#Calculates the Bayes Factor for each prior matrix
for(i in 1:nrow(subset)){
  BF[i] <- (dmvnorm(beta,sigma=(study.mat+prior.mat[i,,]),log=TRUE) - dmvnorm(beta,sigma=(study.mat),log=TRUE))/log(10)
}

prior.prob <- rep(1,nrow(subset))

prior.prob[1] <- 1/(length(prior.prob))

#each model has equal weight
prior.prob[2:length(prior.prob)] <- (1-prior.prob[1])/(length(prior.prob)-1)

#converts BF to posterior prob
prob <- BF + log10(prior.prob/prior.prob[1])
prob <- 10^(prob-(max(prob)+log10(sum(10^(prob-max(prob))))))

prob.pbmb.slc <- prob

#faceted country forest plots

label <- rep(c("Malawi", "Mali", "Meta: Africa", "China"),2)
mean  <- c(beta.country.l1, beta.overall.l1, beta.china.l1, beta.country.slc, beta.overall.slc, beta.china.slc)
lower <- c(beta.country.l1, beta.overall.l1, beta.china.l1, beta.country.slc, beta.overall.slc, beta.china.slc)-(1.96*c(se.country.l1, se.overall.l1, se.china.l1, se.country.slc, se.overall.slc, se.china.slc))
upper <- c(beta.country.l1, beta.overall.l1, beta.china.l1, beta.country.slc, beta.overall.slc, beta.china.slc)+(1.96*c(se.country.l1, se.overall.l1, se.china.l1, se.country.slc, se.overall.slc, se.china.slc))
locus <- c(rep("LACC1:rs3764147", 4), rep("SLC29A3:rs780668", 4))

df <- data.frame(label, mean, lower, upper, locus)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:4)])


cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
country.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[8], cols2[c(3)], cols[c(3,3)])) + scale_colour_manual(values = c(cols[8], cols2[c(3)], cols[c(3,3)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) + scale_y_continuous(breaks = c(-1.0,-0.5,0,0.5,1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0"), limits = c(-1,1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ locus, ncol=3) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
country.fp

#use meta to ask if heterogeneity for this...

library("meta")

or.fem = metagen(beta.country.l1, se.country.l1, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.4441 cochran Q - for dq

or.fem = metagen(beta.country.l2, se.country.l2, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.8383 cochran Q - for dq

or.fem = metagen(beta.country.slc, se.country.slc, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.1589 cochran Q - for dq

#then FP for PB, MB, Leprosy per se
label <- rep(c("MB", "PB", "Overall"),2)
mean  <- c(beta.pbmb.l1, beta.overall.l1, beta.pbmb.slc, beta.overall.slc)
lower <- c(beta.pbmb.l1, beta.overall.l1, beta.pbmb.slc, beta.overall.slc)-(1.96*c(se.pbmb.l1, se.overall.l1, se.pbmb.slc, se.overall.slc))
upper <- c(beta.pbmb.l1, beta.overall.l1, beta.pbmb.slc, beta.overall.slc)+(1.96*c(se.pbmb.l1, se.overall.l1, se.pbmb.slc, se.overall.slc))
locus <- c(rep("LACC1:rs3764147", 3), rep("SLC29A3:rs780668", 3))

df <- data.frame(label, mean, lower, upper, locus)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])



cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
pbmb.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(3,2,1)])) + scale_colour_manual(values = c(cols[c(3,2,1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) + scale_y_continuous(breaks = c(-1.0,-0.5,0,0.5,1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0"), limits = c(-1,1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ locus, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
pbmb.fp


#faceted barplots for LACC1 and SLC28A3 associations with leprosy and subtypes

df <- data.frame(
  prob = c(prob.pbmb.l1[c(1,3,4,2)], prob.pbmb.slc[c(1,3,4,2)]),
  model = rep(c("Null","MB","PB","Both"),2),
  label = c(round(prob.pbmb.l1,2), round(prob.pbmb.slc,2)),
  locus <- c(rep("LACC1", 4), rep("SLC29A3", 4))
)

df$model <- factor(df$model, levels = c("Null","MB","PB","Both"))
colnames(df)[4] <- "locus"

prob_plot <- ggplot(df, aes(x = model, y = prob))+
  ylim(NA,1)+
  geom_col(aes(fill = model), width = 0.7)+
  scale_fill_manual(values = cols[c(8,1,2,3)]) +
  scale_colour_manual(values = cols[c(8,1,2,3)]) +
  theme_bw() + ylab("probability") +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title.x = element_blank(), 
        axis.title=element_text(size=14)) +
  facet_wrap( ~ locus, ncol=3) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))

prob_plot


p1 <- (country.fp)/(pbmb.fp|prob_plot)


ggsave(
  "Fig5.jpg",
  width = 8,
  height = 7,
  dpi = 300
) 

