rm(list = ls(all=TRUE))

library("mvtnorm")
library("ggplot2")
library("RColorBrewer")

#order is MB, PB, malawi, mali
#hla-dqb1
beta.pbmb.dq <- c(0.498093,0.803104)
se.pbmb.dq <- c(0.283809,0.199246)
cor.dq <- 0.0206884
beta.overall.dq <- 0.740509
se.overall.dq <- 0.185808
beta.country.dq <- c(0.722436, 0.781431)
se.country.dq <- c(0.223096, 0.335707)


#hla-b
beta.pbmb.b <- c(-1.07313, -1.56597)
se.pbmb.b <- c(0.572997, 0.550447)
cor.b <- 0.0248532
beta.overall.b <- -1.41162	
se.overall.b <- 0.41148
beta.country.b <- c(-1.60026, -1.28288)	
se.country.b <- c(0.646067, 0.533734)

beta <- beta.pbmb.dq
se <- se.pbmb.dq
cor <- cor.dq

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

prob.pbmb.dq <- prob




#assign at top
beta <- beta.pbmb.b
se <- se.pbmb.b
cor <- cor.b

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

prob.pbmb.b <- prob

#faceted per country plot
label <- rep(c("Malawi", "Mali", "Meta"),2)
mean  <- c(beta.country.dq, beta.overall.dq, beta.country.b, beta.overall.b)
lower <- c(beta.country.dq, beta.overall.dq, beta.country.b, beta.overall.b)-(1.96*c(se.country.dq, se.overall.dq, se.country.b, se.overall.b))
upper <- c(beta.country.dq, beta.overall.dq, beta.country.b, beta.overall.b)+(1.96*c(se.country.dq, se.overall.dq, se.country.b, se.overall.b))
hla <- c(rep("HLA-DBQ1*04:02", 3), rep("HLA-B*49:01", 3))

df <- data.frame(label, mean, lower, upper, hla)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])
df$hla <- factor(df$hla, levels=rev(df$hla)[c(4,1)])


cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
country.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols2[c(3)], cols[c(3,3)])) + scale_colour_manual(values = c(cols2[c(3)], cols[c(3,3)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ hla, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
country.fp

library("meta")

or.fem = metagen(beta.country.dq, se.country.dq, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.8836 cochran Q - for dq

or.fem = metagen(beta.country.b, se.country.b, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.7049 cochran Q - for dq


#forest plots for PB, MB, Leprosy per se
label <- rep(c("MB", "PB", "Overall"),2)
mean  <- c(beta.pbmb.dq, beta.overall.dq, beta.pbmb.b, beta.overall.b)
lower <- c(beta.pbmb.dq, beta.overall.dq, beta.pbmb.b, beta.overall.b)-(1.96*c(se.pbmb.dq, se.overall.dq, se.pbmb.b, se.overall.b))
upper <- c(beta.pbmb.dq, beta.overall.dq, beta.pbmb.b, beta.overall.b)+(1.96*c(se.pbmb.dq, se.overall.dq, se.pbmb.b, se.overall.b))
hla <- c(rep("HLA-DBQ1*04:02", 3), rep("HLA-B*49:01", 3))

df <- data.frame(label, mean, lower, upper, hla)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])
df$hla <- factor(df$hla, levels=rev(df$hla)[c(4,1)])


cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
pbmb.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(3,2,1)])) + scale_colour_manual(values = c(cols[c(3,2,1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ hla, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
pbmb.fp


df <- data.frame(
  prob = c(prob.pbmb.dq[c(1,3,4,2)], prob.pbmb.b[c(1,3,4,2)]),
  model = rep(c("Null","MB","PB","Both"),2),
  label = c(round(prob.pbmb.dq,2), round(prob.pbmb.b,2)),
  hla = c(rep("DBQ1*04:02", 4), rep("B*49:01", 4))
)

df$model <- factor(df$model, levels = c("Null","MB","PB","Both"))
df$hla <- factor(df$hla, levels = c("DBQ1*04:02","B*49:01"))


#posterior probabilities plot per country for HLA
prob_plot <- ggplot(df, aes(x = model, y = prob))+
  ylim(NA,1)+
  geom_col(aes(fill = model), width = 0.7)+
  scale_fill_manual(values = cols[c(8,1,2,3)]) +
  scale_colour_manual(values = cols[c(8,1,2,3)]) +
  theme_bw() + ylab("probability") +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title.x = element_blank(), 
        axis.title=element_text(size=14)) +
    facet_wrap( ~ hla, ncol=2) +
    theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
      size = 12, color = "white", face = "bold.italic"))
  
prob_plot

  
p1 <- (country.fp|pbmb.fp|prob_plot) + plot_layout(widths = c(1, 1, 1.5))

ggsave(
  "Fig4BCD.jpg",
  width = 9,
  height = 3.5,
  dpi = 300
) 
