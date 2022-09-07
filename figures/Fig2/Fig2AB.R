rm(list = ls(all=TRUE))
library(coloc)
library(dplyr)
library(tidyr)
library(moloc)
library(devtools)
library(hyprcoloc)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(ggplot2)
library(mvtnorm)
library(biomaRt)

lep.ra <- read.table("chr10.leprosy.txt", header = T)
ld.malawi <- read.table("chr10.malawi.ld", header = T)
ld.mali <- read.table("chr10.mali.ld", header = T)

lep.ra$ld.malawi <- ld.malawi$R2[match(lep.ra$rsid, ld.malawi$SNP_B)]
lep.ra$ld.mali <- ld.mali$R2[match(lep.ra$rsid, ld.mali$SNP_B)]

lep.ra$r2 <- lep.ra$ld.malawi
lep.ra$r2.2 <- lep.ra$ld.mali

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"Spectral")

lep.ra$bin_r2 <- 1
lep.ra$bin_r2[which(lep.ra$r2>0.1 & lep.ra$r2 <= 0.3)] <- 2
lep.ra$bin_r2[which(lep.ra$r2>0.3 & lep.ra$r2 <= 0.5)] <- 3
lep.ra$bin_r2[which(lep.ra$r2>0.5 & lep.ra$r2 <= 0.7)] <- 4
lep.ra$bin_r2[which(lep.ra$r2>0.7)] <- 5

lep.ra$bin_r2.2 <- 1
lep.ra$bin_r2.2[which(lep.ra$r2.2>0.1 & lep.ra$r2.2 <= 0.3)] <- 2
lep.ra$bin_r2.2[which(lep.ra$r2.2>0.3 & lep.ra$r2.2 <= 0.5)] <- 3
lep.ra$bin_r2.2[which(lep.ra$r2.2>0.5 & lep.ra$r2.2 <= 0.7)] <- 4
lep.ra$bin_r2.2[which(lep.ra$r2.2>0.7)] <- 5

lep.ra$genotyped <- 0
lep.ra$genotyped[grep("JHU", lep.ra$rsid)] <- 1

lep.ra$annotate <- 0
lep.ra$annotate[c(which(lep.ra$rsid=="rs2015583:G"))] <- 1
cols3 <- brewer.pal(8,"Paired")

actr1a_plot <- ggplot(lep.ra, aes(x=position, y=-log10(p))) + 
  xlim(104050000, 104450000) +
  geom_point(data=subset(lep.ra, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(lep.ra, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(lep.ra, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(lep.ra, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(lep.ra, bin_r2==5), color=cols2[1], size=3) + 
  geom_point(data=subset(lep.ra, bin_r2.2==1), color=cols[8], size=1.5) + 
  geom_point(data=subset(lep.ra, bin_r2.2==2), color=cols2[5], size=1.5) + 
  geom_point(data=subset(lep.ra, bin_r2.2==3), color=cols2[3], size=1.5) + 
  geom_point(data=subset(lep.ra, bin_r2.2==4), color=cols2[2], size=1.5) + 
  geom_point(data=subset(lep.ra, bin_r2.2==5), color=cols2[1], size=1.5) + 
  geom_point(data=subset(lep.ra, genotyped==1), color="black", size=1, shape=3) + 
  ylab("-log P-value") + 
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(NA, 8.5) +
  geom_text_repel( data=subset(lep.ra, annotate==1), aes(label=rsid), size=6, col = c("black"), nudge_y = 0.5) +
  annotate("text", x = 104070000, y = 8.0, label = c("Malawi")) +
  annotate("text", x = 104130000, y = 8.0, label = c("Mali")) +
  annotate("text", x = 104100000, y = 8.0, label = quote(r^2), hjust = 0.5) +
  annotate("point", x = 104070000, y = 7.6, size = 3, colour = cols2[1]) +
  annotate("point", x = 104130000, y = 7.6, size = 1.5, colour = cols2[1]) +
  annotate("text", x = 104100000, y = 7.6, label = c("0.7-1.0"), hjust = 0.5) +
  annotate("point", x = 104070000, y = 7.2, size = 3, colour = cols2[2]) +
  annotate("point", x = 104130000, y = 7.2, size = 1.5, colour = cols2[2]) +
  annotate("text", x = 104100000, y = 7.2, label = c("0.5-0.7"), hjust = 0.5) +
  annotate("point", x = 104070000, y = 6.8, size = 3, colour = cols2[3]) +
  annotate("point", x = 104130000, y = 6.8, size = 1.5, colour = cols2[3]) +
  annotate("text", x = 104100000, y = 6.8, label = c("0.3-0.5"), hjust = 0.5) +
  annotate("point", x = 104070000, y = 6.4, size = 3, colour = cols2[5]) +
  annotate("point", x = 104130000, y = 6.4, size = 1.5, colour = cols2[5]) +
  annotate("text", x = 104100000, y = 6.4, label = c("0.1-0.3"), hjust = 0.5) +
  annotate("point", x = 104070000, y = 6.1, size = 3, colour = cols[8]) +
  annotate("point", x = 104130000, y = 6.1, size = 1.5, colour = cols[8]) +
  annotate("text", x = 104100000, y = 6.1, label = c("<0.1"), hjust = 0.5)


#genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=10
sel.pos=104200000
range=2500000

out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2

genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding" & start_position>104050000)

genes$start <- genes$start_position
genes$end <- genes$end_position

genes$start[which(genes$strand==-1)] <- genes$end_position[which(genes$strand==-1)]
genes$end[which(genes$strand==-1)] <- genes$start_position[which(genes$strand==-1)]

genes$order <- rep(seq(1:3),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(104050000, 104450000) +
  ylim(c(1,5)) +
  geom_segment(data = genes,
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes, aes(x=mid, label=external_gene_name), size=4, col = c("black"),
                   nudge_y =-0.1, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


recomb <- read.table("genetic_map_chr10_b37_combined.txt", header = T)
recomb_actr1a <- subset(recomb, position>104050000 & position<104450000)

recomb_rate <- ggplot(recomb_actr1a, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 10") +
  scale_x_continuous(breaks=c(104100000,104200000,104300000,104400000),
                     labels=c("140.1Mb", "140.2Mb", "140.3Mb", "140.4Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


#order is MB, PB, malawi, mali

beta.pbmb <- c(-0.623069,-0.671437)
se.pbmb <- c(0.172837,0.135604)
cor <- 0.00533819
beta.overall <- -0.676509
se.overall <- 0.117611
beta.country <- c(-0.660892, -0.699955)
se.country <- c(0.151807, 0.186011)

beta <- beta.pbmb
se <- se.pbmb
cor <- cor
n <- length(beta)

#1st faceted country plot


label <- c("Malawi", "Mali", "Meta")
mean  <- c(beta.country, beta.overall)
lower <- c(beta.country, beta.overall)-(1.96*c(se.country, se.overall))
upper <- c(beta.country, beta.overall)+(1.96*c(se.country, se.overall))
snp <- c(rep("rs2015583:G", 3))

df <- data.frame(label, mean, lower, upper, snp)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])
df$snp <- factor(df$snp)


cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
country.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols2[c(3)], cols[c(3,3)])) + scale_colour_manual(values = c(cols2[c(3)], cols[c(3,3)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  ylim(-1.5,0.5)+
  facet_wrap( ~ snp, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
country.fp

library("meta")

or.fem = metagen(beta.country, se.country, sm = "OR", byvar = c("Malawi", "Mali"))

#het p=0.8708 cochran Q - for dq

#then FP for PB, MB, Leprosy per se
label <- c("MB", "PB", "Overall")
mean  <- c(beta.pbmb, beta.overall)
lower <- c(beta.pbmb, beta.overall)-(1.96*c(se.pbmb, se.overall))
upper <- c(beta.pbmb, beta.overall)+(1.96*c(se.pbmb, se.overall))
snp <- rep("rs2015583:G", 3)

df <- data.frame(label, mean, lower, upper, snp)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])
df$snp <- factor(df$snp)


cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
pbmb.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(3,2,1)])) + scale_colour_manual(values = c(cols[c(3,2,1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ snp, ncol=1) +
  ylim(-1.5,0.5)+
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
pbmb.fp

beta <- beta.pbmb
se <- se.pbmb
cor <- cor
n <- length(beta)

#set rho to 1, i.e. fixed effects
rho <- 1

#likely effect size under alternative model
sigma <- 0.2

#study variability matrix
study.mat <- matrix(cor,nrow=n,ncol=n)
diag(study.mat) <- se^2

#prior matrix
prior.mat <- matrix(rho,nrow=n,ncol=n)
diag(prior.mat) <- 1
prior.mat <- sigma^2 * prior.mat

#log10 BF alternative/null
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

#Dummy array of all possible prior matrices
prior.mat <- array(NA,c(nrow(subset),n,n))

#matrix of study variability
study.mat <- matrix(0,nrow=n,ncol=n)
diag(study.mat) <- se^2

#all possible values of prior.mat
for(i in 1:(nrow(subset))){
  prior.mat[i,,] <- related * (sigma * subset[i,]) * rep(sigma * subset[i,], each = n)
}

BF <- rep(1,nrow(subset))

#Calculates the Bayes Factor for each prior matrix
for(i in 1:nrow(subset)){
  BF[i] <- (dmvnorm(beta,sigma=(study.mat+prior.mat[i,,]),log=TRUE) - dmvnorm(beta,sigma=(study.mat),log=TRUE))/log(10)
}

prior.prob <- rep(1,nrow(subset))

#equal weight on all models
prior.prob[1] <- 1/(length(prior.prob))

prior.prob[2:length(prior.prob)] <- (1-prior.prob[1])/(length(prior.prob)-1)

#convert BF to probabilities
prob <- BF + log10(prior.prob/prior.prob[1])
prob <- 10^(prob-(max(prob)+log10(sum(10^(prob-max(prob))))))

prob.pbmb <- prob

df <- data.frame(
  prob = c(prob.pbmb[c(1,3,4,2)]),
  model = c("Null","MB","PB","Both"),
  label = c(round(prob.pbmb,2)),
  snp = rep("rs2015583:G", 4)
)

df$model <- factor(df$model, levels = c("Null","MB","PB","Both"))
df$snp <- factor(df$snp)


prob_plot <- ggplot(df, aes(x = model, y = prob))+
  ylim(NA,1)+
  geom_col(aes(fill = model), width = 0.7)+
  scale_fill_manual(values = cols[c(8,1,2,3)]) +
  scale_colour_manual(values = cols[c(8,1,2,3)]) +
  theme_bw() + ylab("probability") +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title.x = element_blank(), 
        axis.title=element_text(size=14)) +
  facet_wrap( ~ snp, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))

prob_plot

#Fig2A
p1 <- (actr1a_plot/genes.plot/recomb_rate) + plot_layout(heights= c(3, 1, 1))
ggsave(
  "Fig2A.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 

#Fig2B
p1 <- (country.fp/pbmb.fp/prob_plot)

ggsave(
  "Fig2B.jpg",
  width = 3.5,
  height = 7,
  dpi = 300
) 

