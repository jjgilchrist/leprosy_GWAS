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
library(metafor)
library(stringr)
cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")


lep <- read.table("leprosy.chr10.txt", header = T)
lep$SNP2 <- as.character(paste("10:",as.character(lep$position), sep = ""))
lep <- na.omit(lep)
lep <- lep[!duplicated(lep$SNP2),]
lep.list <- list(beta = lep$beta, varbeta = (lep$se)^2, type = "cc", N = 1131, s = 0.4350133, snp = as.character(lep$SNP2))

ibd <- read.table("EUR_IBD.chr10.txt", header = T)
ibd$SNP2 <- as.character(paste("10:",as.character(ibd$POS), sep = ""))
ibd <- ibd[!duplicated(ibd$SNP2),]
ibd$BETA <- log(ibd$OR)
ibd.list <- list(beta = ibd$BETA, varbeta = (ibd$SE)^2, type = "cc", N = 34652, s = 0.3717534, snp = as.character(ibd$SNP2))

uc <- read.table("EUR_UC.chr10.txt", header = T)
uc$SNP2 <- as.character(paste("10:",as.character(uc$POS), sep = ""))
uc <- uc[!duplicated(uc$SNP2),]
uc$BETA <- log(uc$OR)
uc.list <- list(beta = uc$BETA, varbeta = (uc$SE)^2, type = "cc", N = 20883, s = 0.285, snp = as.character(uc$SNP2))

asthma <- read.table("asthma_child.chr10.txt", header = T)
asthma$SNP2 <- as.character(paste("10:",as.character(asthma$base_pair_location), sep = ""))
asthma <- asthma[!duplicated(asthma$SNP2),]
asthma.list <- list(beta = asthma$beta, varbeta = (asthma$standard_error)^2, type = "cc", N = 314633, s = 0.5, snp = as.character(asthma$SNP2))

wcc <- read.table("wcc.chr10.txt", header = T)
wcc$SNP2 <- as.character(paste("10:",as.character(wcc$pos), sep = ""))
wcc <- na.omit(wcc)
wcc <- wcc[!duplicated(wcc$SNP2),]
wcc.list <- list(beta = wcc$beta, varbeta = (wcc$se)^2, type = "quant", N = 350470, MAF = wcc$minor_AF, snp = as.character(wcc$SNP2))

neut <- read.table("neut.chr10.txt", header = T)
neut$SNP2 <- as.character(paste("10:",as.character(neut$pos), sep = ""))
neut <- na.omit(neut)
neut <- neut[!duplicated(neut$SNP2),]
neut.list <- list(beta = neut$beta, varbeta = (neut$se)^2, type = "quant", N = 349856, MAF = neut$minor_AF, snp = as.character(neut$SNP2))


label <- c("IBD", "WCC", "Neut", "Asthma (Child)","UC")
mean  <- c(-ibd$BETA[which(ibd$SNP=="rs2015583")],
           wcc$beta[which(wcc$rsid=="rs2015583")],
           1.65E-02, 
           -log(1.07082),
           -uc$BETA[which(uc$SNP=="rs2015583")])
se  <- c(ibd$SE[which(ibd$SNP=="rs2015583")],
         wcc$se[which(wcc$rsid=="rs2015583")],
         3.36E-03, 
         0.0115906,
         -uc$SE[which(uc$SNP=="rs2015583")])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c("IMD",
           "Haematological traits", 
           "Haematological traits",
           "IMD","IMD")

df <- data.frame(label, mean, lower, upper, facet)
df$label <- factor(df$label, levels = rev(c("IBD", "UC","Asthma (Child)", "WCC", "Neut")))



df$facet <- factor(df$facet, levels = c("IMD", "Haematological traits"))

#plot forest plot of effects at rs2015583

pleio.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[rep(4,8)])) + scale_colour_manual(values = c(cols[rep(4,8)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("beta/log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(-0.1,0,0.1),
                     labels=c("-0.1", "0", "0.1")) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ facet, ncol=1, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
pleio.fp


#colocalisation of signal at leprosy c.f. IMD/haematological traits

lep$SNP3 <- paste0("chr",lep$SNP2)

lep.co <- lep
lep.co$wcc.p <- wcc$pval[match(lep.co$SNP2,wcc$SNP2)]
lep.co$neut.p <- neut$pval[match(lep.co$SNP2,neut$SNP2)]
lep.co$asthma.p <- asthma$p_value[match(lep.co$SNP2,asthma$SNP2)]
lep.co$ibd.p <- ibd$PVAL[match(lep.co$SNP2,ibd$SNP2)]
lep.co$uc.p <- uc$PVAL[match(lep.co$SNP2,uc$SNP2)]
lep.co <- na.omit(lep.co)

#read in LD matrix for colocalisation plot - need to provide token for LDmatrix
ld.out1 <- LDmatrix(na.omit(lep.co$SNP3[c(1:990)]), pop = "CEU", r2d = "r2", token = token, file = FALSE)
rownames(ld.out1) <- ld.out1$RS_number
ld.out1$RS_number <- NULL
lep.co$rsid <- as.character(lep.co$rsid)
lep.co$rsid1 <- NA

for (i in c(1:990)){
  lep.co$rsid1[i] <- strsplit(lep.co$rsid[i], ":")[[1]][1]
}

lep.co$r2 <- ld.out1$rs2015583[match(lep.co$rsid1, rownames(ld.out1))]
lep.co <- na.omit(lep.co)

lep.co$bin_r2 <- 1
lep.co$bin_r2[which(lep.co$r2>0.1 & lep.co$r2 <= 0.3)] <- 2
lep.co$bin_r2[which(lep.co$r2>0.3 & lep.co$r2 <= 0.5)] <- 3
lep.co$bin_r2[which(lep.co$r2>0.5 & lep.co$r2 <= 0.7)] <- 4
lep.co$bin_r2[which(lep.co$r2>0.7)] <- 5

coplot1 <- data.frame(lep.co$p, lep.co$ibd.p, lep.co$SNP2, lep.co$bin_r2, "IBD")
coplot2 <- data.frame(lep.co$p, lep.co$asthma.p, lep.co$SNP2, lep.co$bin_r2, "Asthma")
coplot3 <- data.frame(lep.co$p, lep.co$wcc.p, lep.co$SNP2, lep.co$bin_r2, "WCC")
coplot4 <- data.frame(lep.co$p, lep.co$neut.p, lep.co$SNP2, lep.co$bin_r2, "Neut")
coplot5 <- data.frame(lep.co$p, lep.co$uc.p, lep.co$SNP2, lep.co$bin_r2, "UC")

colnames(coplot1) <- c("lep.p", "trait.p", "SNP", "r2.bin", "trait")
colnames(coplot2) <- c("lep.p", "trait.p", "SNP", "r2.bin", "trait")
colnames(coplot3) <- c("lep.p", "trait.p", "SNP", "r2.bin", "trait")
colnames(coplot4) <- c("lep.p", "trait.p", "SNP", "r2.bin", "trait")
colnames(coplot5) <- c("lep.p", "trait.p", "SNP", "r2.bin", "trait")
co.plot <- rbind(coplot1, coplot2, coplot3, coplot4, coplot5)

#facet with pp4s on annotation - no need for r2 colouring legend - as will be on same plot
coloc.abf(lep.list, asthma.list)
coloc.abf(lep.list, ibd.list)
coloc.abf(lep.list, wcc.list)
coloc.abf(lep.list, neut.list)
coloc.abf(lep.list, uc.list)

pp4_labels = data.frame(trait = c("IBD", 
                                  "WCC",
                                  "Neut",
                                  "Asthma",
                                  "UC"),
                        label1 = c("PP4==0.91", 
                                   "PP4==0.94", 
                                   "PP4==0.92", 
                                   "PP4==0.98",
                                   "PP4==0.92"))

co.plot$trait <- factor(co.plot$trait, levels = c("IBD", "UC", "Asthma", "WCC", "Neut"))

cols2 <- brewer.pal(11,"Spectral")
lep.coplot <- ggplot(co.plot, aes(y=-log10(trait.p), x=-log10(lep.p))) + 
  geom_point(data=subset(co.plot, r2.bin==1), color=cols[8], size=3) + 
  geom_point(data=subset(co.plot, r2.bin==2), color=cols2[5], size=3) + 
  geom_point(data=subset(co.plot, r2.bin==3), color=cols2[3], size=3) + 
  geom_point(data=subset(co.plot, r2.bin==4), color=cols2[2], size=3) + 
  geom_point(data=subset(co.plot, r2.bin==5), color=cols2[1], size=3) + 
  xlab("GWAS Leprosy (-logP)") + 
  ylab("GWAS Trait (-logP)") + 
  ylim(NA,9) +
  facet_wrap( ~ trait, ncol=5, scales = "free_x") + geom_text(x=0, y=9, aes(label=label1), data=pp4_labels, parse=TRUE, inherit.aes=F, size=3, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"), plot.title = element_text(size = 20, face = "bold"))

#Evidence for replication in China and India

#China
label <- c("Overall")
mean  <- c(-0.0234580388513035)
lower <- c(-0.0234580388513035)-(1.96*c(0.043230395779606))
upper <- c(-0.0234580388513035)+(1.96*c(0.043230395779606))
locus <- rep("China", 1)

df <- data.frame(label, mean, lower, upper, locus)

cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
pbmb.fp.china <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(3,2,1)])) + scale_colour_manual(values = c(cols[c(3,2,1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) + scale_y_continuous(breaks = c(-1.0,0,1), labels = c("-1.0", "0", "1.0"), limits = c(-1.5,1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ locus, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
pbmb.fp.china

#India
label <- c("MB", "PB", "Overall")
mean  <- c(-0.626, 0.142, -0.044)
lower <- c(-0.626, 0.142, -0.044)-(1.96*c(0.306, 0.278, 0.212))
upper <- c(-0.626, 0.142, -0.044)+(1.96*c(0.306, 0.278, 0.212))
locus <- rep("India", 3)

df <- data.frame(label, mean, lower, upper, locus)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])

cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
pbmb.fp.india <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(3,2,1)])) + scale_colour_manual(values = c(cols[c(3,2,1)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) + scale_y_continuous(breaks = c(-1.0,0,1), labels = c("-1.0", "0", "1.0"), limits = c(-1.5,1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap( ~ locus, ncol=1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
pbmb.fp.india

p1 <- ((pbmb.fp.china/pbmb.fp.india)|(pleio.fp)|(lep.coplot)) + plot_layout(widths = c(1, 1, 2))


ggsave(
  "Fig2CD.jpg",
  width = 10.5,
  height = 4,
  dpi = 300
) 




