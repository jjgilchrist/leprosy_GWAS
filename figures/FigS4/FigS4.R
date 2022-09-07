rm(list = ls(all=TRUE))
library(lumi)
library(ggplot2)
library(RColorBrewer)
library(car)
library(biomaRt)
library(dplyr)


#read in normalised microarray data describing gene expression in whole blood from leprosy patients. Whole blood was left unstimulated or stimulated with M. leprae sonicate.
#data available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100853
exp <- read.table("norm_expn_matrix.txt", header = T, row.names = 1)
#read in sample info
id <- read.table("id_stim.txt", header = T)

#read in gene labels
genes <- read.table("gene_labels.txt", header = T)
genes$ID <- as.character(genes$ID)
genes$gene <- as.character(genes$gene)

#read in non-normalised microarray data describing gene expression in whole blood from leprosy patients. Whole blood was left unstimulated or stimulated with M. leprae sonicate.
#data available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100853
nonorm_total <- read.table("GSE100853_non-normalized.txt",sep="\t", header = T, row.names = 1)

det.p <- nonorm_total[,grep("Detection.Pval", colnames(nonorm_total))]
nonorm <- nonorm_total[,grep("PS", colnames(nonorm_total))]

#retain probes with at least 3 samples with det p-value<0.05
p.no<- c()
for (i in c(1:dim(det.p)[1])){
  p.no[i] <- sum(det.p[i,]<0.05)
}

p3 <- data.frame(cbind(row.names(det.p), p.no))

p3$V1 <- as.character(p3$V1)
p3$p.no <- as.numeric(as.character(p3$p.no))

exp.qc <- exp[p3$V1[which(p3$p.no>2)],]
nonorm.qc <- nonorm[p3$V1[which(p3$p.no>2)],]


nonorm.qc$total <- rowSums(nonorm.qc)
nonorm.qc <- nonorm.qc[rev(order(nonorm.qc$total)),]
nonorm.qc$gene <- genes$gene[match(rownames(nonorm.qc),genes$ID)]

nonorm.qc <- nonorm.qc[!duplicated(nonorm.qc$gene),]
nonorm.qc<- na.omit(nonorm.qc)
rownames(nonorm.qc) <- nonorm.qc$gene


nonorm.qc$total <- NULL
nonorm.qc$gene <- NULL

#add patient labels
colnames(nonorm.qc) <- id$patient

#split into stimulated and non-stimulated samples
no.stim.nonorm.qc <- nonorm.qc[,which(id$stim==0)]
stim.nonorm.qc <- nonorm.qc[,which(id$stim==1)]

#robust spline normalise data
no.stim.norm.qc <- rsn(as.matrix(no.stim.nonorm.qc))
stim.norm.qc <- rsn(as.matrix(stim.nonorm.qc))

#log2 transform
no.stim.norm.log2 <- log2(no.stim.norm.qc)
stim.norm.log2 <- log2(stim.norm.qc)

#calculate gene expression PCs in each condition

pca.no.stim <- prcomp(t(no.stim.norm.log2),center=TRUE,scale.=TRUE)
pc50.no.stim <- pca.no.stim$x[,c(1:50)]

pca.stim <- prcomp(t(stim.norm.log2),center=TRUE,scale.=TRUE)
pc50.stim <- pca.stim$x[,c(1:50)]

pc50.no.stim <- data.frame(pc50.no.stim)
pc50.stim <- data.frame(pc50.stim)

#add normalised expression values of genes of interest to PC matrix - unstimulated samples
pc50.no.stim$ACTR1A <- t(no.stim.norm.log2["ACTR1A",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["ACTR1A",])))]
pc50.no.stim$TMEM180 <- t(no.stim.norm.log2["TMEM180",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["TMEM180",])))]
pc50.no.stim$FBXL15 <- t(no.stim.norm.log2["FBXL15",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["FBXL15",])))]
pc50.no.stim$NFKB2 <- t(no.stim.norm.log2["NFKB2",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["NFKB2",])))]
pc50.no.stim$PSD <- t(no.stim.norm.log2["PSD",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["PSD",])))]
pc50.no.stim$CUEDC2 <- t(no.stim.norm.log2["CUEDC2",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["CUEDC2",])))]
pc50.no.stim$TRIM8 <- t(no.stim.norm.log2["TRIM8",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["TRIM8",])))]
pc50.no.stim$GBF1 <- t(no.stim.norm.log2["GBF1",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["GBF1",])))]
pc50.no.stim$ELOVL3 <- t(no.stim.norm.log2["ELOVL3",])[match(rownames(pc50.no.stim),colnames(t(no.stim.norm.log2["ELOVL3",])))]

#add normalised expression values of genes of interest to PC matrix - stimulated samples
pc50.stim$ACTR1A <- t(stim.norm.log2["ACTR1A",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["ACTR1A",])))]
pc50.stim$TMEM180 <- t(stim.norm.log2["TMEM180",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["TMEM180",])))]
pc50.stim$FBXL15 <- t(stim.norm.log2["FBXL15",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["FBXL15",])))]
pc50.stim$NFKB2 <- t(stim.norm.log2["NFKB2",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["NFKB2",])))]
pc50.stim$PSD <- t(stim.norm.log2["PSD",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["PSD",])))]
pc50.stim$CUEDC2 <- t(stim.norm.log2["CUEDC2",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["CUEDC2",])))]
pc50.stim$TRIM8 <- t(stim.norm.log2["TRIM8",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["TRIM8",])))]
pc50.stim$GBF1 <- t(stim.norm.log2["GBF1",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["GBF1",])))]
pc50.stim$ELOVL3 <- t(stim.norm.log2["ELOVL3",])[match(rownames(pc50.stim),colnames(t(stim.norm.log2["ELOVL3",])))]

#read in rs2015583 genotyping data
geno <- read.table("rs2015583.txt", header = T)
geno$ID <- as.character(geno$ID)

#add genotyping data to each PC/expression matrix
pc50.no.stim$rs2015583 <- geno$rs2015583[match(rownames(pc50.no.stim),geno$ID)]
pc50.stim$rs2015583 <- geno$rs2015583[match(rownames(pc50.stim),geno$ID)]


#plot out forest plots of rs2015583 association with gene expression in stimulated and unstimulated samples
label <- rep(c("ACTR1A", "TMEM180", "FBXL15", "NFKB2", "PSD", "CUEDC2", "TRIM8", "ELOVL3", "GBF1"),2)
mean  <- c(summary(lm(ACTR1A~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(TMEM180~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(FBXL15~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(NFKB2~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(PSD~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(CUEDC2~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(TRIM8~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(ELOVL3~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(GBF1~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,1],
           summary(lm(ACTR1A~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(TMEM180~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(FBXL15~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(NFKB2~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(PSD~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(CUEDC2~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(TRIM8~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(ELOVL3~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1],
           summary(lm(GBF1~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,1])
se  <- c(summary(lm(ACTR1A~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(TMEM180~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(FBXL15~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(NFKB2~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(PSD~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(CUEDC2~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(TRIM8~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(ELOVL3~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(GBF1~rs2015583+as.matrix(pc50.no.stim[,c(1:7)]), data = pc50.no.stim))$coef[2,2],
         summary(lm(ACTR1A~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(TMEM180~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(FBXL15~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(NFKB2~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(PSD~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(CUEDC2~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(TRIM8~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(ELOVL3~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2],
         summary(lm(GBF1~rs2015583+as.matrix(pc50.stim[,c(1:8)]), data = pc50.stim))$coef[2,2])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("Unstimulated",9), rep("Stimulated",9))

df <- data.frame(label, mean, lower, upper, facet)
df$facet <- factor(df$facet, levels = c("Unstimulated", "Stimulated"))

cols <- brewer.pal(8, "Set2")
lep.wb.eqtl.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + 
  aes(fill = factor(label), col=factor(label)) + scale_fill_manual(values = c(rep(cols[8],9))) + scale_colour_manual(values = c(rep(cols[8],9))) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=10),
                     axis.title=element_text(size=10)) +
  facet_wrap( ~ facet, ncol=2, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
lep.wb.eqtl.fp

ggsave("FigS4.jpg", width = 7, height = 7)




