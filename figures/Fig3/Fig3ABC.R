rm(list = ls(all=TRUE))
library(RColorBrewer)
library(ggplot2)
library(arrow)
library(stringr)
library(coloc)
library(ggplot2)
library(RColorBrewer)
library(coloc)
library(dplyr)
library(tidyr)
library(patchwork)
library(mvtnorm)
library(ggrepel)


cols <- brewer.pal(8,"Set2")

label <- c("C10orf32", 
           "TMEM180",
           "TRIM8", 
           "TRIM8",
           "TMEM180", 
           "TRIM8", 
           "TMEM180", 
           "ACTR1A",
           "ACTR1A")
mean  <- c(-0.0511235,0.031729665, -0.197655482, 0.042918152, 0.057436937, 0.027146163, 0.042502697, -0.0852074, -0.0658636)
pval <- c(0.000898577, 5.03468E-05, 2.4669E-24, 2.45297E-07, 0.000382668, 0.000138765, 0.000769062, 4.69427E-09, 3.68064E-05)
df <- c(414,414,261,322,322,367,367,293,283)
se  <- abs(mean/abs(qt(pval/2, df = df)))
           
           
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("CD14 - naive",2), rep("CD14 - LPS2",1), rep("CD14 - LPS24",2), rep("CD14 - IFN",2), "CD4+ T cells", "CD8+ T cells")

df <- data.frame(label, mean, lower, upper, facet)
df$facet <- factor(df$facet, levels = c("CD14 - naive", "CD14 - LPS2", "CD14 - LPS24", "CD14 - IFN", "CD4+ T cells", "CD8+ T cells"))
df$cols <- factor(seq(1:9))


primary.cell.eqtl.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = cols, col=cols) + scale_fill_manual(values = c(rep(cols[8],7),cols[4], cols[8])) + scale_colour_manual(values = c(rep(cols[8],7),cols[4], cols[8])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=10),
                     axis.title=element_text(size=10)) +
  facet_wrap( ~ facet, ncol=3) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
primary.cell.eqtl.fp


#read in leprosy association statistics at chr10 locus
lep <- read.table("lep_assn.txt", header = T)
lep <- na.omit(subset(lep, position>(104226983-250000) & position<(104226983+250000)))
lep <- lep[!duplicated(lep$ID2),]

#read in LD matrices for malawi and europe
malawi.ld <- as.matrix(read.table("MALAWI_LD.txt.gz", header = T, row.names = 1))
eur.ld <- as.matrix(read.table("EUR_LD.txt", header = T, row.names = 1))

lep <- lep[which(lep$ID2 %in% colnames(malawi.ld)),]
malawi.ld <- malawi.ld[which(colnames(malawi.ld) %in% lep$ID2),which(colnames(malawi.ld) %in% lep$ID2)]

#construct list for coloc
lep.list <- list(beta = lep$beta, varbeta = (lep$se)^2, type = "cc", N = 1131, s = 0.4350133, snp = as.character(lep$ID2),LD=malawi.ld)

#read in gtex data for chr10 region - retain phenotypes with rs2015583 association p<0.00001
data <- read.table("gtex.Skin_Not_Sun_Exposed_Suprapubic.v8.EUR.allpairs.chr10_region.txt", header = T)
rs2015583.data <- subset(data, variant_id=="chr10_102467226_A_G_b38")
trait <- "Skin_Not_Sun_Exposed_Suprapubic"
trait_n <- 517

df <- data.frame(tissue=character(),
                 gene=character(),
                 beta=numeric(),
                 se=numeric(),
                 pval=numeric(),
                 pp4=numeric(),
                 stringsAsFactors=FALSE)

#perform coloc.susie/coloc.abf for each phenotype vs leprosy associaition
for (i in c(1:dim(rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),])[1])){
  pheno1 <- data.frame(subset(data, phenotype_id==rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]))
  pheno1$ID2 <- NA
  for (j in c(1:dim(pheno1)[1])){
    pheno1$ID2[j]  <- paste0(str_split(pheno1$variant_id[j], "_")[[1]][1],"_",str_split(pheno1$variant_id[j],"_")[[1]][2],"_",str_split(pheno1$variant_id[j],"_")[[1]][4])
  }
  pheno1 <- na.omit(pheno1)
  pheno1 <- pheno1[which(pheno1$ID2 %in% colnames(eur.ld)),]
  eur.ld <- eur.ld[which(colnames(eur.ld) %in% pheno1$ID2),which(colnames(eur.ld) %in% pheno1$ID2)]
  
  
  pheno1.list <- list(beta = pheno1$slope, varbeta = (pheno1$slope_se)^2, type = "quant", N = trait_n, MAF = pheno1$maf, snp = as.character(pheno1$ID2),LD=eur.ld)
  E=runsusie(pheno1.list,coverage=0.01)
  G=runsusie(lep.list, coverage=0.01)
  susie.res=coloc.susie(E,G)
  abf.out=coloc.abf(pheno1.list, lep.list)
  
  df[i,1] <- trait
  df[i,2] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]
  df[i,3] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope[i]
  df[i,4] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope_se[i]
  df[i,5] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$pval_nominal[i]
  df[i,6] <- abf.out$summary[6]
  
}
df.out <- df

#recapitulate for sun-exposed skin
data <- read.table("gtex.Skin_Sun_Exposed_Lower_leg.v8.EUR.allpairs.chr10_region.txt", header = T)
rs2015583.data <- subset(data, variant_id=="chr10_102467226_A_G_b38")
trait <- "Skin_Sun_Exposed_Lower_leg"
trait_n <- 605

rs2015583.data <- subset(data, variant_id=="chr10_102467226_A_G_b38")

df <- data.frame(tissue=character(),
                 gene=character(),
                 beta=numeric(),
                 se=numeric(),
                 pval=numeric(),
                 pp4=numeric(),
                 stringsAsFactors=FALSE)

for (i in c(1:dim(rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),])[1])){
  pheno1 <- data.frame(subset(data, phenotype_id==rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]))
  pheno1$ID2 <- NA
  for (j in c(1:dim(pheno1)[1])){
    pheno1$ID2[j]  <- paste0(str_split(pheno1$variant_id[j], "_")[[1]][1],"_",str_split(pheno1$variant_id[j],"_")[[1]][2],"_",str_split(pheno1$variant_id[j],"_")[[1]][4])
  }
  pheno1 <- na.omit(pheno1)
  pheno1 <- pheno1[which(pheno1$ID2 %in% colnames(eur.ld)),]
  eur.ld <- eur.ld[which(colnames(eur.ld) %in% pheno1$ID2),which(colnames(eur.ld) %in% pheno1$ID2)]
  
  
  pheno1.list <- list(beta = pheno1$slope, varbeta = (pheno1$slope_se)^2, type = "quant", N = trait_n, MAF = pheno1$maf, snp = as.character(pheno1$ID2),LD=eur.ld)
  E=runsusie(pheno1.list,coverage=0.01)
  G=runsusie(lep.list, coverage=0.01)
  susie.res=coloc.susie(E,G)
  abf.out=coloc.abf(pheno1.list, lep.list)
  
  df[i,1] <- trait
  df[i,2] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]
  df[i,3] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope[i]
  df[i,4] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope_se[i]
  df[i,5] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$pval_nominal[i]
  df[i,6] <- abf.out$summary[6]
  
}

df.out <- rbind(df.out, df)
#recapitulate for tibial nerve
data <- read.table("gtex.Nerve_Tibial.v8.EUR.allpairs.chr10_region.txt", header = T)
rs2015583.data <- subset(data, variant_id=="chr10_102467226_A_G_b38")
trait <- "Nerve_Tibial"
trait_n <- 532

rs2015583.data <- subset(data, variant_id=="chr10_102467226_A_G_b38")

df <- data.frame(tissue=character(),
                 gene=character(),
                 beta=numeric(),
                 se=numeric(),
                 pval=numeric(),
                 pp4=numeric(),
                 stringsAsFactors=FALSE)

for (i in c(1:dim(rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),])[1])){
  pheno1 <- data.frame(subset(data, phenotype_id==rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]))
  pheno1$ID2 <- NA
  for (j in c(1:dim(pheno1)[1])){
    pheno1$ID2[j]  <- paste0(str_split(pheno1$variant_id[j], "_")[[1]][1],"_",str_split(pheno1$variant_id[j],"_")[[1]][2],"_",str_split(pheno1$variant_id[j],"_")[[1]][4])
  }
  pheno1 <- na.omit(pheno1)
  pheno1 <- pheno1[which(pheno1$ID2 %in% colnames(eur.ld)),]
  eur.ld <- eur.ld[which(colnames(eur.ld) %in% pheno1$ID2),which(colnames(eur.ld) %in% pheno1$ID2)]
  
  
  pheno1.list <- list(beta = pheno1$slope, varbeta = (pheno1$slope_se)^2, type = "quant", N = trait_n, MAF = pheno1$maf, snp = as.character(pheno1$ID2),LD=eur.ld)
  abf.out=coloc.abf(pheno1.list, lep.list)
  #E=runsusie(pheno1.list,coverage=0.01)
  #G=runsusie(lep.list, coverage=0.01)
  #susie.res=coloc.susie(E,G)
  
  df[i,1] <- trait
  df[i,2] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$phenotype_id[i]
  df[i,3] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope[i]
  df[i,4] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$slope_se[i]
  df[i,5] <- rs2015583.data[which(rs2015583.data$pval_nominal<0.00001),]$pval_nominal[i]
  df[i,6] <- abf.out$summary[6]
  
}

df.out <- rbind(df.out, df)




#gtex FPs then co.plots
label <- c("TMEM180",
           "AS3MT", 
           "TRIM8",
           "AS3MT", 
           "C10orf95-AS1",
           "TMEM180",
           "SFXN2",
           "AS3MT")
mean  <- df.out$beta
pval <- df.out$pval
se  <- df.out$se


lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("Skin",2), rep("Skin (sun exp)",2), rep("Tibial nerve",4))

df <- data.frame(label, mean, lower, upper, facet)
#df$facet <- factor(df$facet, levels = c("CD14 - naive", "CD14 - LPS2", "CD14 - LPS24", "CD14 - IFN", "CD4+ T cells", "CD8+ T cells"))
df$cols <- factor(seq(1:8))


gtex.eqtl.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = cols, col=cols) + scale_fill_manual(values = c(cols[4], rep(cols[8],4),cols[4],rep(cols[8],2))) + scale_colour_manual(values = c(cols[4], rep(cols[8],4),cols[4],rep(cols[8],2))) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=10),
                     axis.title=element_text(size=10)) +
  facet_wrap( ~ facet, ncol=3) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
gtex.eqtl.fp

#construct colocalisation plots
#read in CD4 ACTR1A eQTL stats
cd4.actr1a <- read.table("cd4_actr1a.txt", header= T)
#read in positions to harmonise SNP identifiers
pos = read.table(file = "/Users/jamesgilchrist/Documents/Leprosy/leprosy_manuscript/final_analysis/working/eqtl/est_bb_cd4/rs2015583_region.list_pos", header = F, row.names = 1)
pos$SNP2 <- as.character(paste("10:",as.character(pos$V2), sep = ""))

cd4.actr1a$pos2 <- pos$V2[match(cd4.actr1a$rsid, rownames(pos))]
cd4.actr1a$maf <- pos$maf[match(cd4.actr1a$rsid, rownames(pos))]
cd4.actr1a$SNP2 <- as.character(paste("10:",as.character(cd4.actr1a$pos2), sep = ""))
cd4.actr1a <- cd4.actr1a[which(!is.na(cd4.actr1a$pos2)),]

#add cd4 association stats to leprosy GWAS stats
lep$SNP2 <- as.character(paste("10:",as.character(lep$position), sep = ""))
lep$cd4.p <- cd4.actr1a$pvalue[match(lep$SNP2,cd4.actr1a$SNP2)]
lep$r2 <- cd4.actr1a$r2[match(lep$SNP2,cd4.actr1a$SNP2)]
lep$ID2 <- as.character(lep$ID2)

#read in TMEM180 eqtl stats (tibial nerve)
tib_n <- read.table("gtex.Nerve_Tibial.v8.EUR.allpairs.chr10_region.txt", header = T)
tib_n <- subset(tib_n, phenotype_id=="ENSG00000138111.14")
tib_n$ID2 <- NA
for (j in c(1:dim(tib_n)[1])){
  tib_n$ID2[j]  <- paste0(str_split(tib_n$variant_id[j], "_")[[1]][1],"_",str_split(tib_n$variant_id[j],"_")[[1]][2],"_",str_split(tib_n$variant_id[j],"_")[[1]][4])
}
#add to leprosy list
lep$tib.p <- tib_n$pval_nominal[match(lep$ID2,tib_n$ID2)]

#read in TMEM180 eqtl stats (skin)
skin <- read.table("gtex.Skin_Not_Sun_Exposed_Suprapubic.v8.EUR.allpairs.chr10_region.txt", header = T)
skin <- subset(skin, phenotype_id=="ENSG00000138111.14")
skin$ID2 <- NA
for (j in c(1:dim(skin)[1])){
  skin$ID2[j]  <- paste0(str_split(skin$variant_id[j], "_")[[1]][1],"_",str_split(skin$variant_id[j],"_")[[1]][2],"_",str_split(skin$variant_id[j],"_")[[1]][4])
}
#add to leprosy list
lep$skin.p <- skin$pval_nominal[match(lep$ID2,skin$ID2)]


eur.ld <- data.frame(eur.ld)
lep$r2 <- eur.ld$chr10_102467226_G[match(lep$ID2,rownames(eur.ld))]

lep$bin_r2 <- 1
lep$bin_r2[which(lep$r2>0.1 & lep$r2 <= 0.3)] <- 2
lep$bin_r2[which(lep$r2>0.3 & lep$r2 <= 0.5)] <- 3
lep$bin_r2[which(lep$r2>0.5 & lep$r2 <= 0.7)] <- 4
lep$bin_r2[which(lep$r2>0.7)] <- 5

coplot1 <- data.frame(lep$p, lep$bin_r2, lep$cd4.p)
coplot1$trait <- "ACTR1A (CD4)"
colnames(coplot1) <- c("lep.p", "bin_r2", "trait_p", "trait")

coplot2 <- data.frame(lep$p, lep$bin_r2, lep$skin.p)
coplot2$trait <- "TMEM180 (skin)"
colnames(coplot2) <- c("lep.p", "bin_r2", "trait_p", "trait")

coplot3 <- data.frame(lep$p, lep$bin_r2, lep$tib.p)
coplot3$trait <- "TMEM180 (nerve)"
colnames(coplot3) <- c("lep.p", "bin_r2", "trait_p", "trait")

coplots <- na.omit(rbind(coplot1, coplot2, coplot3))

pp4_labels = data.frame(trait = c("ACTR1A (CD4)", 
                                  "TMEM180 (skin)", 
                                  "TMEM180 (nerve)"), 
                        label1 = c("PP4==0.94", 
                                   "PP4==0.92", 
                                   "PP4==0.95"))

cols2 <- brewer.pal(11,"Spectral")
primary.gtex.coplot <- ggplot(coplots, aes(x=-log10(trait_p), y=-log10(lep.p))) + 
  geom_point(data=subset(coplots, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(coplots, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(coplots, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(coplots, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(coplots, bin_r2==5), color=cols2[1], size=3) + 
  ylab("Leprosy GWAS (-logP)") + 
  xlab("eQTL (-logP)") + 
  ylim(NA,9) +
  facet_wrap( ~ trait, ncol=3, scales = "free_x") + geom_text(x=0, y=9, aes(label=label1), data=pp4_labels, parse=TRUE, inherit.aes=F, size=4, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"), plot.title = element_text(size = 20, face = "bold"))
primary.gtex.coplot


p1 <- (primary.cell.eqtl.fp)/(gtex.eqtl.fp)/(primary.gtex.coplot)

ggsave(
  "Fig3ABC.jpg",
  width = 7,
  height = 10.5,
  dpi = 300
) 
