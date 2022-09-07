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



#read in rab32 region summary statistics
total <- read.table("rab32.region", header = T, sep = "\t")
#read in local LD structure
ld <- read.table("rab32.malawi.ld", header = T, sep = "\t")
ld.mali <- read.table("rab32.mali.ld", header = T, sep = "\t")

total$ld <- ld$R2[match(total$rsid2, ld$SNP_B)]
total$dp <- ld$DP[match(total$rsid2, ld$SNP_B)]
total$ld2 <- ld.mali$R2[match(total$rsid2, ld.mali$SNP_B)]
total$dp2 <- ld.mali$DP[match(total$rsid2, ld.mali$SNP_B)]
total$r2 <- total$ld
total$r2.2 <- total$ld2


cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"Spectral")

#split LD to peak SNP into bins
total$bin_r2 <- 1
total$bin_r2[which(total$r2>0.1 & total$r2 <= 0.3)] <- 2
total$bin_r2[which(total$r2>0.3 & total$r2 <= 0.5)] <- 3
total$bin_r2[which(total$r2>0.5 & total$r2 <= 0.7)] <- 4
total$bin_r2[which(total$r2>0.7)] <- 5

total$bin_r2.2 <- 1
total$bin_r2.2[which(total$r2.2>0.1 & total$r2.2 <= 0.3)] <- 2
total$bin_r2.2[which(total$r2.2>0.3 & total$r2.2 <= 0.5)] <- 3
total$bin_r2.2[which(total$r2.2>0.5 & total$r2.2 <= 0.7)] <- 4
total$bin_r2.2[which(total$r2.2>0.7)] <- 5

#highlight genotyped markers
total$genotyped <- 0
total$genotyped[grep("JHU", total$rsid2)] <- 1

#highlight peak SNP
total$annotate <- 0
total$annotate[c(which(total$rsid=="rs34271799"))] <- 1

#plot regional association plot at rab32 region
cols3 <- brewer.pal(8,"Paired")
rab32_plot <- ggplot(total, aes(x=position, y=-log10(FixedEffectMetaAnalysis.pvalue))) + 
  xlim(146700000, 147150000) +
  geom_point(data=subset(total, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(total, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(total, bin_r2==5), color=cols2[1], size=3) + 
  geom_point(data=subset(total, bin_r2.2==1), color=cols[8], size=1.5) + 
  geom_point(data=subset(total, bin_r2.2==2), color=cols2[5], size=1.5) + 
  geom_point(data=subset(total, bin_r2.2==3), color=cols2[3], size=1.5) + 
  geom_point(data=subset(total, bin_r2.2==4), color=cols2[2], size=1.5) + 
  geom_point(data=subset(total, bin_r2.2==5), color=cols2[1], size=1.5) + 
  geom_point(data=subset(total, genotyped==1), color="black", size=1, shape=3) + 
  ylab("-log P-value") + 
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(NA, 6) +
  geom_text_repel( data=subset(total, annotate==1), aes(label=rsid2), size=6, col = c("black"), nudge_y = 0.5) +
  annotate("text", x = 146720000, y = 6.0, label = c("Malawi")) +
  annotate("text", x = 146780000, y = 6.0, label = c("Mali")) +
  annotate("text", x = 146750000, y = 6.0, label = quote(r^2), hjust = 0.5) +
  annotate("point", x = 146720000, y = 5.7, size = 3, colour = cols2[1]) +
  annotate("point", x = 146780000, y = 5.7, size = 1.5, colour = cols2[1]) +
  annotate("text", x = 146750000, y = 5.7, label = c("0.7-1.0"), hjust = 0.5) +
  annotate("point", x = 146720000, y = 5.4, size = 3, colour = cols2[2]) +
  annotate("point", x = 146780000, y = 5.4, size = 1.5, colour = cols2[2]) +
  annotate("text", x = 146750000, y = 5.4, label = c("0.5-0.7"), hjust = 0.5) +
  annotate("point", x = 146720000, y = 5.1, size = 3, colour = cols2[3]) +
  annotate("point", x = 146780000, y = 5.1, size = 1.5, colour = cols2[3]) +
  annotate("text", x = 146750000, y = 5.1, label = c("0.3-0.5"), hjust = 0.5) +
  annotate("point", x = 146720000, y = 4.8, size = 3, colour = cols2[5]) +
  annotate("point", x = 146780000, y = 4.8, size = 1.5, colour = cols2[5]) +
  annotate("text", x = 146750000, y = 4.8, label = c("0.1-0.3"), hjust = 0.5) +
  annotate("point", x = 146720000, y = 4.5, size = 3, colour = cols[8]) +
  annotate("point", x = 146780000, y = 4.5, size = 1.5, colour = cols[8]) +
  annotate("text", x = 146750000, y = 4.5, label = c("<0.1"), hjust = 0.5)


#genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=6
sel.pos=146900000
range=2500000

listAttributes(gene.ensembl)

out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2

genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding")

genes$start <- genes$start_position
genes$end <- genes$end_position

genes$start[which(genes$strand==-1)] <- genes$end_position[which(genes$strand==-1)]
genes$end[which(genes$strand==-1)] <- genes$start_position[which(genes$strand==-1)]

plot.range <- c(146700000, 147150000)
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(146700000, 147150000) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes,
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[c(7:9),], aes(x=mid, label=external_gene_name), size=4, col = c("black"),
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


recomb <- read.table("genetic_map_chr6_b37.txt", header = T)
recomb_rab32 <- subset(recomb, position>146700000 & position<147150000)

recomb_rate <- ggplot(recomb_rab32, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 6") +
  scale_x_continuous(breaks=c(146800000,146900000,147000000,147100000),
                     labels=c("146.8Mb", "146.9Mb", "147.0Mb", "147.1Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


p1 <- (rab32_plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 0.5, 0.5))

ggsave(
  "FigS7C.jpg",
  width = 7,
  height = 5,
  dpi = 300
) 



#plot forest plot of rs34271799 association in Malawi, Mali and in meta-analysis
cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.beta_1.add.cc.1

label <- rep(c("Malawi", "Mali", "Meta"),3)
mean  <- c(total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$FixedEffectMetaAnalysis.beta_1)
lower <- c(total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$FixedEffectMetaAnalysis.beta_1)-(1.96*c(total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.se_1, total[which(total$rsid2=="rs34271799:CT"),]$cohort.1.se_1, total[which(total$rsid2=="rs34271799:CT"),]$FixedEffectMetaAnalysis.se_1))
upper <- c(total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid2=="rs34271799:CT"),]$FixedEffectMetaAnalysis.beta_1)+(1.96*c(total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.se_1, total[which(total$rsid2=="rs34271799:CT"),]$cohort.1.se_1, total[which(total$rsid2=="rs34271799:CT"),]$FixedEffectMetaAnalysis.se_1))

df <- data.frame(label, mean, lower, upper)
df$facet <- "RAB32:rs34271799"

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:3)])
rab32.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  ylim(-2,1) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols2[c(3)], cols[c(3,3)])) + scale_colour_manual(values = c(cols2[c(3)], cols[c(3,3)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15))+
  facet_wrap(~facet, ncol=1)+
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
rab32.fp


#plot forest plot of rs2275606 association in Malawi, Mali and in meta-analysis, and in China
cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")
total[which(total$rsid2=="rs34271799:CT"),]$cohort.2.beta_1.add.cc.1

z=-0.862+sqrt((0.743-2.404*log(3.94e-14)))
beta=log(1.3)
se=beta/z

label <- rep(c("Malawi", "Mali", "Meta: Africa", "China"),3)
mean  <- c(total[which(total$rsid=="rs2275606"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$FixedEffectMetaAnalysis.beta_1, beta)
lower <- c(total[which(total$rsid=="rs2275606"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$FixedEffectMetaAnalysis.beta_1, beta)-(1.96*c(total[which(total$rsid=="rs2275606"),]$cohort.2.se_1, total[which(total$rsid=="rs2275606"),]$cohort.1.se_1, total[which(total$rsid=="rs2275606"),]$FixedEffectMetaAnalysis.se_1, se))
upper <- c(total[which(total$rsid=="rs2275606"),]$cohort.2.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$cohort.1.beta_1.add.cc.1, total[which(total$rsid=="rs2275606"),]$FixedEffectMetaAnalysis.beta_1, beta)+(1.96*c(total[which(total$rsid=="rs2275606"),]$cohort.2.se_1, total[which(total$rsid=="rs2275606"),]$cohort.1.se_1, total[which(total$rsid=="rs2275606"),]$FixedEffectMetaAnalysis.se_1, se))

df <- data.frame(label, mean, lower, upper)
df$facet <- "RAB32:rs2275606"

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:4)])
rab32.fp.orig <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  ylim(-1.5,1.5) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[8], cols2[c(3)], cols[c(3,3)])) + scale_colour_manual(values = c(cols[8], cols2[c(3)], cols[c(3,3)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("log(OR)") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15))+
  facet_wrap(~facet, ncol=1)+
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))
rab32.fp.orig


p1 <- (rab32.fp.orig|rab32.fp)

ggsave(
  "FigS7AB.jpg",
  width = 7,
  height = 3.5,
  dpi = 300
) 
