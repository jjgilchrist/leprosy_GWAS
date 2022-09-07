rm(list = ls(all=TRUE))

library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(patchwork)
library(edgeR)


#read in RNASeq data describing gene expression in whole blood in healthy household contacts and leprosy patients pre- and post-diagnosis
#data available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163498

exp <- read.table("geo_fragments_per_gene.tsv", header = T, row.names = 1)
samples <- read.table("geo_samples.txt", header = T)


exp2 <- exp[-c(1:5),]

samples$mapped_reads <- colSums(exp[c(6:58740),])
samples$unmapped_reads <- colSums(exp[c(1:5),])
samples$total_reads <- colSums(exp[c(1:5),])+colSums(exp[c(6:58740),])


#exclusions
exp.data <- DGEList(counts=exp2)
cpm <- cpm(exp.data)
lcpm <- cpm(exp.data, log=TRUE)

L <- mean(exp.data$samples$lib.size) * 1e-6
M <- median(exp.data$samples$lib.size) * 1e-6


#remove low expressed genes
keep.exprs <- filterByExpr(exp.data)
exp.data <- exp.data[keep.exprs,, keep.lib.sizes=FALSE]


#brings down to 14460 genes

#Normalise gene expression
exp.data <- calcNormFactors(exp.data, method = "TMM")
exp.data$samples$norm.factors

lcpm <- cpm(exp.data, log=TRUE)



lcpm.order <- lcpm[,as.character(samples$Sample_name)]

samples$actr1a <- lcpm.order["ENSG00000138107",]
samples$tmem180 <- lcpm.order["ENSG00000138111",]

#test normality
shapiro.test(samples$actr1a)
shapiro.test(samples$tmem180)

res.aov <- aov(actr1a ~ Group, data = samples)
summary(res.aov)
kruskal.test(tmem180 ~ Group, data = samples)


samples$Group <- factor(samples$Group, levels = c("HHC", "First", "Second"))
samples$facet1 <- "ACTR1A"
samples$facet2 <- "TMEM180"
p_label1 = data.frame(Group = c("Group"), label = c("NS"))

rnaseq.actr1a.wb = ggplot(samples, aes(x=factor(Group), y=actr1a)) +
  geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(Group), col=factor(Group)) + scale_fill_manual(values = cols[c(8,3,3)]) + scale_colour_manual(values = cols[c(8,3,3)]) +
  ylab("Normalized gene expression") +
  #scale_y_continuous(breaks=c(-0.4,-0.2,0.0,0.2), labels=c(-0.4,-0.2,0.0,0.2), limits = c(NA,0.5)) +
  scale_x_discrete(breaks=c("HHC", "First", "Second"), labels=c("HHC", "Leprosy\npre-Dx", "Leprosy\npost-Dx")) +
  ylim(NA,8.2) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  geom_text(x = 2, y = 8, vjust=-1, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
  facet_wrap( ~ facet1, ncol=1, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic")) #+
#geom_segment(x = 1, y = 0.48, xend = 4, yend = 0.48, colour = "black") +
#geom_segment(x = 1, y = 0.42, xend = 3, yend = 0.42, colour = "black") +
#geom_segment(x = 1, y = 0.36, xend = 2, yend = 0.36, colour = "black") +
#geom_segment(x = 2, y = 0.3, xend = 3, yend = 0.3, colour = "black") +
#geom_segment(x = 2, y = 0.24, xend = 4, yend = 0.24, colour = "black") +
#geom_segment(x = 3, y = 0.18, xend = 4, yend = 0.18, colour = "black") +    
#geom_text(x = 2.5, y = 0.48, vjust=-0.3, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 2, y = 0.42, vjust=-0.3, aes(label=label), data=p_label2, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 1.5, y = 0.36, vjust=-0.3, aes(label=label), data=p_label3, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 2.5, y = 0.3, vjust=-0.3, aes(label=label), data=p_label4, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 3, y = 0.24, vjust=-0.3, aes(label=label), data=p_label5, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 3.5, y = 0.18, vjust=-0.3, aes(label=label), data=p_label6, parse=FALSE, inherit.aes=F, size = 5)
rnaseq.actr1a.wb


rnaseq.tmem180.wb = ggplot(samples, aes(x=factor(Group), y=tmem180)) +
  geom_dotplot(binaxis="y", binwidth=0.1, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(Group), col=factor(Group)) + scale_fill_manual(values = cols[c(8,3,3)]) + scale_colour_manual(values = cols[c(8,3,3)]) +
  ylab("Normalized gene expression") +
  #scale_y_continuous(breaks=c(-0.4,-0.2,0.0,0.2), labels=c(-0.4,-0.2,0.0,0.2), limits = c(NA,0.5)) +
  scale_x_discrete(breaks=c("HHC", "First", "Second"), labels=c("HHC", "Leprosy\npre-Dx", "Leprosy\npost-Dx")) +
  ylim(0,6.5) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  geom_text(x = 2, y = 6.2, vjust=-0.5, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
  facet_wrap( ~ facet2, ncol=1, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic")) #+
#geom_segment(x = 1, y = 0.48, xend = 4, yend = 0.48, colour = "black") +
#geom_segment(x = 1, y = 0.42, xend = 3, yend = 0.42, colour = "black") +
#geom_segment(x = 1, y = 0.36, xend = 2, yend = 0.36, colour = "black") +
#geom_segment(x = 2, y = 0.3, xend = 3, yend = 0.3, colour = "black") +
#geom_segment(x = 2, y = 0.24, xend = 4, yend = 0.24, colour = "black") +
#geom_segment(x = 3, y = 0.18, xend = 4, yend = 0.18, colour = "black") +    
#geom_text(x = 2.5, y = 0.48, vjust=-0.3, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 2, y = 0.42, vjust=-0.3, aes(label=label), data=p_label2, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 1.5, y = 0.36, vjust=-0.3, aes(label=label), data=p_label3, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 2.5, y = 0.3, vjust=-0.3, aes(label=label), data=p_label4, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 3, y = 0.24, vjust=-0.3, aes(label=label), data=p_label5, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 3.5, y = 0.18, vjust=-0.3, aes(label=label), data=p_label6, parse=FALSE, inherit.aes=F, size = 5)
rnaseq.tmem180.wb



samples

p1 <- (rnaseq.actr1a.wb|rnaseq.tmem180.wb)

ggsave(
  "FigS3.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 