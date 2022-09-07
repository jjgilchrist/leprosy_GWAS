rm(list = ls(all=TRUE))

library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(patchwork)
library(edgeR)


#read in leprosy skin Bx RNASeq data (Montoya et al dataset)
data <- read.table("actr1a_skin_bx_rnaseq.txt", header = T)
#test normality of expression data
shapiro.test(log2(data$ACTR1A))

#remove reversal reactions
data.no.RR <- data[-which(data$disease=="RR"),]
t.test(log2(ACTR1A) ~ disease, data = data.no.RR)
 
#construct box plot of ACTR1A expression (quantified by RNASeq) association with leprosy subtype
cols <- brewer.pal(8,"Set2")
p_label1 = data.frame(disease = c("disease"), label = c("p=0.00068"))
rnaseq.actr1a.box = ggplot(data.no.RR, aes(x=factor(disease), y=log2(ACTR1A))) +
  geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(disease), col=factor(disease)) + scale_fill_manual(values = cols[c(1,2,3)]) + scale_colour_manual(values = cols[c(1,2,3)]) +
  ylab("Normalized ACTR1A expression") +
  scale_y_continuous(breaks=c(11.8,12.0, 12.2,12.4), labels=c(11.8,12.0, 12.2,12.4), limits = c(NA,12.5)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  geom_text(x = 1.5, y = 12.4, vjust=-1, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5)
rnaseq.actr1a.box

#read in leprosy skin Bx microarray data (Belone et al dataset)
data <- read.table("lep_array.txt", header = T)

data$facet1 <- "ACTR1A"
data$facet2 <- "TMEM180"

#test normality
shapiro.test(data$ACTR1A)
shapiro.test(data$TMEM180)
#ACTR1A normally distributed, TMEM180 not

#plot ACTR1A healthy control vs leprosy association
t.test(ACTR1A~comp1, data = data)
p_label1 = data.frame(comp1 = c("comp1"), label = c("p=0.0004"))
array.actr1a.comp1 = ggplot(data, aes(x=factor(comp1), y=ACTR1A)) +
  geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(comp1), col=factor(comp1)) + scale_fill_manual(values = cols[c(8,3)]) + scale_colour_manual(values = cols[c(8,3)]) +
  ylab("Normalized gene expression") +
  scale_y_continuous(breaks=c(-0.4,-0.2,0,0.2), labels=c(-0.4,-0.2,0,0.2), limits = c(NA,0.25)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  #geom_segment(x = 1, y = 12.6, xend = 3, yend = 12.6, colour = "black") +
  #geom_segment(x = 1, y = 12.5, xend = 2, yend = 12.5, colour = "black") +
  #geom_segment(x = 3, y = 12.4, xend = 2, yend = 12.4, colour = "black") +
  geom_text(x = 1.5, y = 0.2, vjust=-0.5, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
  facet_wrap( ~ facet1, ncol=1, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic")) #+
  #geom_text(x = 2, y = 12.6, vjust=-1, aes(label=label), data=p_label2, parse=FALSE, inherit.aes=F, size = 5) +
  #geom_text(x = 2.5, y = 12.4, vjust=-1, aes(label=label), data=p_label3, parse=FALSE, inherit.aes=F, size = 5)
array.actr1a.comp1

#plot TMEM180 healthy control vs leprosy association
wilcox.test(TMEM180~comp1, data = data)
p_label1 = data.frame(comp1 = c("comp1"), label = c("NS"))
array.tmem180.comp1 = ggplot(data, aes(x=factor(comp1), y=TMEM180)) +
  geom_dotplot(binaxis="y", binwidth=0.03, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(comp1), col=factor(comp1)) + scale_fill_manual(values = cols[c(8,3)]) + scale_colour_manual(values = cols[c(8,3)]) +
  ylab("Normalized gene expression") +
  #scale_y_continuous(breaks=c(-0.4,-0.2,0,0.2), labels=c(-0.4,-0.2,0,0.2), limits = c(NA,NA)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  #geom_segment(x = 1, y = 12.6, xend = 3, yend = 12.6, colour = "black") +
  #geom_segment(x = 1, y = 12.5, xend = 2, yend = 12.5, colour = "black") +
  #geom_segment(x = 3, y = 12.4, xend = 2, yend = 12.4, colour = "black") +
  geom_text(x = 1.5, y = 1.18, vjust=-0.5, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
  facet_wrap( ~ facet2, ncol=1, scales="free_y") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic")) #+
#geom_text(x = 2, y = 12.6, vjust=-1, aes(label=label), data=p_label2, parse=FALSE, inherit.aes=F, size = 5) +
#geom_text(x = 2.5, y = 12.4, vjust=-1, aes(label=label), data=p_label3, parse=FALSE, inherit.aes=F, size = 5)
array.tmem180.comp1

#test for association of ACTR1A expression with leprosy subtype - again in Belone et al dataset
res.aov <- aov(ACTR1A ~ comp5, data = data)
summary(res.aov)
TukeyHSD(res.aov)

p_label1 = data.frame(disease = c("disease"), label = c("p=0.00851"))
p_label2 = data.frame(disease = c("disease"), label = c("p=0.00027"))
p_label3 = data.frame(disease = c("disease"), label = c("NS"))
p_label4 = data.frame(disease = c("disease"), label = c("p=0.00001"))
p_label5 = data.frame(disease = c("disease"), label = c("p=0.00018"))
p_label6 = data.frame(disease = c("disease"), label = c("NS"))

#plot association with leprosy subtype in Belone et al array data
data1 <- data[which(!is.na(data$comp5)),]
data1$comp5 <- factor(data1$comp5, levels = c("HC", "LL", "B", "TT"))
array.actr1a.comp5 = ggplot(data1, aes(x=factor(comp5), y=ACTR1A)) +
  geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  aes(fill = factor(comp5), col=factor(comp5)) + scale_fill_manual(values = cols[c(8,1,7,2)]) + scale_colour_manual(values = cols[c(8,1,7,2)]) +
  ylab("Normalized ACTR1A expression") +
  scale_y_continuous(breaks=c(-0.4,-0.2,0.0,0.2), labels=c(-0.4,-0.2,0.0,0.2), limits = c(NA,0.55)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title.y=element_text(size=15), axis.title.x=element_blank()) +
  geom_segment(x = 1, y = 0.48, xend = 4, yend = 0.48, colour = "black") +
  geom_segment(x = 1, y = 0.38, xend = 3, yend = 0.38, colour = "black") +
  geom_segment(x = 2, y = 0.28, xend = 3, yend = 0.28, colour = "black") +
  geom_segment(x = 2, y = 0.18, xend = 4, yend = 0.18, colour = "black") +
  geom_text(x = 2.5, y = 0.48, vjust=-0.3, aes(label=label), data=p_label1, parse=FALSE, inherit.aes=F, size = 5) +
  geom_text(x = 2, y = 0.38, vjust=-0.3, aes(label=label), data=p_label2, parse=FALSE, inherit.aes=F, size = 5) +
  geom_text(x = 2.5, y = 0.28, vjust=-0.3, aes(label=label), data=p_label4, parse=FALSE, inherit.aes=F, size = 5) +
  geom_text(x = 3, y = 0.18, vjust=-0.3, aes(label=label), data=p_label5, parse=FALSE, inherit.aes=F, size = 5)
array.actr1a.comp5




p1 <- (array.actr1a.comp1|array.tmem180.comp1)/(array.actr1a.comp5)/(rnaseq.actr1a.box)

ggsave(
  "Fig3DEF.jpg",
  width = 7,
  height = 10.5,
  dpi = 300
) 

