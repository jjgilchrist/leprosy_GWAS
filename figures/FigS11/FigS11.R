rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(patchwork)

#read in data - Malawi
pca = read.table("malawi_pca.evec", header = T)
cols <- brewer.pal(12, "Set3")

eth.pc12 = ggplot(pca, aes(PC2, PC1)) + 
  geom_point(data=subset(pca, eth==1), color=cols[1], size=4) +
  geom_point(data=subset(pca, eth==2), color=cols[2], size=4) +
  geom_point(data=subset(pca, eth==3), color=cols[3], size=4) +
  geom_point(data=subset(pca, eth==4), color=cols[4], size=4) +
  geom_point(data=subset(pca, eth==5), color=cols[5], size=4) +
  geom_point(data=subset(pca, eth==6), color=cols[6], size=4) +
  theme_bw() +
  xlab("PC2") +
  ylab("PC1") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
eth.pc34 = ggplot(pca, aes(PC4, PC3)) + 
  geom_point(data=subset(pca, eth==1), color=cols[1], size=4) +
  geom_point(data=subset(pca, eth==2), color=cols[2], size=4) +
  geom_point(data=subset(pca, eth==3), color=cols[3], size=4) +
  geom_point(data=subset(pca, eth==4), color=cols[4], size=4) +
  geom_point(data=subset(pca, eth==5), color=cols[5], size=4) +
  geom_point(data=subset(pca, eth==6), color=cols[6], size=4) +
  theme_bw() +
  xlab("PC4") +
  ylab("PC3") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
eth.pc56 = ggplot(pca, aes(PC6, PC5)) + 
  geom_point(data=subset(pca, eth==1), color=cols[1], size=4) +
  geom_point(data=subset(pca, eth==2), color=cols[2], size=4) +
  geom_point(data=subset(pca, eth==3), color=cols[3], size=4) +
  geom_point(data=subset(pca, eth==4), color=cols[4], size=4) +
  geom_point(data=subset(pca, eth==5), color=cols[5], size=4) +
  geom_point(data=subset(pca, eth==6), color=cols[6], size=4) +
  theme_bw() +
  xlab("PC6") +
  ylab("PC5") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
eth.pc78 = ggplot(pca, aes(PC8, PC7)) + 
  geom_point(data=subset(pca, eth==1), color=cols[1], size=4) +
  geom_point(data=subset(pca, eth==2), color=cols[2], size=4) +
  geom_point(data=subset(pca, eth==3), color=cols[3], size=4) +
  geom_point(data=subset(pca, eth==4), color=cols[4], size=4) +
  geom_point(data=subset(pca, eth==5), color=cols[5], size=4) +
  geom_point(data=subset(pca, eth==6), color=cols[6], size=4) +
  theme_bw() +
  xlab("PC8") +
  ylab("PC7") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
eth.pc910 = ggplot(pca, aes(PC10, PC9)) + 
  geom_point(data=subset(pca, eth==1), color=cols[1], size=4) +
  geom_point(data=subset(pca, eth==2), color=cols[2], size=4) +
  geom_point(data=subset(pca, eth==3), color=cols[3], size=4) +
  geom_point(data=subset(pca, eth==4), color=cols[4], size=4) +
  geom_point(data=subset(pca, eth==5), color=cols[5], size=4) +
  geom_point(data=subset(pca, eth==6), color=cols[6], size=4) +
  theme_bw() +
  xlab("PC10") +
  ylab("PC9") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

cc.pc12 = ggplot(pca, aes(PC2, PC1)) + 
  geom_point(data=subset(pca, cc=="CASE"), color=cols[7], size=4) +
  geom_point(data=subset(pca, cc=="CONTROL"), color=cols[8], size=4) +
  theme_bw() +
  xlab("PC2") +
  ylab("PC1") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
cc.pc34 = ggplot(pca, aes(PC4, PC3)) + 
  geom_point(data=subset(pca, cc=="CASE"), color=cols[7], size=4) +
  geom_point(data=subset(pca, cc=="CONTROL"), color=cols[8], size=4) +
  theme_bw() +
  xlab("PC4") +
  ylab("PC3") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
cc.pc56 = ggplot(pca, aes(PC6, PC5)) + 
  geom_point(data=subset(pca, cc=="CASE"), color=cols[7], size=4) +
  geom_point(data=subset(pca, cc=="CONTROL"), color=cols[8], size=4) +
  theme_bw() +
  xlab("PC6") +
  ylab("PC5") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
cc.pc78 = ggplot(pca, aes(PC8, PC7)) + 
  geom_point(data=subset(pca, cc=="CASE"), color=cols[7], size=4) +
  geom_point(data=subset(pca, cc=="CONTROL"), color=cols[8], size=4) +
  theme_bw() +
  xlab("PC8") +
  ylab("PC7") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
cc.pc910 = ggplot(pca, aes(PC10, PC9)) + 
  geom_point(data=subset(pca, cc=="CASE"), color=cols[7], size=4) +
  geom_point(data=subset(pca, cc=="CONTROL"), color=cols[8], size=4) +
  theme_bw() +
  xlab("PC10") +
  ylab("PC9") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

plots.out <- (eth.pc12|eth.pc34|eth.pc56|eth.pc78|eth.pc910)/(cc.pc12|cc.pc34|cc.pc56|cc.pc78|cc.pc910)

ggsave(
  "FigS11.jpg",
  width = 21,
  height = 10.5,
  dpi = 300
) 