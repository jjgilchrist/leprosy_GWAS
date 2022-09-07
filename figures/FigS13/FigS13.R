rm(list = ls(all=TRUE))
library(ggplot2)
library(patchwork)
library(RColorBrewer)

cols <- brewer.pal(8,"Set2")

#read in genotype PCs - India samples
pcs <- read.table("india.pc.eigenvec", header = T, row.names = 1)
pcs <- na.omit(pcs)

p.pc12 = ggplot(pcs, aes(PC1, PC2)) + 
  geom_point(aes(colour = factor(cc)), size = 2) + scale_colour_manual(name = "cc", values =c(cols[4], cols[5])) +
  theme_bw() +
  ylab("PC2") +
  xlab("PC1") +
  xlim(-0.5,0.5) +
  ylim(-0.5,0.5) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.position = "none",
        plot.title=element_text(size=50, face = "bold"))
p.pc12

p.pc34 = ggplot(pcs, aes(PC3, PC4)) + 
  geom_point(aes(colour = factor(cc)), size = 2) + scale_colour_manual(name = "cc", values =c(cols[4], cols[5])) +
  theme_bw() +
  ylab("PC4") +
  xlab("PC3") +
  xlim(-0.5,0.5) +
  ylim(-0.5,0.5) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.position = "none",
        plot.title=element_text(size=50, face = "bold"))
p.pc34

p.pc56 = ggplot(pcs, aes(PC5, PC6)) + 
  geom_point(aes(colour = factor(cc)), size = 2) + scale_colour_manual(name = "cc", values =c(cols[4], cols[5])) +
  theme_bw() +
  ylab("PC6") +
  xlab("PC5") +
  xlim(-0.5,0.5) +
  ylim(-0.5,0.5) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.position = "none",
        plot.title=element_text(size=50, face = "bold"))
p.pc56

p1 <- (p.pc12 + p.pc34+ p.pc56) 


ggsave(
  "FigS13.jpg",
  width = 7,
  height = 3.5,
  dpi = 300
) 

