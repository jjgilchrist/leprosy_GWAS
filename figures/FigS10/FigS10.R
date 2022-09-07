rm(list = ls(all=TRUE))

library(ggplot2)
library(patchwork)
library(RColorBrewer)

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#read in data
pca <- read.table("malawi_mali_agvp.txt", header = TRUE)

#plot first 2 PCs of genotyping data - highlighting study samples and outliers
p = ggplot(pca, aes(PC1, PC2))+ 
  geom_point(data=subset(pca, study=="AGVP"), size=4, col = cols[8]) +
  geom_point(data=subset(pca, outlier==1), size=6, col = "black") +
  geom_point(data=subset(pca, study=="Mali"), size=4, col = cols[1]) +
  geom_point(data=subset(pca, study=="Malawi"), size=4, col = cols[2]) +
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25),
        legend.position = "none")
p

ggsave(
  "FigS10.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 