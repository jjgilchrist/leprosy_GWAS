rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(ggbio)
library(data.table)
library(GenomicRanges)

#read in data - all SNPs with p<10e-5, remainder of SNPs down-sampled to 100,000 - full summary statistics are deposited with NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics; accession codes: Malawi, GCST90129399; Mali, GCST90129400; meta-analysis, GCST90129401).
mali.add <- read.table("mali.sub.txt", header = T)
malawi.add <- read.table("malawi.sub.txt", header = T)

mali.add$chr <- factor(mali.add$chr, levels = c(1:22))
malawi.add$chr <- factor(malawi.add$chr, levels = c(1:22))

levels(mali.add$chr) <- c(1:22)
levels(malawi.add$chr) <- c(1:22)

mali.add$bp <- as.numeric(as.character(mali.add$bp))
malawi.add$bp <- as.numeric(as.character(malawi.add$bp))
mali.add$pval <- as.numeric(as.character(mali.add$pval))
malawi.add$pval <- as.numeric(as.character(malawi.add$pval))


don.2 <- mali.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(mali.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

don.3 <- malawi.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(malawi.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)


axisdf = don.2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

mali.manh <- ggplot(don.2, aes(x=BPcum, y=-log10(pval))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","2","4","6","8","10"), breaks = c(0,2,4,6,8,10), expand = c(0, 0), limits= c(0,10)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  #title
  ggtitle("Mali") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20),
    plot.title = element_text(size = 25, face = "bold")
  )

malawi.manh <- ggplot(don.3, aes(x=BPcum, y=-log10(pval))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","2","4","6","8","10"), breaks = c(0,2,4,6,8,10), expand = c(0, 0), limits= c(0,10)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = don.3, linetype="dashed", color = "red") +
  geom_segment(aes(x = 10177, y = 5, xend = 2879943885, yend = 5), data = don.3, linetype="dotted", color = "red") +
  
  
  #title
  ggtitle("Malawi") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20),
    plot.title = element_text(size = 25, face = "bold")
  )


p1 <- (malawi.manh/mali.manh)

ggsave(
  "FigS2.jpg",
  width = 14,
  height = 7,
  dpi = 300
)




