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
meta.add <- read.table("Fig1B_data.txt", header = T)

meta.add$chr <- factor(meta.add$chr, levels = c(1:22))
levels(meta.add$chr) <- c(1:22)

#highlight HLA and chr10:
meta.add$label <- NA
meta.add$label[which(meta.add$rsid=="rs2015583:G")] <- "10q24.32"
meta.add$label[which(meta.add$rsid=="rs9270296:T")] <- "HLA"

meta.add$anno <- NA
meta.add$anno[which(meta.add$rsid=="rs2015583:G")] <- 1
meta.add$anno[which(meta.add$rsid=="rs9270296:T")] <- 1

don.1 <- meta.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(meta.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = don.1 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

meta.manh <- ggplot(don.1, aes(x=BPcum, y=-log10(pval))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  
  # Add test using ggrepel to avoid overlapping
  geom_text_repel( data=subset(don.1, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = don.1, linetype="dashed", color = "red") +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","2","4","6","8","10"), breaks = c(0,2,4,6,8,10), expand = c(0, 0), limits= c(0,10)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

p1 <- meta.manh

ggsave(
  "Fig1B.jpg",
  width = 9.5,
  height = 4.5,
  dpi = 300
)
