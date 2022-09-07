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


#read in sub-sampled association statistics for Malawi GWAS, Mali GWAS and meta-analysis. Full summary statistics are deposited with NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics; accession codes: Malawi, GCST90129399; Mali, GCST90129400; meta-analysis, GCST90129401).
#n.b. lambdas printed on these figures refer to those derived from the complete dataset. The QQ plots represent the downsampled set of summary statistics.

meta.add <- read.table("meta.sub.txt", header = T)
mali.add <- read.table("mali.sub.txt", header = T)
malawi.add <- read.table("malawi.sub.txt", header = T)


cols3 <- brewer.pal(8,"Paired")
df <- data.frame(observed = -log10(sort(na.omit(meta.add$pval))),
    expected = -log10(ppoints(length(na.omit(meta.add$pval)))))

lambda <- paste("lambda == ", 1.0123)
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  meta.qq <- ggplot(df) +
    geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
    geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
    xlab(log10Pe) +
    ylab(log10Po) +
    ggtitle("Meta-analysis") +
    annotate("text", x = 1, y = 9, label = lambda, parse = TRUE, size = 10) +
    theme_bw() + ylim(NA, 10) + 
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
      axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))


df <- data.frame(observed = -log10(sort(na.omit(malawi.add$pval))),
                 expected = -log10(ppoints(length(na.omit(malawi.add$pval)))))

lambda <- paste("lambda == ", 1.0333)

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

malawi.qq <- ggplot(df) +
  geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
  geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
  xlab(log10Pe) +
  ylab(log10Po) +
  ggtitle("Malawi") +
  annotate("text", x = 1, y = 9, label = lambda, parse = TRUE, size = 10) +
  theme_bw() + ylim(NA, 10) + 
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))

df <- data.frame(observed = -log10(sort(na.omit(as.numeric(as.character(mali.add$pval))))),
                 expected = -log10(ppoints(length(na.omit((as.numeric(as.character(mali.add$pval))))))))


lambda <- paste("lambda == ", 1.0497)

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

mali.qq <- ggplot(df) +
  geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
  geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
  xlab(log10Pe) +
  ylab(log10Po) +
  ggtitle("Mali") +
  annotate("text", x = 1, y = 9, label = lambda, parse = TRUE, size = 10) +
  theme_bw() + ylim(NA, 10) + 
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))


total.qq <- (malawi.qq|mali.qq|meta.qq)


ggsave(
  "FigS1.jpg",
  width = 21,
  height = 7,
  dpi = 300
)


