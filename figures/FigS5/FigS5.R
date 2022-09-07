rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(biomaRt)

#download regional gene positions from biomart
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=6
sel.pos=31500000
range=2500000

out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2

hla.genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding")
hla.genes2 <- hla.genes[c(grep("HLA", hla.genes$external_gene_name), 
                          which(hla.genes$external_gene_name=="TNF"), 
                          which(hla.genes$external_gene_name=="HCG26"), 
                          which(hla.genes$external_gene_name=="C4A"), 
                          which(hla.genes$external_gene_name=="C4B"),
                          which(hla.genes$external_gene_name=="C2")), ]


hla.genes2$start <- hla.genes2$start_position
hla.genes2$end <- hla.genes2$end_position

hla.genes2$start[which(hla.genes2$strand==-1)] <- hla.genes2$end_position[which(hla.genes2$strand==-1)]
hla.genes2$end[which(hla.genes2$strand==-1)] <- hla.genes2$start_position[which(hla.genes2$strand==-1)]

plot.range <- c(29500000,33500000)
hla.genes2$order <- rep(seq(1:6),100)[c(1:length(hla.genes2$end_position))]

#read in leprosy association statistics at HLA (includes stats conditionned on HLA-DQB1)
total <- read.table("add_total.txt", header = T)
#read in 1000GP genetic map
recomb <- read.table("genetic_map_chr6_b37.txt", header = T)
#read in Mali LD structure at HLA - total already contains Malawi LD structure to peak (r2)
mali.ld2 <- read.table("hla_cond1.mali.ld", header = T)
mali.ld2$SNP_B <- as.character(mali.ld2$SNP_B)
total$rsid <- as.character(total$rsid)
total$mali.ld2 <- mali.ld$R2[match(total$rsid, mali.ld2$SNP_B)]


cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"Spectral")

#SNPs/alleles to highlight
total$annotate <- 0
total$annotate[c(which(total$rsid=="rs2516438:C"), which(total$rsid=="HLAB_4901"))] <- 1

total$bin_r2_cond <- 1
total$bin_r2_cond[which(total$r2_cond>0.1 & total$r2_cond <= 0.3)] <- 2
total$bin_r2_cond[which(total$r2_cond>0.3 & total$r2_cond <= 0.5)] <- 3
total$bin_r2_cond[which(total$r2_cond>0.5 & total$r2_cond <= 0.7)] <- 4
total$bin_r2_cond[which(total$r2_cond>0.7)] <- 5

total$bin_r2_cond2 <- 1
total$bin_r2_cond2[which(total$mali.ld2>0.1 & total$mali.ld2 <= 0.3)] <- 2
total$bin_r2_cond2[which(total$mali.ld2>0.3 & total$mali.ld2 <= 0.5)] <- 3
total$bin_r2_cond2[which(total$mali.ld2>0.5 & total$mali.ld2 <= 0.7)] <- 4
total$bin_r2_cond2[which(total$mali.ld2>0.7)] <- 5

total$genotyped <- 0

#highlight genotyped SNPs
total$genotyped[grep("JHU", total$rsid)] <- 1


cols3 <- brewer.pal(8,"Paired")

#plot conditional HLA association
cond_hla <- ggplot(total, aes(x=position, y=-log10(pval_cond))) + 
  xlim(29900000, 33500000) +
  geom_point(data=subset(total, bin_r2_cond==1), color=cols[8], size=3) + 
  geom_point(data=subset(total, bin_r2_cond==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total, bin_r2_cond==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total, bin_r2_cond==4), color=cols2[2], size=3) + 
  geom_point(data=subset(total, bin_r2_cond==5), color=cols2[1], size=3) + 
  geom_point(data=subset(total, bin_r2_cond2==1), color=cols[8], size=1.5) + 
  geom_point(data=subset(total, bin_r2_cond2==2), color=cols2[5], size=1.5) + 
  geom_point(data=subset(total, bin_r2_cond2==3), color=cols2[3], size=1.5) + 
  geom_point(data=subset(total, bin_r2_cond2==4), color=cols2[2], size=1.5) + 
  geom_point(data=subset(total, bin_r2_cond2==5), color=cols2[1], size=1.5) + 
  geom_point(data=subset(total, genotyped==1), color="black", size=1, shape=3) + 
  geom_point(data=subset(total, hla==1), color="black", size=4, shape=5, stroke = 1.2) + 
  geom_point(data=subset(total, hla==1 & pval_cond<0.0005), color=cols3[2], size=4, shape=5, stroke = 1.2) + 
  ylab("-log P-value") + 
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(NA, 7) +
  geom_text_repel( data=subset(total, annotate==1), aes(label=variant), size=6, col = c("black"), nudge_y = 0.5) +
  annotate("text", x = 29900000, y = 7.0, label = c("Malawi")) +
  annotate("text", x = 30500000, y = 7.0, label = c("Mali")) +
  annotate("text", x = 30200000, y = 7.0, label = quote(r^2), hjust = 0.5) +
  annotate("point", x = 29900000, y = 6.6, size = 3, colour = cols2[1]) +
  annotate("point", x = 30500000, y = 6.6, size = 1.5, colour = cols2[1]) +
  annotate("text", x = 30200000, y = 6.6, label = c("0.7-1.0"), hjust = 0.5) +
  annotate("point", x = 29900000, y = 6.2, size = 3, colour = cols2[2]) +
  annotate("point", x = 30500000, y = 6.2, size = 1.5, colour = cols2[2]) +
  annotate("text", x = 30200000, y = 6.2, label = c("0.5-0.7"), hjust = 0.5) +
  annotate("point", x = 29900000, y = 5.8, size = 3, colour = cols2[3]) +
  annotate("point", x = 30500000, y = 5.8, size = 1.5, colour = cols2[3]) +
  annotate("text", x = 30200000, y = 5.8, label = c("0.3-0.5"), hjust = 0.5) +
  annotate("point", x = 29900000, y = 5.4, size = 3, colour = cols2[5]) +
  annotate("point", x = 30500000, y = 5.4, size = 1.5, colour = cols2[5]) +
  annotate("text", x = 30200000, y = 5.4, label = c("0.1-0.3"), hjust = 0.5) +
  annotate("point", x = 29900000, y = 5.1, size = 3, colour = cols[8]) +
  annotate("point", x = 30500000, y = 5.1, size = 1.5, colour = cols[8]) +
  annotate("text", x = 30200000, y = 5.1, label = c("<0.1"), hjust = 0.5) +
  annotate("text", x = 33500000, y = 7.0, label = c("Conditioned: HLA-DQB1"), hjust = 1, size = 7)

#plot location of genes
genes <- ggplot(hla.genes2, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(29900000, 33500000) +
  ylim(c(0,7.5)) +
  geom_segment(data = hla.genes2,
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = hla.genes2, aes(label=external_gene_name), size=2, col = c("black"),
                   nudge_x = 30000, segment.color = NA) +
  annotate("segment", x = 29909037, xend = 31324965, y = 0.5, yend = 0.5,
           colour = "black") +
  annotate("text", x = (31324965-29909037)/2+29909037, y = 0, label = c("MHC Class I")) +
  annotate("segment", x = 31543344, xend = 32003195, y = 0.5, yend = 0.5,
           colour = "black") +
  annotate("text", x = (32003195-31543344)/2+31543344, y = 0, label = c("MHC Class III")) +
  annotate("segment", x = 32407619, xend = 33054978, y = 0.5, yend = 0.5,
           colour = "black") +
  annotate("text", x = (33054978-32407619)/2+32407619, y = 0, label = c("MHC Class II")) +
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



#plot recombination
recomb_hla <- subset(recomb, position>29900000 & position<33500000)

recomb_rate <- ggplot(recomb_hla, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 6") +
  scale_x_continuous(breaks=c(30000000,31000000,32000000,33000000),
                   labels=c("30.0Mb", "31.0Mb", "32.0Mb", "33.0Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

p1 <- (cond_hla/genes/recomb_rate) + plot_layout(heights = c(3, 1, 1))

ggsave(
  "FigS5.jpg",
  width = 9,
  height = 7,
  dpi = 300
)


