# leprosy_GWAS

**Genome-wide association study of leprosy in Malawi and Mali.** (PLoS Pathogens, 2022)

https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010312

Code and source data for GWAS of leprosy in Malawi and Mali.

**Dataset**
612 Malawian participants: 284 leprosy cases, 328 controls.
350 Malian participants: 208 leprosy cases, 142 controls.

**Abstract**
Leprosy is a chronic infection of the skin and peripheral nerves caused by *Mycobacterium leprae*. Despite recent improvements in disease control, leprosy remains an important cause of infectious disability globally. Large-scale genetic association studies in Chinese, Vietnamese and Indian populations have identified over 30 susceptibility loci for leprosy. There is a significant burden of leprosy in Africa, however it is uncertain whether the findings of published genetic association studies are generalizable to African populations. To address this, we conducted a genome-wide association study (GWAS) of leprosy in Malawian and Malian individuals. In that analysis, we replicated four risk loci previously reported in China, Vietnam and India; MHC Class I and II, *LACC1* and *SLC29A3*. We further identified a novel leprosy susceptibility locus at 10q24 (rs2015583; combined p=8.81x10<sup>-9</sup>; OR=0.51 [95% CI 0.40-0.64]). Using publicly-available data we characterise regulatory activity at this locus, identifying *ACTR1A* as a candidate mediator of leprosy risk. This locus shows evidence of recent positive selection and demonstrates pleiotropy with established risk loci for inflammatory bowel disease and childhood-onset asthma. A shared genetic architecture for leprosy and inflammatory bowel disease has been previously described. We expand on this, strengthening the hypothesis that selection pressure driven by leprosy has shaped the evolution of autoimmune and atopic disease in modern populations. More broadly, our data highlights the importance of defining the genetic architecture of disease across genetically diverse populations, and that disease insights derived from GWAS in one population may not translate to all affected populations.

**Overview of repository**
* Figure folders: contains R script and source data to reproduce each main and supplementary Figure from the manuscript.
* Leprosy association analysis: contains scripts and genotype/phenotype data to perform association analysis and meta-analysis to reproduce the study's main findings. Association analysis was performed using SNPTEST v2.5.6: [Marchini J, *et al*. "A new multipoint method
for genome-wide association studies by imputation of genotypes." *Nature Genetics.* 2009.](https://doi.org/10.1038/ng2088) Meta-analysis was performed using bingwa v2.0.6: [Malaria Genomic Epidemiology Network. "A novel locus of resistance to severe malaria in a region of ancient balancing selection." *Nature.* 2015.](https://doi.org/10.1038/nature15390)

**Data availability and sources**
Complete GWAS summary statistics are available through the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics; accession codes: Malawi, GCST90129399; Mali, GCST90129400; meta-analysis, GCST90129401). Patient level genotype and phenotype data for additional MalariaGEN control samples (Mali) are available are available via the European Genome-Phenome Archive, with accession code EGAS00001001311.

Colocalisation analysis makes use of publically-available eQTL mapping data:
* Genotype-Tissue Expression (GTEx) Project V8 [GTEx Consortium. "The GTEx Consortium atlas of genetic regulatory effects across human tissues." *Science.* 2020.](https://doi.org/10.1126/science.aaz1776)
* T cells [Kasela S, *et al*. "Pathogenic implications for autoimmune mechanisms derived by comparative eQTL analysis of CD4+ versus CD8+ T cells." *PLoS Genetics.* 2017.](https://doi.org/10.1371/journal.pgen.1006643)
* Monocytes [Fairfax BP, *et al*. "Innate immune activity conditions the effect of regulatory variants upon monocyte gene expression." *Science*. 2014.](https://doi.org/10.1126/science.1246949)
* B cells [Fairfax BP, *et al*. "Genetics of gene expression in primary immune cells identifies cell type-specific master regulators and roles of HLA alleles." *Nature Genetics*. 2012.](https://doi.org/10.1038/ng.2205)
* NK cells [Gilchrist JJ, *et al*. "Natural Killer cells demonstrate distinct eQTL and transcriptome-wide disease associations, highlighting their role in autoimmunity." *Nature Communications.* 2022.](https://doi.org/10.1038/s41467-022-31626-4)
* Neutrophils [Naranbhai V, *et al*. "Genomic modulators of gene expression in human neutrophils." *Nature Communications*. 2015.](https://doi.org/10.1038/ncomms8545)

Colocalisation analysis also makes use of GWAS summary statistics from UK Biobank samples (http://www.nealelab.is/uk-biobank/) and immune-mediated diseases from the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics).

Figure 3 and Supplementary Figures 3 and 4 make use of publically-available gene expression and eQTL mapping data in leprosy patients and healthy controls:
* [Tio-Coma M, *et al*. "Blood RNA signature RISK4LEP predicts leprosy years before clinical onset." *EBioMedicine.* 2021.](https://doi.org/10.1016/j.ebiom.2021.103379)
* [Montoya DJ, *et al*. "Dual RNA-Seq of Human Leprosy Lesions Identifies Bacterial Determinants Linked to Host Immune Response." *Cell Reports.* 2019.](https://doi.org/10.1016/j.celrep.2019.02.109)
* [Manry J, *et al*. "Deciphering the genetic control of gene expression following Mycobacterium leprae antigen stimulation." *PLoS Genetics.* 2017.](https://doi.org/10.1371/journal.pgen.1006952)
* [Belone AdFF, *et al*. "Genome-Wide Screening of mRNA Expression in Leprosy Patients." *Frontiers in Genetics.* 2015.](https://doi.org/10.3389/fgene.2015.00334). 

