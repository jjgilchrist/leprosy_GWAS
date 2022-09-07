# leprosy_GWAS

**Genome-wide association study of leprosy in Malawi and Mali.** (PLoS Pathogens, 2022)

https://www.medrxiv.org/content/10.1101/2022.01.31.22270046v1

Code and source data for GWAS of leprosy in Malawi and Mali.

**Dataset**
612 Malawian participants: 284 leprosy cases, 328 controls.
350 Malian participants: 208 leprosy cases, 142 controls.

**Abstract**
Leprosy is a chronic infection of the skin and peripheral nerves caused by *Mycobacterium leprae*. Despite recent improvements in disease control, leprosy remains an important cause of infectious disability globally. Large-scale genetic association studies in Chinese, Vietnamese and Indian populations have identified over 30 susceptibility loci for leprosy. There is a significant burden of leprosy in Africa, however it is uncertain whether the findings of published genetic association studies are generalizable to African populations. To address this, we conducted a genome-wide association study (GWAS) of leprosy in Malawian and Malian individuals. In that analysis, we replicated four risk loci previously reported in China, Vietnam and India; MHC Class I and II, *LACC1* and *SLC29A3*. We further identified a novel leprosy susceptibility locus at 10q24 (rs2015583; combined p=8.81x10<sup>-9</sup>; OR=0.51 [95% CI 0.40-0.64]). Using publicly-available data we characterise regulatory activity at this locus, identifying *ACTR1A* as a candidate mediator of leprosy risk. This locus shows evidence of recent positive selection and demonstrates pleiotropy with established risk loci for inflammatory bowel disease and childhood-onset asthma. A shared genetic architecture for leprosy and inflammatory bowel disease has been previously described. We expand on this, strengthening the hypothesis that selection pressure driven by leprosy has shaped the evolution of autoimmune and atopic disease in modern populations. More broadly, our data highlights the importance of defining the genetic architecture of disease across genetically diverse populations, and that disease insights derived from GWAS in one population may not translate to all affected populations.

**Overview of repository**
* Figure folders: contains R script and source data to reproduce each main and supplementary Figure from the manuscript.
* Leprosy association analysis: contains scripts and sample genotype, phenotype and covariate data to perform association analysis at leprosy-associated SNPs (p<1x10<sup>-5</sup>, n=182) identified in the discovery GWAS (Malawi) and meta-analysis, and leprosy-associated SNPs (n=25) identified in previous studies. Association analysis was performed using SNPTEST v2.5.1: [Marchini J, *et al*. "A new multipoint method
for genome-wide association studies by imputation of genotypes." *Nature Genetics.* 2009.](https://doi.org/10.1038/ng2088). Meta-analysis was performed using bingwa v2.0: [Malaria Genomic Epidemiology Network. "A novel locus of resistance to severe malaria in a region of ancient balancing selection." *Nature.* 2015.](https://doi.org/10.1038/nature15390)

**Data availability and sources**


GWAS summary statistics are available through the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics; accession codes: Malawi, GCST90129399; Mali, GCST90129400; meta-analysis, GCST90129401). Patient level genotype and phenotype data for additional MalariaGEN control samples (Mali) are available are available via the European Genome-Phenome Archive, with accession code EGAS00001001311.


