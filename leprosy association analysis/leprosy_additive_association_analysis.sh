#!/bin/bash

#Additive association analysis of leprosy risk at leprosy-associated loci in this study (Malawi)

./snptest_v2.5.1 \
-data malawi.gwas_snps.gen.gz \
malawi.sample \
-o malawi.gwas_snps.snptest.gz \
-frequentist add \
-method newml \
-pheno cc \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 \
-log malawi.gwas_snps.txt

#Additive association analysis of leprosy risk at leprosy-associated loci in this study (Mali)

./snptest_v2.5.1 \
-data mali.gwas_snps.gen.gz \
mali.sample \
-o mali.gwas_snps.snptest.gz \
-frequentist add \
-method newml \
-pheno cc \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 \
-log mali.gwas_snps.txt

#Additive association analysis of leprosy risk at loci associated with leprosy in previous studies (Malawi)

./snptest_v2.5.1 \
-data malawi.prior_snps.gen.gz \
malawi.sample \
-o malawi.prior_snps.snptest.gz \
-frequentist add \
-method newml \
-pheno cc \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 \
-log malawi.prior_snps.txt

#Additive association analysis of leprosy risk at loci associated with leprosy in previous studies (Mali)

./snptest_v2.5.1 \
-data mali.prior_snps.gen.gz \
mali.sample \
-o mali.prior_snps.snptest.gz \
-frequentist add \
-method newml \
-pheno cc \
-cov_names PC1 PC2 PC3 PC4 PC5 PC6 \
-log mali.prior_snps.txt

#Meta-analysis of leprosy risk at leprosy-associated loci in this study (Malawi & Mali)

./bingwa_v2.0-dev-linux-x86_64/bingwa \
-data malawi.gwas_snps.snptest.gz \
mali.gwas_snps.snptest.gz \
-flat-file \
-flip-alleles \
-o gwas_snps.meta.txt

#Meta-analysis of leprosy risk at loci associated with leprosy in previous studies (Malawi & Mali)

./bingwa_v2.0-dev-linux-x86_64/bingwa \
-data malawi.prior_snps.snptest.gz \
mali.prior_snps.snptest.gz \
-flat-file \
-flip-alleles \
-o prior_snps.meta.txt
