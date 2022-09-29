#!/bin/bash

#Multinomial association analysis of leprosy risk (stratified by leprosy subtype) at the ACTR1A locus and SNPs previously associated with leprosy in other populations (Malawi)

./snptest_v2.5.6 \
-data malawi.bgen \
malawi.sample \
-o malawi.mnom.snptest.tsv \
-frequentist add \
-method newml \
-pheno pbmb \
-baseline_phenotype control \
-log malawi.mnom.log

#Multinomial association analysis of leprosy risk (stratified by leprosy subtype) at the ACTR1A locus and SNPs previously associated with leprosy in other populations (Mali)

./snptest_v2.5.6 \
-data mali.bgen \
mali.sample \
-o mali.mnom.snptest.tsv \
-frequentist add \
-method newml \
-pheno pbmb \
-baseline_phenotype control \
-log mali.mnom.log

#Meta-analysis of multinomial leprosy risk (Malawi & Mali)

./bingwa_v2.0.6-osx/bingwa \
-data malawi.mnom.snptest.tsv \
mali.mnom.snptest.tsv \
-flip-alleles \
-o meta.mnom.sqlite

#inspect sqlite output in R

library(RSQLite)

db = dbConnect( dbDriver( "SQLite" ), "meta.mnom.sqlite" )
D = dbGetQuery( db, "SELECT * FROM BingwaMetaView" )
D[order(D[,23]),]

#association analysis excludes additional MalariaGEN control samples and principal components as covariate data - effect estimates are therefore subtly different to those in the ful analysis

