#!/bin/bash

#Additive association analysis of leprosy risk at the ACTR1A locus and SNPs previously associated with leprosy in other populations (Malawi)

./snptest_v2.5.6 \
-data malawi.bgen \
malawi.sample \
-o malawi.snptest.tsv \
-frequentist add \
-method newml \
-pheno cc \
-log malawi.log

#Additive association analysis of leprosy risk at the ACTR1A locus and SNPs previously associated with leprosy in other populations (Mali)

./snptest_v2.5.6 \
-data mali.bgen \
mali.sample \
-o mali.snptest.tsv \
-frequentist add \
-method newml \
-pheno cc \
-log mali.log

#Meta-analysis of additive leprosy risk (Malawi & Mali)

./bingwa_v2.0.6-osx/bingwa \
-data malawi.snptest.tsv \
mali.snptest.tsv \
-flip-alleles \
-o meta.sqlite

#inspect sqlite output in R

library(RSQLite)

db = dbConnect( dbDriver( "SQLite" ), "meta.sqlite" )
D = dbGetQuery( db, "SELECT * FROM BingwaMetaView" )
D[order(D[,19]),]

#association analysis excludes additional MalariaGEN control samples and principal components as covariate data - effect estimates are therefore subtly different to those in the ful analysis

