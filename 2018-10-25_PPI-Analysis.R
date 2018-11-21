# 2018-10-25_PPI-Analysis.R
#
# Purpose:  R code for the submission of the PPI Analysis task
#
# Version:  1.0
#
# Date:     2017  10  25
# Author:   Gabriela Morgenshtern (g.morgenshtern@mail.utoronto.ca)
#
# Versions:
#
#           1.0    First draft version, mostly notes from learning unit

if (!require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}

if (!require(biomaRt, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("biomaRt")
  library(biomaRt)
}

## Script setup
load("./data/STRINGedges.RData")

head(STRINGedges)

# =    4  Task for submission  =================================================


# Write a loop that will go through your personalized list of Ensemble IDs and
#    for each ID:
#    --  print the ID,
#    --  print the first row's HGNC symbol,
#    --  print the first row's wikigene description.
#    --  print the first row's phenotype.
#
# (Hint, you can structure your loop in the same way as the loop that
# created CPdefs. )

# Place the R code for this loop and its output into your report if you are
# submitting a report for this unit. Please read the requirements carefully.

(ENSPsel <- names(V(gSTR)[BCsel]))

# creating the personalized ENSPsel list:
set.seed(1000703679)
(ENSPsel <- sample(ENSPsel))

# specifying the database:
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

CPdefs <- list()
for (ID in ENSPsel) {
  CPdefs[[ID]] <- getBM(filters = "ensembl_peptide_id",
                        attributes = c("hgnc_symbol",
                                       "wikigene_description",
                                       "interpro_description",
                                       "phenotype_description"),
                        values = ID,
                        mart = myMart)
}

charSpacer <- rep('*', 40)
for (ID in ENSPsel) {

  cat("ENSP ID:", ID,"\n")
  cat("hgnc_symbol:", CPdefs[[ID]][1,][1,"hgnc_symbol"], "\n")
  cat("wikigene_description:", CPdefs[[ID]][1,][1,"wikigene_description"], "\n")
  cat("interpro_description:", CPdefs[[ID]][1,][1,"interpro_description"], "\n")
  cat("phenotype_description:", CPdefs[[ID]][1,][1,"phenotype_description"], "\n")
  cat(charSpacer, "\n")
}

# [END]
