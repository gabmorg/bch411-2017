# Mutation_Impact_Oral_Exam
#
# Purpose:  A Bioinformatics Course:
#             R code accompanying the oral exam preparations for the Mutation Impact unit
#
# Version:  1.1
#
# Date:     2018  11 19  - 2018  11 21
# Author:   Gabriela Morgenshtern (g.morgenshtern@mail.utoronto.ca)
#
# Versions:
#           1.1    loop that keeps track of frequency of mutation types only
#           1.0    loop that keeps track of the mutation locations

# Testing for data availability and loading necessary packages:
expect_true(file.exists("./data/ABC-INT-Mutation_impact.RData"))
load(file = "./data/ABC-INT-Mutation_impact.RData")

# Load genetic code tables from the Biostrings package
if (! require(Biostrings, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biostrings")
  library(Biostrings)
}



# Overall Method: adopting the approach used in the Genetic Code Optimality unit

######## Helper function definitions ########
# To mutate, we split a random codon into it's three nucleotides,
# then also randomly replace one of the three with another nucleotide.
randMutSingle <- function(vC) {
  # Parameter:
  #    vC   chr     a vector of codons
  # Value:  vector  a 2 element vector: a chr codons with a
  #                 single point mutation from vC, and the num index of mutation in vC

  nuc <- c("A", "C", "G", "T")
  imutCodon <- sample(seq_along(vC), 1)      # select a random index from vC (codon vector)
  triplet <- unlist(strsplit(vC[imutCodon], ""))         # split into three nucl.
  iNuc <- sample(1:3, 1)                         # choose one of the three
  mutNuc <- sample(nuc[nuc != triplet[iNuc]], 1) # choose a mutated nucleotide
  triplet[iNuc] <- mutNuc                        # replace the original
  vC[imutCodon] <- paste0(triplet, collapse = "")        # collapse it to a codon

  retVect <- c(vC[imutCodon], imutCodon)
  return(retVect)
}

# To translate, we pass a single codon and the genetic code. In the future,
# this will allow us to avoid nested loops and decrease runtime
traForSingle <- function(codon, GC) {
  # Parameters:
  #      codon   chr   a single codon
  #      GC      chr   a genetic code
  # Value:
  #      A vector of a single amino acid
  vAA <- character(length(codon))
  vAA <- GC[codon]         # translate and store
  return(vAA)

}

traFor <- function(vC, GC) {
  # Parameters:
  #      vC   chr   a codon vector
  #      GC   chr   a genetic code
  # Value:
  #      A vector of amino acids
  vAA <- character(length(vC))

  for (i in seq_along(vC)) {
    vAA[i] <- GC[vC[i]]         # translate and store
  }
  return(vAA)

}

evalMutSingle <- function(nat, mut, iMut){
  # Parameters:
  #      nat      chr   a vector of the original AA seq
  #      mut      chr   the mutated AA
  #      iMut     num   the index of the mutated AA in the nat
  # Value:
  #      chr      "sil" = "silent", "mis"="missense", "non"= "nonsense"
  #               for the type of mutation observed

  mutType <- character(3)

  if(mut == "*"){
    mutType <- "non"
  }
  else if(nat[iMut] == mut){
    mutType <- "sil"
  }
  else if(nat[iMut] != mut && nat[iMut] != "*"){
    mutType <- "mis"
  }

  return(mutType)
}
#####

######## Stochastic Environment Simulation - Experimental Setup ########
# We keep track of the index where the mutation occured to allow comparative
# analysis of the mutation locations, which in stochastic simulation are expected to
# be random, whereas in IntOGen are expected to be specific

# 1. Drop the stop codons:
KRascodonsExp <- KRascodons[-length(KRascodons)]
OR1A1codonsExp <- OR1A1codons[-length(OR1A1codons)]
PTPN11codonsExp <- PTPN11codons[-length(PTPN11codons)]

# 2. Get the original AA sequences of each gene:
KRAS_AA <- traFor(KRascodonsExp, GENETIC_CODE)
OR1A1_AA <- traFor(OR1A1codonsExp, GENETIC_CODE)
PTPN11_AA <- traFor(PTPN11codonsExp, GENETIC_CODE)

# 3. GLOBAL counters of the observed mutations:
#KRAS mutation counts:

KRASmutations <- numeric(3)
names(KRASmutations) <- c("silent", "missense", "nonsense")

#OR1A1 mutation counts:
OR1A1mutations <- numeric(3)
names(OR1A1mutations) <- c("silent", "missense", "nonsense")

#PTPN11 mutation counts:
PTPN11mutations <- numeric(3)
names(PTPN11mutations) <- c("silent", "missense", "nonsense")

######## Experimentation ########
# 1. Create random point mutations N times:
# Loop plan:
# copy each gene's sequence into new var each loop
# randomly mutate a codon in each of the genes
# compare mutated codon to the original codon, classify and add count to appropriate category
# keep count of each gene's silent (same AA), missense (new AA), and nonsense (stop codon)

# Create a point mutation in each of the 3 genes (KRas, PTPN11 and OR1A1) N times
# Keep track of the number of missense, silent ("synonymous"),
# and nonsense ("truncating")" mutations you find.

set.seed(1000703679)

N <- 100000
countMuts <- 0

for (i in 1:N) {
  # mutate all 3 genes and store the mutation index
  mutKRas <- randMutSingle(KRascodonsExp)
  imutKRas <- as.numeric(mutKRas[2])
  mutOR1A1 <- randMutSingle(OR1A1codonsExp)
  imutOR1A1 <- as.numeric(mutOR1A1[2])
  mutPTPN11 <- randMutSingle(PTPN11codonsExp)
  imutPTPN11 <- as.numeric(mutPTPN11[2])

  # translate all 3 genes
  mutKRas <- traForSingle(mutKRas[1], GENETIC_CODE)
  mutOR1A1 <- traForSingle(mutOR1A1[1], GENETIC_CODE)
  mutPTPN11 <- traForSingle(mutPTPN11[1], GENETIC_CODE)

  # evaluate all 3 genes
  mutTypeKRas <- evalMutSingle(KRAS_AA, mutKRas, imutKRas)
  mutTypeOR1A1 <- evalMutSingle(OR1A1_AA, mutOR1A1, imutOR1A1)
  mutTypePTPN11 <- evalMutSingle(PTPN11_AA, mutPTPN11, imutPTPN11)

  # update mutation counters (repetitive conditionals kept for readers' clarity)
  # KRas
  if(mutTypeKRas == "non"){
    KRASmutations[3] <- KRASmutations[3] + 1
  }
  else if(mutTypeKRas == "mis"){
    KRASmutations[2] <- KRASmutations[2] + 1
  }
  else if(mutTypeKRas == "sil"){
    KRASmutations[1] <- KRASmutations[1] + 1
  }

  # OR1A1
  if(mutTypeOR1A1 == "non"){
    OR1A1mutations[3] <- OR1A1mutations[3] + 1
  }
  else if(mutTypeOR1A1 == "mis"){
    OR1A1mutations[2] <- OR1A1mutations[2] + 1
  }
  else if(mutTypeOR1A1 == "sil"){
    OR1A1mutations[1] <- OR1A1mutations[1] + 1
  }

  # PTPN11
  if(mutTypePTPN11 == "non"){
    PTPN11mutations[3] <- PTPN11mutations[3] + 1
  }
  else if(mutTypePTPN11 == "mis"){
    PTPN11mutations[2] <- PTPN11mutations[2] + 1
  }
  else if(mutTypePTPN11 == "sil"){
    PTPN11mutations[1] <- PTPN11mutations[1] + 1
  }

  # sanity check
  countMuts <- countMuts + 1
}


######## Analysis ########
# Contrast that with the relative frequency of the mutations in each
# category reported on the IntOGen Web page for each of the three genes.
#KRas
KRASSilentPerC <- KRASmutations[1]/N *100
KRASMissensePerC <- KRASmutations[2]/N *100
KRASNonsensePerC <- KRASmutations[3]/N *100

#OR1A1
OR1A1SilentPerC <- OR1A1mutations[1]/N *100
OR1A1MissensePerC <- OR1A1mutations[2]/N *100
OR1A1NonsensePerC <- OR1A1mutations[3]/N *100

#PTPN11
PTPN11SilentPerC <- PTPN11mutations[1]/N *100
PTPN11MissensePerC <- PTPN11mutations[2]/N *100
PTPN11NonsensePerC <- PTPN11mutations[3]/N *100

# Print relative frequencies from stochastic simulation:
cat("************ KRas ************\n")
cat("KRas stochastic environment simulated:", KRASSilentPerC,
"% silent mutations\n",
    "KRas stochastic environment simulated:", KRASMissensePerC,
    "% missense mutations\n",
    "KRas stochastic environment simulated:", KRASNonsensePerC,
    "% nonsense mutations\n")

cat("************ OR1A1 ************\n")
cat("OR1A1 stochastic environment simulated:", OR1A1SilentPerC,
    "% silent mutations\n",
    "OR1A1 stochastic environment simulated:", OR1A1MissensePerC,
    "% missense mutations\n",
    "OR1A1 stochastic environment simulated:", OR1A1NonsensePerC,
    "% nonsense mutations\n")

cat("************ PTPN11 ************\n")
cat("PTPN11 stochastic environment simulated:", PTPN11SilentPerC,
    "% silent mutations\n",
    "PTPN11 stochastic environment simulated:", PTPN11MissensePerC,
    "% missense mutations\n",
    "PTPN11 stochastic environment simulated:", PTPN11NonsensePerC,
    "% nonsense mutations\n")

######## Output ########
# ************ KRas ************
# KRas stochastic environment simulated: 20.682 % silent mutations
# KRas stochastic environment simulated: 74.197 % missense mutations
# KRas stochastic environment simulated: 5.121 % nonsense mutations
# ************ OR1A1 ************
# OR1A1 stochastic environment simulated: 23.406 % silent mutations
# OR1A1 stochastic environment simulated: 73.21 % missense mutations
# OR1A1 stochastic environment simulated: 3.384 % nonsense mutations
# ************ PTPN11 ************
# PTPN11 stochastic environment simulated: 21.044 % silent mutations
# PTPN11 stochastic environment simulated: 74.259 % missense mutations
# PTPN11 stochastic environment simulated: 4.697 % nonsense mutations
###### FROM UNIT NOTES ######
# Passenger mutations are expected to be randomly
# distributed throughout the genome, driver mutations
# are expected to have either a gain of function or loss
# of function effect.

# Gain of function mutations: very specific, targeting only a small
# number of amino acids in a defined region of the protein.
# We actually expect purifying selection against mutations
# elsewhere.
# Loss of function mutations: should be enriched in
# missense and nonsense mutations relative to
# silent mutations.

# The task of this unit is to analyze the relative
# frequencies of neutral, missense and nonsense
# mutations in a gene, and contrast that with the
# frequencies one would expect if the distribution
# of mutations was purely due to chance.
##############################



# In the loop experiment above, how many of the mutations were of each category for each gene?
# KRas:
# missense: 20.7%
# silent ("synonymous"): 74.2%
# nonsense ("truncating"): 5.1%
#
# PTPN11:
# missense: 21%
# silent ("synonymous"): 74.3%
# nonsense ("truncating"): 4.7%
#
# OR1A1:
# missense: 23.4%
# silent ("synonymous"): 73.2%
# nonsense ("truncating"): 3.4%
#
# In Intogen, what is the frequency of mutations of each category for each gene?
# KRAS:
# https://www.intogen.org/search?gene=KRAS
# missense: 99%
# silent ("synonymous"): 1%
# nonsense ("truncating"): 0%
#
# PTPN11:
# https://www.intogen.org/search?gene=PTPN11#frequency
# missense: 82%
# silent ("synonymous"): 15%
# nonsense ("truncating"): 4%
#
# OR1A1:
# https://www.intogen.org/search?gene=OR1A1
# missense: 59%
# silent ("synonymous"): 35%
# nonsense ("truncating"): 6%

# (FOR SUBMISSION, MAKE ABOVE INTO COMPARISON TABLE)

# Describe whether you think there is an important difference between
# the expected categories of mutations (i.e. the stochastic background
# that you simulated), and categories of mutations that were observed in
# cancer genomes.

###### FROM UNIT NOTES ######
# system can go terribly wrong if Ras gets mutated in a
# way that damages its catalytic activity and prevent
# GTP hydrolysis
#
# An interesting new development therefore was the
# recent discovery that a phosphatase - PTPN11 -
#   somehow works synergistically with Ras to
# facilitate its activation of effectors: inhibition
# of PTPN11 suppressed oncogenesis[2]. If this is a
# pathophysiologically relevant effect, we expect
# cancer mutations to spare PTPN11, or even to
# deregulate it to enhance its activity. Do they?
##############################

###### FROM UNIT NOTES ######
# evaluate mutations of the KRas gene, a known cancer
# driver, an olfactory receptor (OR1A1), most likely
# not involved in cancer, and the PTPN11 phosphatase,
# a gene of interest whose role in cancer we would like
# to understand better.
#############################

## Difference may lie in the fact that the stochastic background simulated
## has no underlying motivation for creating a specific mutation, whereas the
## recorded cancer mutations have to do with cancer survival, too

# Document your activities and results in your Journal. Add a brief
# conclusion / interpretation.
## were the findings signifcantly different for each category?
