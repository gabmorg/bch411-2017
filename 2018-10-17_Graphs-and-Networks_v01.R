# FND-MAT-Graphs_and_networks.R
#
# Purpose:  R code accompanying the FND-MAT-Graphs_and_networks short report assignment (Option B - Centrality)
#
# Version:  2.0
#
# Date:     2017  10  20
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.0   Second final (functional) version. Including testing comments
#           1.0    First final version for learning units.
#           0.1    First code copied from 2016 material.


### Installing script requirements
if (! require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}

if (!require(RColorBrewer, quietly=TRUE)) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

### Setting variables
sccdata <- scCCnet

# Steps #
## Informative overview plot setup:
# 1) get list of genes:
allGenes <- unique(c(sccdata$protein1, sccdata$protein2))

# 2) create adjacency matrix:  1 if tuple of prot1/prot2 exists in sccdata, 0 otw
# Empty adjacency matrix (AM)
N <- length(allGenes)
AM <- matrix(numeric(N * N), ncol = N)
rownames(AM) <- allGenes
colnames(AM) <- allGenes

# Fill AM values based on existing rows in sccdata
for (iRow in 1:(length(sccdata$protein1))) {
  protein1 <- as.character(sccdata[iRow,]["protein1"])
  protein2 <- as.character(sccdata[iRow,]["protein2"])
  AM[protein1, protein2] <- 1

  # TESTING ITERATIONS
  # print(iRow)
  # print(protein1)
  # print(protein2)
  # print(AM[protein1, protein2])
}

# 3) create plot of adjacency matrix: myG
myG <- graph_from_adjacency_matrix(AM, mode = "undirected")

# 4) calculate layout coordinates: OR myGxy <- layout_with_graphopt(myG)
myGxy <- layout_with_fr(myG)
# myGxy <- layout_with_graphopt(myG, niter = 40, charge = 0.0012, spring.length = 0, mass = 20)
# myGxy <- layout_with_lgl(myG)

## Betweeness centrality to scale and colour nodes:
# 5) plot! using betweeness centrality as in class example:
bC <- centr_betw(myG)
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw + 1)) + 1


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = myGxy*1.05,
     rescale = FALSE,
     xlim = c(min(myGxy[,1]) * 0.99, max(myGxy[,1]) * 1.01),
     ylim = c(min(myGxy[,2]) * 0.99, max(myGxy[,2]) * 1.01),
     vertex.color=brewer.pal(max(nodeBetw), "RdYlGn")[nodeBetw],
     vertex.size = 15 + (5 * degree(myG)),
     vertex.label = V(myG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.4,
     label.dist = 1,
     label.degree = 0,
     edge.width = 0.5
     )
par(oPar)
title("Gene interaction network based on the scCCnet dataset", line = -2)
# Create a legend
Group <- gl(10, 1, 10, labels = c(seq_len(10)))
legend("topright",bty = "n",
       legend=levels(Group),
       fill=brewer.pal(max(nodeBetw), "RdYlGn"),
       border=NA,
       title="Centrality scale")

# Shortest path and diameter information:
diameter(myG)
get_diameter(myG)
lines(myGxy[get_diameter(myG),]*1.05, lwd=8, col="#ff63a788")

## Analyzing graph
# 6) order nodes according to their centrality score (code example needs playing with):
centralityScores <- sort(components(myG)$membership, decreasing = TRUE)

# 7) get top  10 nodes by centrlity score
topCentrality10 <- centralityScores[1:10]
topCentrality10

## Researching analysis findings --> See report writeup
# 8) Interpret the results. What are the highest-centrality nodes? Do genes with high centrality have anything in common?
# look up in: https://yeastmine.yeastgenome.org/yeastmine/bag.do

# Genes:          YLR330W YOR299W YLR084C YOR301W YJR072C YOR262W YLR243W YGR245C YLR074C YEL029C
# Centralities:   11      11      10      10       9       9       9       8       8       7


