# 2018-10-20_Graphs_and_networks.R
#
# Purpose:  R code accompanying the FND-MAT-Graphs_and_networks short report assignment (Option B - Centrality)
#
# Version:  2.0
#
# Date:     2017  10  20
# Author:   Gabriela Morgenshtern (g.morgenshtern@mail.utoronto.ca)
#
# Versions:
#           2.1   Final (functional) version. Ready for presenting as deliverable
#           2.0   Second draft (functional) version. Including testing comments
#           1.0    First draft version, mostly notes from learning unit


##################################################################
###### Installing script requirements ############################
##################################################################

if (! require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}

if (!require(RColorBrewer, quietly=TRUE)) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

##################################################################
###### B.1 Informative overview plot ############################
##################################################################

# 1) get list of genes:
allGenes <- unique(c(scCCnet$protein1, scCCnet$protein2))

# 2) create adjacency matrix:  1 if tuple of prot1/prot2 exists in scCCnet, 0 otw
# Empty adjacency matrix (AM)
N <- length(allGenes)
AM <- matrix(numeric(N * N), ncol = N)
rownames(AM) <- allGenes
colnames(AM) <- allGenes

# Fill AM values based on existing rows in scCCnet
for (iRow in 1:(length(scCCnet$protein1))) {
  protein1 <- as.character(scCCnet[iRow,]["protein1"])
  protein2 <- as.character(scCCnet[iRow,]["protein2"])
  AM[protein1, protein2] <- 1
}

# 3) create plot of adjacency matrix: myG
myG <- graph_from_adjacency_matrix(AM, mode = "undirected")

# 4) calculate layout coordinates:
myGxy <- layout_with_fr(myG)

# 5) Plot the network; color and scale the nodes according to their centrality score:
bC <- centr_betw(myG)
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw + 1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = myGxy*1.0,
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

# Create a legend and title
title("Gene interaction network based on the scCCnet dataset", line = -1)

Group <- gl(10, 1, 10, labels = c(seq_len(10)))
legend("bottomright",bty = "n",
       legend=levels(Group),
       fill=brewer.pal(max(nodeBetw), "RdYlGn"),
       border=NA,
       title="Centrality scale")

# Shortest path and diameter information:
diameter(myG)
get_diameter(myG)
lines(myGxy[get_diameter(myG),]*1.0, lwd=8, col="#ff63a788")

##################################################################
###### B.2: Betweeness centrality to scale and colour nodes ######
##################################################################

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
