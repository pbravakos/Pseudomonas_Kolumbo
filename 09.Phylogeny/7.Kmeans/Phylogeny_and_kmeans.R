library("ape")

MLTree <- read.tree("Pseudomonas_RaxML.raxml.bestTree")


#The matrix edge contains the beginning and ending node number for all the nodes and tips in the tree. By convention, the tips of the tree are numbered 1 through n for n tips; and the nodes are numbered n + 1 through n + m for m nodes. m = n - 1 for a fully bifurcating tree. This is just to keep track of which nodes are internal and which are leaves.
MLTree.edge <- MLTree$edge


#The vector tip.label contains the labels for all the tips in the tree. The order of tip.label is the order of the tips numbered 1 through n in edge.
MLTree$tip.label


#The integer Nnode contains the number of internal nodes in the tree, including the root of the tree if the tree is rooted.
MLTree$Nnode

plot(MLTree, cex = 0.4)
nodelabels()
tiplabels()

#"ape" and most other phylogenetics packages are equipped to handle phylogenies that are binary or multifurcating; however not all functions will be. We can easily check if our tree is binary, and convert between binary & multifurcating trees.
is.binary.tree(MLTree)


## Calculate patristic distances in a matrix where the values in the matrix are the sum of the branch lengths separating each pair of species. Thus, if two species are far apart on the phylogeny their distance will be larger than that for two closely related species.

p.dist.mat_ML <- cophenetic(MLTree)

## To calculate a patristric distance between two OTUs

MLTree$tip.label
p.dist.mat_ML["Pseudomonas_taetrolens_BS3652", "Pseudomonas_agarici_NCPPB_2472"]
p.dist.mat_ML["Pseudomonas_taetrolens_BS3652", "Pseudomonas_frederiksbergensis_ERGS4"]
p.dist.mat_ML["PSR1JGgm14", "IKP2SMP32"]
max(p.dist.mat_ML)
min(p.dist.mat_ML) #it is zero because it is the distance whith itself!!
summary(p.dist.mat_ML)

# Kmeans analysis
library("cluster")
kmean<-kmeans(p.dist.mat_ML, centers=16, iter.max=1000)
clusplot(p.dist.mat_ML, kmean$cluster, main = "Kmeans clustering", color=FALSE, shade=FALSE, labels=2, lines=3, plotchar = FALSE, span = TRUE, col.p = "black", col.txt = "blue", col.clus ="red", cex.txt = 0.05, cex = 0.09)

