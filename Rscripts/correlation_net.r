args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(qgraph))
suppressMessages(library(vegan))
suppressMessages(library(MCL))
suppressMessages(library(SpiecEasi))

read.table(args[1], header =T, row.names =1, sep = "\t") -> inFile
t(inFile) -> tFile
tFile * 1000 -> aboveZero
round(aboveZero) -> forNet
args[2] -> outName



##### Function 1: Construct microbiome network using permutation
build.cor.net <- function(MB, method, num_perms, sig_level) {
  taxa <- dim(MB)[2]
  MB.mat <- array(0, dim = c(taxa, taxa, num_perms + 1))
  # Perform permutation
  MBperm <- permatswap(MB, "quasiswap", times = num_perms)
  # Convert to relative abundance
  MB.relative <- MB / rowSums(MB)
  MB.mat[,,1] <- as.matrix(cor(MB.relative, method = method))
  for(p in 2:num_perms) {
    MBperm.relative <- MBperm$perm[[p-1]] / rowSums(MBperm$perm[[p-1]])
    MB.mat[, , p] <- as.matrix(cor(MBperm.relative, method = method))
  }
  # Get p-values
  pvals <- sapply(1:taxa,
                  function(i) sapply(1:taxa, function(j)
                    sum(MB.mat[i, j, 1] > MB.mat[i, j, 2:num_perms])))
  pvals <- pvals / num_perms
  # p-value correction
  pvals_BH <- array(p.adjust(pvals, method = "BH"),
                    dim=c(nrow(pvals), ncol(pvals)))
  # Build adjacency matrix
  adj.mat <- ifelse(pvals_BH >= (1 - sig_level), 1, 0)
  # Add names to rows & cols
  rownames(adj.mat) <- colnames(MB)
  colnames(adj.mat) <- colnames(MB)
  # Build and return the network
  graph <- graph.adjacency(adj.mat, mode = "undirected", diag = FALSE)
}

# Execute this command after running Function 1
cor.net.2 <- build.cor.net(forNet,
                           method = 'pearson',
                           num_perms = 100,
                           sig_level = 0.01)
#To find hubs (large number of links compared to other nodes) - can be thought of as keystone species

# Use cor.net.2 for the rest of the method
#cor.net.2
# Hub detection
cor.net.2.cn <- closeness(cor.net.2)
cor.net.2.bn <- betweenness(cor.net.2)
cor.net.2.pr <- page_rank(cor.net.2)$vector
cor.net.2.hs <- hub_score(cor.net.2)$vector

#The net.hs is a vector of OTUs with value from 0-1, with values closer to 1 being the more likley hubs

# Sort the species based on hubbiness score
cor.net.2.hs.sort <- sort(cor.net.2.hs, decreasing = TRUE)
# Choose the top 5 keystone species
cor.net.2.hs.top5 <- head(cor.net.2.hs.sort, n = 5)

#Make clusters - i don't really understand the details here

# Get clusters
cor.2.wt <- walktrap.community(cor.net.2)
cor.2.ml <- multilevel.community(cor.net.2)
# Get membership of walktrap clusters
#membership(wt)
# Get clusters using MCL method
cor.2.adj <- as_adjacency_matrix(cor.net.2)
cor.2.mc <- mcl(cor.2.adj, addLoops = TRUE)

#Basic features of the network

# Network features
cor.2.nodes <- V(cor.net.2)
cor.2.edges <- V(cor.net.2)
cor.2.node.names <- V(cor.net.2)$name
cor.2.num.nodes <- vcount(cor.net.2)
cor.2.num.edges <- ecount(cor.net.2)


# Function 2: Plot network with node size scaled to hubbiness
plot.net <- function(net, scores, outfile, title) {
  # Convert node label from names to numerical IDs.
  features <- V(net)$name
  col_ids <- seq(1, length(features))
  V(net)$name <- col_ids
  node.names <- features[V(net)]
  # Nodes' color.
  V(net)$color <- "white"
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex = 1)
  title(title, cex.main = 4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  #legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}

plot.net(cor.net.2, cor.net.2.hs, outfile = paste0("/Users/dave/Desktop/microbiom_R_network/", outName, "_cor.2_uncolored"), title = outName)

cor.2.AP <- articulation.points(cor.net.2)

outName_degree.dist <- degree_distribution(cor.net.2, mode = "all", cumulative = F)
as.data.frame(outName_degree.dist) -> outDF
write.table(outDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/", outName, "_cor.2_degree.txt"), quote =F,  col.names = F, row.names = T)

# Function 3: Plot network with clusters and node size scaled to hubbiness
plot.net.cls <- function(net, scores, cls, AP, outfile, title) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  colors <- sample(colours(), length(cls))
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colors[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  col_ids <- seq(1, length(node.names))
  V(net)$name <- col_ids
  # To draw a halo around articulation points.
  #AP <- lapply(names(AP), function(x) x)
  # marks <- lapply(1:length(AP), function(x) which(node.names == AP[[x]]))
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Customized layout to avoid nodes overlapping.
  e <- get.edgelist(net)
  class(e) <- "numeric"
  l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
                                         area=8*(vcount(net)^2),
                                         repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex=0.9,
       vertex.label.color = "black",
       mark.border="black",
       #mark.groups = marks,
       mark.col = "white",
       mark.expand = 10,
       mark.shape = 1,
       layout=l)
  title(title, cex.main=4)
  # Plot legend containing OTU names.
  #labels = paste(as.character(V(net)), node.names, sep = ") ")
  #legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}
# Execute this command after running Function 3
plot.net.cls(cor.net.2, cor.net.2.hs, cor.2.wt,
             outfile =  paste0("/Users/dave/Desktop/microbiom_R_network/", outName, "_cor.2_colored"), title = outName)
