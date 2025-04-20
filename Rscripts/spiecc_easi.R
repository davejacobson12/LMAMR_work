args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(qgraph))
suppressMessages(library(vegan))
suppressMessages(library(MCL))
suppressMessages(library(SpiecEasi))

read.table(args[1], header =T, row.names =1, sep = "\t") -> inFile
t(inFile) -> tFile
tFile * 1000 -> forNet
args[2] -> outName


SpiecEasi.matrix <- spiec.easi(forNet, method = "glasso", lambda.min.ratio = 1e-2, nlambda = 20, icov.select.params = list(rep.num = 50))
print(SpiecEasi.matrix$refit)
# Add OTU names to rows and columns
rownames(SpiecEasi.matrix$refit) <- colnames(forNet)
# Build network from adjacency
SpiecEasi.net <- graph.adjacency(SpiecEasi.matrix$refit,
                              mode = "undirected",
                              diag = FALSE)

#To find hubs (large number of links compared to other nodes) - can be thought of as keystone species

# Use SpiecEasi.net for the rest of the method
#SpiecEasi.net
# Hub detection
SpiecEasi.net.cn <- closeness(SpiecEasi.net)
SpiecEasi.net.bn <- betweenness(SpiecEasi.net)
SpiecEasi.net.pr <- page_rank(SpiecEasi.net)$vector
SpiecEasi.net.hs <- hub_score(SpiecEasi.net)$vector

#The net.hs is a vector of OTUs with value from 0-1, with values closer to 1 being the more likley hubs

# Sort the species based on hubbiness score
SpiecEasi.net.hs.sort <- sort(SpiecEasi.net.hs, decreasing = TRUE)
# Choose the top 5 keystone species
SpiecEasi.net.hs.top5 <- head(SpiecEasi.net.hs.sort, n = 5)

as.data.frame(SpiecEasi.net.hs.top5) -> hubs
hubs$pop <- rep(outName, 5)

write.table(hubs, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_5hubs.txt"), quote =F,  col.names = F, row.names = T)

#Make clusters - i don't really understand the details here

# Get clusters
SpiecEasi.wt <- walktrap.community(SpiecEasi.net)
SpiecEasi.ml <- multilevel.community(SpiecEasi.net)
# Get membership of walktrap clusters
membership(SpiecEasi.wt) -> member
print(member)
#print(max(member))
# Get clusters using MCL method
SpiecEasi.adj <- as_adjacency_matrix(SpiecEasi.net)
SpiecEasi.mc <- mcl(SpiecEasi.adj, addLoops = TRUE)

connectVec <- vector(mode = "numeric", length = length(V(SpiecEasi.net)$name))
nodeNames <- V(SpiecEasi.net)$name
for(i in 1:length(V(SpiecEasi.net)$name)){
    connectVec[i] <- length(neighbors(SpiecEasi.net, nodeNames[i]))
}
as.data.frame(connectVec) -> connectDF
connectDF$taxa <- nodeNames
zero <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec <= 0)))]
one5 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 1 & connectDF$connectVec <= 5)))]
six10 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 6 & connectDF$connectVec <= 10)))]
eleven15 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 11 & connectDF$connectVec <= 15)))]
sixteen20 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 16 & connectDF$connectVec <= 20)))]
twentyone25 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 21 & connectDF$connectVec <= 25)))]
twentysix30 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 26 & connectDF$connectVec <= 30)))]
thirtyone35 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 31 & connectDF$connectVec <= 35)))]
thirtysix40 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 36 & connectDF$connectVec <= 40)))]
fortyAbove <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 41)))]
as.data.frame(zero) -> zeroDF
as.data.frame(one5) -> one5DF
as.data.frame(six10) -> six10DF
as.data.frame(eleven15) -> eleven15DF
as.data.frame(sixteen20) -> sixteen20DF
as.data.frame(twentyone25) -> twentyone25DF
as.data.frame(twentysix30) -> twentysix30DF
as.data.frame(thirtyone35) -> thirtyone35DF
as.data.frame(thirtysix40) -> thirtysix40DF
as.data.frame(fortyAbove) -> fortyAboveDF
zeroDF$pop <- rep(outName, length(rownames(zeroDF)))
one5DF$pop <- rep(outName, length(rownames(one5DF)))
six10DF$pop <- rep(outName, length(rownames(six10DF)))
eleven15DF$pop <- rep(outName, length(rownames(eleven15DF)))
sixteen20DF$pop <- rep(outName, length(rownames(sixteen20DF)))
twentyone25DF$pop <- rep(outName, length(rownames(twentyone25DF)))
twentysix30DF$pop <- rep(outName, length(rownames(twentysix30DF)))
thirtyone35DF$pop <- rep(outName, length(rownames(thirtyone35DF)))
thirtysix40DF$pop <- rep(outName, length(rownames(thirtysix40DF)))
fortyAboveDF$pop <- rep(outName, length(rownames(fortyAboveDF)))

write.table(zeroDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_zero.txt"), quote =F,  col.names = F, row.names = F)
write.table(one5DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_one5.txt"), quote =F,  col.names = F, row.names = F)
write.table(six10DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_six10.txt"), quote =F,  col.names = F, row.names = F)
write.table(eleven15DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_eleven15.txt"), quote =F,  col.names = F, row.names = F)
write.table(sixteen20DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_sixteen20.txt"), quote =F,  col.names = F, row.names = F)
write.table(twentyone25DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_twentyone25.txt"), quote =F,  col.names = F, row.names = F)
write.table(twentysix30DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_twentysix30.txt"), quote =F,  col.names = F, row.names = F)
write.table(thirtyone35DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_thirtyone35.txt"), quote =F,  col.names = F, row.names = F)
write.table(thirtysix40DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_thirtysix40.txt"), quote =F,  col.names = F, row.names = F)
write.table(fortyAboveDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_fortyAbove.txt"), quote =F,  col.names = F, row.names = F)

write.table(connectDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_connections_", outName, "_rawConnnects.txt"), quote =F,  col.names = F, row.names = F)

names(member[grep("^1$", member)]) -> cluster1
names(member[grep("^2$", member)]) -> cluster2
names(member[grep("^3$", member)]) -> cluster3
names(member[grep("^4$", member)]) -> cluster4
names(member[grep("^5$", member)]) -> cluster5
names(member[grep("^6$", member)]) -> cluster6
names(member[grep("^7$", member)]) -> cluster7

write.table(cluster1, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster1.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster2, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster2.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster3, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster3.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster4, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster4.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster5, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster5.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster6, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster6.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster7, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_cluster7.txt"), quote =F,  col.names = F, row.names = T)


#Basic features of the network

# Network features
SpiecEasi.nodes <- V(SpiecEasi.net)
SpiecEasi.edges <- V(SpiecEasi.net)
SpiecEasi.node.names <- V(SpiecEasi.net)$name
SpiecEasi.num.nodes <- vcount(SpiecEasi.net)
SpiecEasi.num.edges <- ecount(SpiecEasi.net)


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

plot.net(SpiecEasi.net, SpiecEasi.net.hs, outfile = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_uncolored"), title = outName)

SpiecEasi.AP <- articulation.points(SpiecEasi.net)
#print(SpiecEasi.AP)

outName_degree.dist <- degree_distribution(SpiecEasi.net, mode = "all", cumulative = F)
as.data.frame(outName_degree.dist) -> outDF
write.table(outDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_degree.txt"), quote =F,  col.names = F, row.names = T)

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
  AP <- lapply(names(AP), function(x) x)
  marks <- lapply(1:length(AP), function(x) which(node.names == AP[[x]]))
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
       mark.groups = marks,
       mark.col = "white",
       mark.expand = 10,
       mark.shape = 1,
       layout=l)
  title(title, cex.main=4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  #legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}
# Execute this command after running Function 3
plot.net.cls(SpiecEasi.net, SpiecEasi.net.hs, SpiecEasi.wt, SpiecEasi.AP,
             outfile =  paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/SpiecEasi_", outName, "_colored"), title = outName)
