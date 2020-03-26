# Title     : TODO
# Objective : TODO
# Created by: Emanuele
# Created on: 26/03/2020

library(tidyverse)
library(reshape2)

nclust <- read.table(file = "analysis/histo-time.dat", header = FALSE, dec = ".")

  # V1 = timestep
  # V2 = number of monomers
  # V3 = number of dimers
  # V4 = ...
time <- nclust[,1]
nclust[,1] <- NULL
nclust[,1] <- NULL

s <- 'size'
clust_size <- seq(from = 1, to = (ncol(nclust)))
clust_size <- clust_size + 1
#clust_size <- replace(clust_size, clust_size==0, 1)
nclust <- t(t(nclust)*clust_size)
clust_size <- paste0(s, "_", clust_size)
colnames(nclust) <- clust_size
#colnames(nclust)[1] <- "time"

matrix_nclust <- melt(nclust)

ggplot(data = matrix_nclust, aes(x = Var1, y = Var2)) + geom_raster(aes(fill = value))


#nclust_matrix = data.matrix(nclust)

# Prova grafico con due variabili (anche se Cri ha detto che sono 3)

#nclust
#ggplot(data = nclust) + geom_raster(mapping = aes(fill = nclust, x = "time", y = "size_1"))


clust_size




