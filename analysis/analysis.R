# Title     : TODO
# Objective : TODO
# Created by: Emanuele
# Created on: 26/03/2020

library(tidyverse)
library(reshape2)
library(hexbin)

nclust <- read.table(file = "analysis/histo-time.dat", header = FALSE, dec = ".")

  # V1 = timestep
  # V2 = number of monomers
  # V3 = number of dimers
  # V4 = ...
#time <- nclust[,1]
nclust[,1] <- NULL
nclust[,1] <- NULL

s <- 'size'
clust_size <- seq(from = 1, to = (ncol(nclust)))
clust_size <- clust_size + 1
nclust <- t(t(nclust)*clust_size)
clust_size <- paste0(s, "_", clust_size)
colnames(nclust) <- clust_size

# Subset to create an histogram and compare with wet lab results
histogram <- subset(nclust, select = size_2:size_25)
is_a_fibril <- subset(nclust, select = size_25:ncol(nclust))
histogram <- cbind(histogram, total = rowSums(is_a_fibril))
total <- cbind(total = rowSums(is_a_fibril))
total_to_plot <- melt(total, varnames = c("time", "fibril"))



# Fibril histogram
ggplot(data = total_to_plot, aes(x = time, y = value)) + geom_col(aes(fill = value)) +
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", name = "Fibril MW") +
  theme_bw() +
  scale_y_discrete(name = "Fibril MW")


# Matrix of the entire clustsize to make the heatmap
matrix_nclust <- melt(nclust, varnames = c("time", "cluster_size"))
colnames(matrix_nclust)[3] <- "cluster_amount_MW"
matrix_nclust[matrix_nclust == 0] <- NA

# QUESTA è UN'IDEA MOLTO CARINA MA STO AVENDO UN ATTIMO DI PROBLEMI
ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) + geom_point(aes(size = cluster_amount_MW,
                                                                               colour = cluster_amount_MW),
                                                                           alpha = 1/50) +
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600", "size_700"),
                   name = "Molecules in a Cluster") +
  scale_radius(guide = "none") +
  scale_colour_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", name = "Amount", guide = "none") +
  theme_bw() +
  expand_limits(y = 700)


# QUESTO CI PIACE, AL MOMENTO LO LASCIO COSì. SI PUò TOGLIERE IL COLORE IN BASE ALLA DIMENSIONE PERCHé TANTO è 2D
ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) + geom_raster(aes(fill = cluster_amount_MW)) + # ok
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", guide = "none") + # Questo va bene
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules in a Cluster") +
  labs(title = "Heat map of fibrils elongation", fill = "Amount") +
  expand_limits(y = 700) +
  theme_bw()


# MERGE
ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) + geom_raster(aes(fill = cluster_amount_MW)) + # ok
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", guide = "none") + # Questo va bene
  labs(title = "Heat map of fibrils elongation") +
  geom_point(aes(size = cluster_amount_MW, colour = cluster_amount_MW), alpha = 1/40) +
  scale_radius(guide = "none") +
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules in a Cluster") +
  scale_colour_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", guide = "none") +
  expand_limits(y = 700) +
  theme_bw()