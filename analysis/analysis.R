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
  scale_fill_gradient(low = "white", high = "dark blue")


# ggplot(data = total_to_plot, aes(x = time, y = value)) + geom_col(aes(fill = value)) +
#   scale_fill_gradient(low = "white", high = "dark blue") # questo ci piace!!

# Matrix of the entire clustsize to make the heatmap
matrix_nclust <- melt(nclust, varnames = c("time", "cluster_size"))
colnames(matrix_nclust)[3] <- "cluster_amount_MW"
matrix_no0 <- matrix_nclust[matrix_nclust$cluster_amount_MW != 0, ]

# QUESTA è UN'IDEA MOLTO CARINA MA STO AVENDO UN ATTIMO DI PROBLEMI
ggplot(data = matrix_no0, aes(x = time, y = cluster_size)) + geom_point(aes(size = cluster_amount_MW,
                                                                            colour = cluster_amount_MW)) +
  #SIZE E COLOR PER CLASSI!!
  #scale_fill_gradient2(low = "white", high = "dark blue") + # Questo va bene
  scale_y_discrete(breaks = c("size_0", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules in a Cluster") +
  scale_x_discrete(name = "Time") +
  labs(title = "Heat map of fibrils elongation") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.y = element_text(size = 7)) +
  scale_size(guide = "none") +
  expand_limits(x = 2, y = 50)

# QUESTO CI PIACE, AL MOMENTO LO LASCIO COSì. SI PUò TOGLIERE IL COLORE IN BASE ALLA DIMENSIONE PERCHé TANTO è 2D
ggplot(data = matrix_no0, aes(x = time, y = cluster_size)) + geom_raster(aes(fill = cluster_amount_MW)) + # ok
  #scale_fill_gradient2(low = "white", high = "dark blue", mid = "light blue", midpoint = 0.8) +
  scale_fill_gradient(low = "white", high = "dark blue") + # Questo va bene
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules in a Cluster") +
  #scale_x_discrete(name = "Time") +
  labs(title = "Heat map of fibrils elongation", fill = "Amount") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.y = element_text(size = 7))


#library(ggplot2)
#ggplot(mtcars, aes(wt, mpg)) +
#  geom_point(aes(colour = cut(qsec, c(-Inf, 17, 19, Inf))),
#             size = 5) +
#  scale_color_manual(name = "qsec",
#                     values = c("(-Inf,17]" = "black",
#                                  "(17,19]" = "yellow",
#                                  "(19, Inf]" = "red"),
#                     labels = c("<= 17", "17 < qsec <= 19", "> 19"))




#ggplot(data = matrix_no0, aes(x = time, y = cluster_size)) + geom_point(aes(size = cluster_amount_MW,
#                                                                            colour = cut(cluster_amount_MW, c (100, 200)))) +
#scale_color_manual(name = "cluster_amount_MW",
#                    values = c("(100,200]" = "yellow"),
#                    labels = c("<=100", "17 <x<19", "<19")) +