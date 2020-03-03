# Title     : TODO
# Objective : TODO
# Created by: emanuele
# Created on: 2/13/20


library(tidyverse)
library(plot3D)
#sample_file <- here::here("data/sample.csv")

dat <- read.csv("analysis/histo-time.dat", sep = '', header = FALSE)
s <- 's'
cluster_size <- seq(from = 0, to = (ncol(dat) -1))
cluster_mass <- replace(cluster_size, cluster_size==0, 1)
cluster_mass[2]=0
dat_index <- seq(from = 0, to = (ncol(dat) -1))


#cluster_mass
cluster_size <- paste0(s, "_", cluster_size)


colnames(dat) <- cluster_size
dat <- t(t(dat)*cluster_mass)
#colnames(dat)[s_0] <- "frame"
#dat # fin qui tutto bene

#dat$s_fib <- rowSums(dat[,25:ncol(dat)])
#cbind(dat, total = rowSums(dat))


#not_a_fibril <- seq(1, by = 1, length = ncol(dat) / 2)
not_a_fibril <- subset(dat, select = s_0:s_24)
#not_a_fibril # ok

is_a_fibril <- subset(dat, select = s_25:ncol(dat))
dat <- cbind(dat, total = rowSums(is_a_fibril))

#is_a_fibril

#dat <- cbind(is_a_fibril$total)


#is_a_fibril # ok

#fibril_2 <- rowSums(fibril)
#fibril_2
#length(fibril_2)

#not_a_fibril$s_fib <- list(fibril_2)

#not_a_fibril




#colnames(dat) <- cluster_size
#nmr_molecules
z <- data.matrix(dat)
#colnames(z) <- NULL
#rownames(z) <- NULL
#z <- t(z)*dat_index
#z
#z[,-2]
#dat_index
#z <- z[,-20]
#ncol(z)
#nrow(z)
hist3D(z = z, bty = "g", phi = 30, space = 0.3, ticktype = "detailed", d = 2, theta = 45)
#text3D(x = 1:117, y = 1:117, z = 1:117, labels = rownames(z), add = TRUE)
#length(dat_index)

for_plotly <- rbind(z, dat_index)
for_plotly
write.csv(for_plotly, "analysis/histo-matrix.csv")
#dat

#ncol(dat)
#nrow(dat)
# pesato sulla dimensioni
