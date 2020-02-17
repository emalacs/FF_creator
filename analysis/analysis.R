# Title     : TODO
# Objective : TODO
# Created by: emanuele
# Created on: 2/13/20


library(tidyverse)
library(plot3D)
#sample_file <- here::here("data/sample.csv")

dat <- read.csv("analysis/histo-time.dat", sep = '', header = FALSE)

nmr_molecules <- seq(from = 0, to = (ncol(dat) -1))
size <- 's'
#nmr_molecules <- paste0(size, "_", nmr_molecules)
#nmr_molecules <- paste(nmr_molecules, "_size", sep = "")
dat_index <- seq(from = 0, to = (ncol(dat) -1))




colnames(dat) <- nmr_molecules
#nmr_molecules
dat
z <- data.matrix(dat)
#colnames(z) <- NULL
#rownames(z) <- NULL
#z
dat_index[2]=0
z <- t(z)*dat_index
z
#z[,-2]
#dat_index
#z <- z[,-20]
#ncol(z)
#nrow(z)
hist3D(z = z, theta = -180)

write.csv(z, "analysis/histo-matrix.csv")
#dat

#ncol(dat)
#nrow(dat)
# pesato sulla dimensioni
