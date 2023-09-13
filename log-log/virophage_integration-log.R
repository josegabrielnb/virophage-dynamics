setwd("~/Desktop/JoseGabriel/PLoS_CompBio/r_figures/")

library(tidyverse)
library(plot.matrix)

sputnik <- read_csv("sputnik-logspace_outcomes-1e-5.csv",col_names = c('l','f','o'))

glimpse(sputnik)

sputnik$f <- signif(sputnik$f, 2)

sputnik_tibble <- sputnik %>% pivot_wider(names_from = l, values_from = o) %>% arrange(-row_number()) %>% column_to_rownames(var = "f") 

sputnik_matrix <- as.matrix(sputnik_tibble)

mavirus <- read_csv("mavirus-logspace_outcomes-1e-5.csv",col_names = c('l','f','o'))

glimpse(mavirus)

mavirus$f <- signif(mavirus$f,2)

mavirus_tibble <- mavirus %>% pivot_wider(names_from = l, values_from = o) %>% arrange(-row_number()) %>% column_to_rownames(var = "f")

mavirus_matrix <- as.matrix(mavirus_tibble)

par(mfcol=c(1,2))

plot(mavirus_matrix,col=c("#000000","#0072B2","#E69F00"), main="Mavirus mechanism",xlab=expression(paste("Integration rate (",lambda,")")),ylab="Inhibition parameter (f)",key=NULL,border=NA)

plot(sputnik_matrix,col=c("#000000","#0072B2"), main="Sputnik mechanism",xlab=expression(paste("Integration rate (",lambda,")")),ylab="Inhibition parameter (f)",key=NULL,border=NA)
