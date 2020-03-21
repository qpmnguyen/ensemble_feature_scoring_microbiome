library(profvis)
source("R/laplace_score.R")
profvis(laplace_score(test[,1], weights))
