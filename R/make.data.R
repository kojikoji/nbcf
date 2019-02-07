
library(tidyverse)
library(Matrix)
change.point <- c(0.2, 0.7)
change.var <- list(C1=c(1, 2, 3), C2=c(3, 4, 5, 6))
p.list <- list(
  I1=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  I2=c(0.3, 0.1, 0.9, 0.5, 0.5, 0.5),
  I3=c(0.3, 0.1, 0.2, 0.6, 0.1, 0.8))
r.val <- 30
generate.sample <- function(p.vec) purrr::map(
                                            seq(length(p.vec)),
                                            ~ list(label=.x, val=rnbinom(1, 30, p.vec[.x])))
sim.df <- tibble(t=seq(100)/100) %>% 
  dplyr::mutate(
           I=ifelse(
             t < 0.2,
             "I1",
             ifelse(
               t < 0.7,
               "I2",
               "I3"))) %>%
  dplyr::mutate(
           sample = purrr::map(
                             I,
                             ~ generate.sample(p.list[[.x]]))) %>%
  unnest() %>%
  dplyr::mutate(
           label = purrr::map(
                            sample,
                            ~ .x[[1]]),
           val = purrr::map(
                            sample,
                            ~ .x[[2]])) %>%
  select(-sample)  %>%
  unnest() %>%
  tidyr::spread(label, val)
sim.data <- list(
  t.vec = sim.df$t,
  count.mat = Matrix(as.matrix(sim.df[,3:8]), sparse = T))
"sim.data"

