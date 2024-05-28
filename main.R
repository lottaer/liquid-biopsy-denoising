#     Estimation of purity and subclonal ratio from liquid biosy sequencing. 
#     Inspired by Lakatos et. al., <https://github.com/elakatos/liquidCNA>.
#     Copyright (C) 2024  Lotta Eriksson & Linnea Hallin
#                    lottaer@chalmers.se   hallinl@chalmers.se
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ----- Setup ------------------------------------------------------------------
## ----- Load packages ---------------------------------------------------------

library(QDNAseq)
library(ggplot2)
library(magrittr)
library(rstudioapi)
library(cowplot)
library(forcats)
library(pbapply)
library(dplyr)

# set the working directory to source file location, and print location
getActiveDocumentContext()$path %>% dirname %>% setwd %T>% print 

##  ----- Settings -------------------------------------------------------------
# plot settings
theme_set(theme_bw(base_size = 14)) # set default theme and text size
cb.palette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", 
  "#D55E00", "#CC79A7"
) # color-blind friendly palette
options(ggplot2.discrete.colour = cb.palette) # set cb.palette as default

# progress bar settings
pboptions(char = "=") 

# dplyr settings 
options(dplyr.summarise.inform = FALSE) # suppress output in dplyr summarise

## ----- Load required functions -----------------------------------------------

# source all R functions
source("handle_data.R") 
source("denoising.R")   
source("simulate_data.R")
source("segmentation.R")
source("estimation_functions.R")
source("autotuning.R")
source("evaluation.R")
# source functions for constructing denoising model
reticulate::use_virtualenv("./env")
reticulate::source_python("models.py")

## ----- Load data about segments and chromosome position ----------------------

genome.info <- QDNAseq::getBinAnnotations(500)@data

# ----- Data exploration (outlier detection and removal) -----------------------
# load the data
cn.samples <- files.to.dataset("example_data/cell-line", remove.blacklist = F) %>% 
  remove.missing

# find outliers below the median, plot
outliers.lower <- find.outliers.df(cn.samples, region = "lower")
plot.outliers(outliers.lower) + theme(text = element_text(size = 8))
plot.outliers(outliers.lower, 1500:2000) # zoom in on problematic region

# find the outliers and their position in the original genome (before removal 
# of NA's)
blacklist.id <- 1737:1742
outliers.names <- colnames(cn.samples[, blacklist.id])
blacklist <- which(colnames(cn.samples) %in% outliers.names)

# look at the data again after removing the outliers
cn.samples <- cn.samples[, -blacklist.id]
outliers.blacklist <- find.outliers.df(cn.samples, region = "lower")

# there are no more obvious outlier patterns :)
plot.outliers(outliers.blacklist)

# ----- Tuning -----------------------------------------------------------------
## ----- Rolling median tuning -------------------------------------------------

# create the tuning data
tuning.data <- simulate.dataset(
  20, 1, seq(0.05, 0.5, length.out = 10), N.segments = 4320,
  noise.mean = c(0.1, 0.2, 0.4, 0.8), noise.sd = c(0.01, 0.02, 0.04, 0.08),
  cna.dep.noise = T, downsample = T
)

k.seq <- c(21, 31, 41, 51, 61) # size of sliding window to try
roll.results <- lapply(k.seq, \(k) {
  # repeat the same procedure as for the CNN
  print(sprintf("[INFO] Running tuning with value k = %i", k))
  res.df <- compute.mse.pipeline(tuning.data, compute.segmentation = T, k = k)
  res.df$k <- k
  res.df$noise <- tuning.data$noise.mean
  res.df
}) %>% do.call(rbind, .)
# plot the results and select an appropriate model
plot.roll.tuning(roll.results) 

## ---- Denoising autoencoder tuning -------------------------------------------

# create the tuning data
tuning.data.denoiser <- simulate.dataset(
  200, 1, seq(0.05, 0.5, length.out = 10), N.segments = 4320,
  noise.mean = c(0.1), noise.sd = c(0.01),
  cna.dep.noise = T, downsample = T
)

lambdas <- c(0, 0.025, 0.05, 0.075, 0.1)
for(penalty in lambdas) {
  results <- run.tuning(tuning.data.denoiser, penalty)
  # save the results since the tuning process is very time consuming
  saveRDS(results, paste0("data/results_", penalty, ".RData"))
}

# read all the tuning files from disk
tuning.files <- list.files("data", pattern = "results", full.names = TRUE) 
# concatenate all results into data frame
results <- lapply(tuning.files, readRDS) %>% do.call(rbind, .)

# OPTIONAL: plot the validation loss for each penalty value

limits <- tapply(results$val.loss, results$penalty, 
                 \(lambda) quantile(lambda, c(0, 0.95)))
# plot the validation loss for each penalty value
lapply(1:length(lambdas), \(i) {
  plot.models.custom(results, ylim = c(0, 0.01), lambda = lambdas[i],
                     remove.outlier = T)
})
# the models are not directly comparable, since a larger penalty
# will result in a larger validation loss. Hence, we select the best model for
# each penalty value

## ----- Train the best models again on the entire dataset ---------------------

# fit and save the new models
fit.best.models(results, tuning.data.denoiser, save.models = T, n.models = 1)

# load the models from disk
models <- load.best.models(n.best = 1)
sprintf("Number of models to tune: %i", length(models))

# select model based on mse
df.mse <- run.penalty.tuning(models, tuning.data, lambdas = lambdas,
                             pipeline.fun = "mse", compute.seg = T)

# select the best model based on plots
plot.penalty.tuning(df.mse)

## ----- Train selected model again on more diverse data -----------------------

# the final model was retrained on more, and more diverse data
data.retrain <- simulate.dataset(
  N.each = 200, N.replicates = 1, 
  cancer.prop = seq(0.05, 0.3, length.out = 26), 
  N.segments = 4320, min.length = 12, max.length = 130,
  cna.dep.noise = T, noise.mean = 0.1
) %>% scale.dataset %>% shuffle.dataset

# the default values corresponds to the model selected in our work
model <- run.conv1d.denoiser(
  data.retrain$segment, data.retrain$CN
)

# ----- Denoising autoencoder example ------------------------------------------

# print the model summary and plot the loss
print.model.summary(model$model)
plot.loss(model$history) # plot the training and validation loss

new.observation <- 
  simulate.dataset(N.each = 1, N.replicates = 1, cancer.prop = 0.3, 
                   N.segments = 4320, min.length = 12, max.length = 130,
                   cna.dep.noise = T, noise.mean = 0.2, noise.sd = 0,
                   downsample = T) %>% 
  scale.dataset

denoised <- run.prediction(
  model$model, new.observation$CN[1,], edge.smoothing = T
)

# plot the prediction without segmentation
plot.reconstruction(
  new.observation$CN[1,], denoised, new.observation$segment[1,]
)

# plot the prediction with segmentation
plot.reconstruction(
  new.observation$CN[1,], denoised %>% compute.segmentation(show.plot = F) %>% 
    get.profile, new.observation$segment[1,]
)

# ----- BCP example ------------------------------------------------------------
# simulate new observation
new.observation <- 
  simulate.dataset(N.each = 1, N.replicates = 1, cancer.prop = 0.3, 
                   N.segments = 4320, min.length = 12, max.length = 130,
                   cna.dep.noise = T, noise.mean = 0.2, noise.sd = 0,
                   downsample = T) 

# compute bcp and the postprocessing and plot the diagnostics
bcp.obj <- compute.segmentation.bcp(new.observation$CN[1,])
bcp.post <- post.process.bcp(bcp.obj, cutoff = 0.1)
plot.bcp.diagnostics(bcp.obj, 0.1, new.observation$segment,
                     bcp.post %>% get.profile)

# ----- Evaluation -------------------------------------------------------------
# simalate data for evaluation
eval.data <- simulate.dataset(
  5, 1, seq(0.05, 0.5, by = 0.05), 4320, 
  noise.mean = c(0.1), noise.sd = c(0),
  cna.dep.noise = T, downsample = T
) %>% shuffle.dataset

# our model
model <- model$model
denoised.nn <- 
  apply(eval.data$CN, 1, 
        \(sample) run.prediction(model, sample, T, T)[,,1] %>%
          rescale.denoised(sample)) %>% t
# segment the data
segmented.nn <- 
  pbapply(denoised.nn, 1, \(row) compute.segmentation(row, show.plot = F) %>% 
            get.profile) %>% t
purities.nn <- 
  pbapply(denoised.nn, 1, \(s) compute.purity(s)[2, ])
nn.corrected <- correct.eval(segmented.nn, purities.nn)

# rolling median 
denoised.rollmed <- apply(eval.data$CN, 1, roll, k = 21) %>% t
segmented.rollmed <- 
  pbapply(denoised.rollmed, 1, \(row) compute.segmentation(row, show.plot = F) %>% 
            get.profile) %>% t
purities.rollmed <- 
  pbapply(denoised.rollmed, 1, \(s) compute.purity(s)[2, ])
rollmed.corrected <- correct.eval(segmented.rollmed, purities.rollmed)

# bcp
bcp.obj <- pbapply(eval.data$CN %>% t, 2, compute.segmentation.bcp)
denoised.bcp <- lapply(bcp.obj, \(obj) obj$segment) %>% do.call(cbind, .) %>% t
purity.bcp <- pbapply(denoised.bcp, 1, \(s) compute.purity(s)[2, ])
segmented.bcp <- 
  lapply(bcp.obj, \(o) post.process.bcp(o, show.plots = F) %>% get.profile) %>% 
  do.call(rbind, .)
bcp.corrected <- correct.eval(segmented.bcp, purity.bcp)

# metadata
purity <- eval.data$purity
noise <- eval.data$noise.mean
true.data <- eval.data$segment
true.corr <- 
  sapply(1:nrow(true), \(i) correct.profile(true[i, ], purity[i])) %>% t

corrected.list <- list(nn.corrected, rollmed.corrected, bcp.corrected)
error.vec <- 
  lapply(corrected.list, \(mat) compute.mse.matrix(mat, true.corr)) %>% unlist
N <- nrow(true.corr)
df.eval <- data.frame(
  error = error.vec,
  purity = cut(purity, seq(0, 0.5, by = 0.1)),
  noise = eval.data$noise.mean,
  method = rep(c("autoencoder", "rolling median", "bcp"), each = N)
)
df.eval$method <- 
  factor(df.eval$method, levels = c("bcp", "autoencoder", "rolling median"))

# ----- Purity estimation ------------------------------------------------------
## ----- Simulated data --------------------------------------------------------
# simulate data
data.sim <- simulate.dataset(
  N.each = 200, N.replicates = 1, 
  cancer.prop = seq(0.05, 0.3, length.out = 26), 
  N.segments = 4320, min.length = 12, max.length = 130,
  cna.dep.noise = T, noise.mean = 0.1
)

# autoencoder
model <- load.denoiser()
denoised.nn <- apply(data.sim$CN, 1, 
                     \(sample) run.prediction(model, sample, T, T)[,,1] %>%
                       rescale.denoised(sample)) %>% t
purities.nn <- pbapply(denoised.nn, 1, \(s) compute.purity(s)) %>% 
  do.call(cbind, .) %>% t %>% as.data.frame

# rolling median 
denoised.rollmed <- apply(data.sim$CN, 1, roll, k = 21) %>% t
purities.rollmed <- pbapply(denoised.rollmed, 1, \(s) compute.purity(s)) %>% 
  do.call(cbind, .) %>% t %>% as.data.frame

# bcp
bcp.obj <- apply(data.sim$CN %>% t, 2, compute.segmentation.bcp)
denoised.bcp <- lapply(bcp.obj, \(obj) obj$segment) %>% do.call(cbind, .) %>% t
purity.bcp <- pbapply(denoised.bcp, 1, \(s) compute.purity(s)) %>% 
  do.call(cbind, .) %>% t %>% as.data.frame

# plot the results
purity <- data.sim$purity
plot.purity(data.frame(purities.nn, true = purity))
plot.purity(data.frame(purities.rollmed, true = purity))
plot.purity(data.frame(purities.bcp, true = purity))

## ----- Insilico data ---------------------------------------------------------
# import the insilico data and metadata
insilico <- readRDS("data/insilico.RData") %>% as.data.frame
params.ins <- readRDS("data/insilico_metadata.RData")
methods = c("nn", "rollmed", "bcp")

# compute the purities and plot
pur.ins <- lapply(methods, \(m) run.purity.insilico(insilico, params.ins, m, 
                                                    model = model)) %>%
  do.call(rbind, .)
plot.purity.insilico(pur.ins)

# ----- Subclonal tracking -----------------------------------------------------
## ----- Simulated data --------------------------------------------------------
# simulate parameters, choose noise levels
params.sim <- lapply(1:200, \(i) sample.params.long(5))
noise.sim <- c(1, 2, 4)

# simulate data, preprocess and run the pipeline for all noise levels
ratios.sim <- lapply(noise.sim, \(sigma) {
  data <- lapply(params.sim, \(p) {
    simulate.long(4320, p$purity, p$ratio, noise.mean = 0.1*sigma, 
                  noise.sd = 0, cna.dep.noise = T, downsample = T)
  })
  pre <- pblapply(data, preprocess.subcl, insilico = F) 
  rs <- subcl.tracking.pipeline(data, pre, cutoff.opt = T)
  # optional: save data and pre so you don't have to rerun everything
  # saveRDS(data, paste0("path/data_", sigma, ".RData"))
  # saveRDS(pre, paste0("path/pre_", sigma, ".RData"))
  rs
})

# plot the results
lapply(1:length(noise.sim), \(i) {
  plot.ratio(ratios.sim[[i]]) + ggtitle(paste0("noise level: ", noise.sim[i]))
})

## ----- Insilico data ---------------------------------------------------------
# import the insilico data and metadata
insilico <- readRDS("example_data/insilico.RData")
params.ins <- readRDS("example_data/insilico_metadata.RData")
noise.ins <- c(50, 20, 10, 5)

# simulate data, preprocess and run the pipeline for all noise levels
ratios.ins <- lapply(1:length(noise.ins), \(sigma) {
  data <- sample.dataset.insilico(insilico, params.ins, 200, noise.id = sigma)
  pre <- pblapply(data, preprocess.subcl, insilico = T)
  rs <- subcl.tracking.pipeline(data, pre, cutoff.opt = T,
                                adj = c(seq(0.25, 2.5, 0.25)))
  # optional: save data and pre so you don't have to rerun everything
  # saveRDS(data, paste0("path/data_", sigma, ".RData"))
  # saveRDS(pre, paste0("path/pre_", sigma, ".RData"))
  rs
})

# plot the results
lapply(1:length(noise.ins), \(i) {
  plot.ratio(ratios.ins[[i]]) + 
  ggtitle(paste0("min # of reads (million): ", noise.ins[i]))
})
