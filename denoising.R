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

#  ----- Load packages --------------------------------------------------------

library(dplyr)
library(reticulate)
library(zoo)
library(bcp)
library(rlist)

# ----- Preprocessing of data --------------------------------------------------

# Normalized one sample into [0,1]
scale.data <- \(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Centers one sample around center.value
center.data <- \(x, center.value = 0.5) {
  x + center.value - mean(x)
}

# Normalize a sample in training dataset into [0,1]
scale.train.test <- \(true.data, noisy.data) {
  min.value <- min(noisy.data)
  max.value <- max(noisy.data)
  list(segment = (true.data - min.value) / (max.value - min.value),
       CN = (noisy.data - min.value) / (max.value - min.value))
}

# Rescale a denoised to the original scale
rescale.denoised <- \(denoised, original) {
  min <- min(original)
  max <- max(original)
  (max-min)*denoised+min
}

# Rescale a denoised dataset to the original scale
rescale.denoised.df <- \(denoised.df, original.df) {
  denoised.df <- sapply(1:nrow(denoised.df), \(i) 
                     rescale.denoised(denoised.df[i,] %>% as.matrix, 
                                      original.df[i,] %>% as.matrix)) %>% 
    t %>% as.data.frame
  rownames(denoised.df) <- rownames(original.df)
  return(denoised.df)
}

# Scale each sample in dataset into [0,1]
scale.dataset <- \(simulated.data, center = FALSE) {
  true <- simulated.data$segment
  noisy <- simulated.data$CN
  if (center) {
    noisy <- apply(noisy, 1, center.data)
    true <- apply(true, 1, center.data) 
  }
  for (i in 1:nrow(true)) {
    scaled.data <- scale.train.test(true[i, ], noisy[i, ])
    true[i, ] <- scaled.data$segment
    noisy[i, ] <- scaled.data$CN
  }
  list(segment = true, CN = noisy, purity = simulated.data$purity)
}

scale.mode <- \(cn.data, mode = 2) {
  density.est <- density(cn.data %>% as.matrix %>% t)
  modes <- .compute.modes(density.est)
  modes.max <- which.max(density.est$y) # find index of highest peak
  max.x <- density.est$x[modes.max]
  cn.data / max.x * mode
}

# shuffles the order of the rows in the dataset
shuffle.dataset <- \(dataset) {
  CN <- dataset$CN 
  segments <- dataset$segment 
  N <- nrow(CN) # each column is a sample
  order <- sample(1:N, size = N)
  list(CN = CN[order, ], segment = segments[order, ],
       purity = dataset$purity[order])
}

# ----- Plots ------------------------------------------------------------------

# plot simulated data, one sample
plot.simulated.data <- \(dataset, sample = 1) {
  N.features <- ncol(dataset$CN)
  data.frame(
    signal = c(dataset$CN[sample, ], dataset$segment[sample, ]),
    type = rep(c("noisy", "true"), each = N.features),
    index = rep(1:N.features, 2)
  ) %>% ggplot(aes(x = index, y = signal, color = type)) + geom_line() +
    scale_color_manual(values = c("lightgray", "red")) + labs(color = NULL) +
    theme(legend.position = "top")
}

# plot multiple samples in the same plot
plot.multiple.samples <- \(data, samples.to.plot = seq(1,3)) {
  plots <- lapply(samples.to.plot, \(sample) {
    plot.simulated.data(data, sample) + theme(legend.position = "none")
  }) 
  prow <- plot_grid(plotlist = plots, ncol = 1, align = "vh")
  legend <- get_legend(
    plots[[1]] + guides(color = guide_legend(ncol = 1)) + 
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend, rel_widths = c(3, .4))
}

# plot reconstruction of denoised signal from model
plot.reconstruction <- \(signal, reconstruction, true.segment = NULL) {
  signal %<>% as.numeric
  reconstruction %<>% as.numeric
  df <- data.frame(
    true.signal = c(reconstruction, signal), 
    index = rep(1:length(signal), 2),
    label = rep(c("reconstruction", "signal") %>% 
                  as.factor, each = length(signal)) %>%
      relevel(ref = "signal")
  )
  colors <- c("lightgray", "red")
  if(!is.null(true.segment)) {
    df.true <- data.frame(
      true.signal = true.segment, 
      index = 1:length(true.segment),
      label = "true segmentation"
    )
    df <- rbind(df, df.true)
    df$label <- factor(df$label) %>% reorder(., rep(c(3,1,2), 
                                                    each = length(signal))) 
    colors <- c("lightgray", "azure4", "red")
  }
  gg <- ggplot(data = df, aes(x = index, y = true.signal, color = label)) + 
    geom_line() + scale_color_manual(values = colors) + labs(color = NULL) +
    theme(legend.position = "top")
  gg
}

# plots the validation and training loss from training the network
plot.loss <- \(list.losses) {
  training.loss <- list.losses$loss
  validation.loss <- list.losses$val_loss
  data.frame(
    loss = c(training.loss, validation.loss),
    loss.type = rep(c("training", "validation"), each = length(training.loss)),
    epoch = rep(1:length(training.loss), 2)
  ) %>% ggplot(aes(x = epoch, y = loss, color = loss.type)) + geom_line() +
        geom_point() + labs(color = NULL) + theme(legend.position = "top")
}

# ----- Code for running CNN ---------------------------------------------------

# the model should be an object returned by the run.autoencoder function
run.prediction <- \(model, new.data, scale = FALSE, edge.smoothing = F, 
                    n.smooth = 32) {
  if(scale) new.data %<>% scale.data
  pred <- py$predict_observations(model, new.data)
  if(edge.smoothing) {
    len <- length(pred)
    pred[1:n.smooth] <- pred[n.smooth]
    pred[(len-n.smooth):len] <- pred[(len-n.smooth)]
  }
  return(pred)
}

# ---- Train model -------------------------------------------------------------

# Run conv1d training
# original.data: the ground truth (copy number profiles)
# noisy.data: observed data
# penalty: penalty parameter to use in loss function
# pooling.method: pooling method to use
# ker.size: size of kernel
# save.path: path to save the trained model in
run.conv1d.denoiser <- \(original.data, noisy.data, penalty = 0.025, 
                         feature.size = 128, n.layers = 2, 
                         pooling.method = "max", ker.size = 7,
                         save.path = "denoising_model") {
  n.samples <- nrow(original.data)
  py$run_conv1d_denoiser(
    original.data,
    noisy.data,
    penalty = penalty,
    feature_size = feature.size,
    n_layers = n.layers,
    pooling_method = pooling.method,
    kernel_size = ker.size,
    save_path = save.path
  )
}

# load a model saved to disk to avoid training
load.denoiser <- \(path = "denoising_model") {
  py$load_model(path)
}

# print information about the model, e.g. number of parameters
print.model.summary <- \(model) {
  py$print_model_summary(model)
}

# pipeline for running the conv1d model
# dataset: data to train on
# model: model function to use
# pre.process: how to preprocess the data before training
run.pipeline <- \(dataset, model = run.conv1d.denoiser, 
                  pre.process = scale.dataset, ...) {
  if (!is.null(pre.process)) {
    dataset %<>% pre.process
  }
  original.data <- dataset$segment %>% as.matrix
  noisy.data <- dataset$CN %>% as.matrix
  model(original.data, noisy.data, ...)
}

# ----- BCP --------------------------------------------------------------------

# computes the bcp segmentation of an object and returns the posterior mean and
# the posterior probabilities, along with the data
# data.t: a matrix or dataframe with samples as columns and genomic bins as rows
#         NOTE that this is transposed compared to the usual format
compute.segmentation.bcp <- \(data.t, ...) {
  mcmc.object <- bcp(data.t, ...)
  mcmc.object %$% 
    list(segment = posterior.mean, probs = posterior.prob, data = data[, -1])
}

# post process the bcp segmentation
# mcmc.object: the object returned by compute.segmentation.bcp
# cutoff: the cutoff for the posterior probability
# min.bin: the minimum number of bins in a segment
# trim: the fraction of bins to trim from the start and end of a segment
# show.plots: whether to show the plots
post.process.bcp <- \(mcmc.object, cutoff = 0.05, min.bin = 12, trim = 0.1, 
                      show.plots = FALSE) {
  is.break <- (mcmc.object$probs > cutoff) * 1
  segments <- segment.df(is.break) # find start and en
  cna.profile <- mcmc.object$segment
  # begin by just taking the median
  for (i in 1:nrow(segments)) {
    # # TODO : compute the value
    inds <- segments[i, ] %$% start:end
    n.bins <- length(inds)
    n.to.remove <- floor(trim * n.bins) # floor to handle short segments
    inds.trimmed <- (inds[1] + n.to.remove):(inds[n.bins] - n.to.remove)
    segments[i, ]$value <- median(cna.profile[inds.trimmed])
    # the value should be the trimmed median 
  }
  segments %<>% postprocess.segment #%<>% filter.segments
  if(show.plots) {
    plot.bcp.diagnostics(mcmc.object, cutoff = cutoff, 
                         post.processed.profile = segments %>% get.profile) %>% 
      print
  }
  segments.adj <- seg.from.bps(mcmc.object$data, segments)
  if(show.plots) {
    plot.segments(mcmc.object$data, segmentation = segments.adj) %>% print# +
  }
  segments.adj
  # segments
}

# plot the diagnostic plots from a bcp object
# mcmc.obj: the object returned by compute.segmentation.bcp
# cutoff: the cutoff for the posterior probability
# true.segment: the true segment to plot
# post.processed.profile: the post processed profile to plot
plot.bcp.diagnostics <- \(mcmc.obj, cutoff = 0.05, true.segment = NULL, 
                          post.processed.profile = NULL) {
  data <- mcmc.obj$data
  x <- 1:length(data)
  if (is.data.frame(data)) {
    data %<>% as.numeric
  }
  df.point <- data.frame(
    x = x,
    y = data
  )
  p1 <- data.frame(
    x = x,
    y = mcmc.obj$segment %>% as.numeric
  ) %>% ggplot(aes(x, y)) +
    geom_point(data = df.point, aes(x, y), color = "lightgray") +
    geom_line(color = "red") + labs(x = "genomic bin", y = "posterior mean")
  if (!is.null(true.segment)) {
    df.true <- data.frame(x = x, y = true.segment %>% as.numeric)
    p1 <- p1 + geom_line(data = df.true, aes(x = x, y = y), 
                         linetype = "dashed")
  }
  p2 <- data.frame(
    x = x,
    y = mcmc.obj$prob %>% as.numeric
  ) %>% ggplot(aes(x, y)) + geom_line() + 
    geom_hline(yintercept = cutoff, color = "red") + 
    labs(x = "genomic bin", y = "posterior probability")
  plot.list <- list(p1, p2)
  if (!is.null(post.processed.profile)) {
    p3 <- data.frame(
      x = x,
      y = post.processed.profile %>% as.numeric
    ) %>% ggplot(aes(x, y)) +
      geom_point(data = df.point, aes(x, y), color = "lightgray") +
      geom_line(color = "red") + 
      labs(x = "genomic bin", y = "post processed profile")
    if (!is.null(true.segment)) {
      df.true <- data.frame(x = x, y = true.segment %>% as.numeric)
      p3 <- p3 + geom_line(data = df.true, aes(x = x, y = y), 
                           linetype = "dashed")
    }
    plot.list <- list.append(plot.list, p3)
  }
  plot_grid(plotlist = plot.list, ncol = 1)
}

# ----- Rolling average/median -------------------------------------------------

roll <- \(noisy, k = 121, stop = 0.001, max.iters = 100,
          edge.smoothing = TRUE, verbose = FALSE, method = "median") {
  if(method == "median") {
    rollfun <- rollmedian
  } else if(method == "mean" | method == "avg") {
    rollfun <- rollmean
  } else {
    stop("invalid 'method'. options are 'median' or 'mean'.")
  }
  ix <- floor(k/2)
  denoised <- noisy
  i <- 0
  for(i in 1:max.iters) {
    i <- i+1
    rolled <- rollfun(noisy, k = k)
    if(mean(denoised[(ix+1):(length(noisy)-ix)] - rolled) < stop) {
      denoised[(ix+1):(length(noisy)-ix)] <- rolled
      break
    }
    denoised[(ix+1):(length(noisy)-ix)] <- rolled
  }
  if(edge.smoothing) {
    denoised[1:ix] %<>% median
    denoised[(length(noisy)-ix+1):length(noisy)] %<>% median
  }
  if(verbose) {
    print(sprintf("number of iterations: %i", i))
  }
  return(denoised)
} 

# computes rolling median or mean on dataset
# noisy.df: data to denoise
# k: size of sliding window to use
# stop: stopping criteria for denoising
# max.iters: number of iterations to run denoising
# edge.smoothing: TRUE if edges should be padded
# method: method to use, either median or mean
roll.df <- \(noisy.df, k = 121, stop = 0.001, max.iters = 100,
              edge.smoothing = TRUE, verbose = FALSE, method = "median") {
  apply(noisy.df, 1, \(noisy) roll(noisy, k, stop, max.iters,
                                   edge.smoothing, verbose, method)) %>% 
    t %>% as.data.frame
}
