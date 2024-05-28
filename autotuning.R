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

# ----- Load packages ----------------------------------------------------------

library(gridExtra)
library(ggpubr)

# ----- Tuning functions -------------------------------------------------------

# function for making folds to use in K-fold CV
make.folds <- \(data, K) {
  folds_raw <- seq(1, nrow(data)) %>%
    .[sample(length(.))] %>%
    cut(breaks = K, labels = FALSE)
  lapply(1:K, \(i) which(folds_raw == i, arr.ind = TRUE))
}

# define total variation loss
tv.mse <- \(true, pred, penalty){
  mse = mean((true - pred)^2)
  tv = mean(abs(diff(pred)))
  return(mse + penalty * tv)
}

# tuning function to be run on ozzy
run.tuning <- \(data, penalty, K = 10,
                     n.layers = seq(2,3), 
                     initial.size = c(32, 64, 128), 
                     pooling.method = c('max'), 
                     ker.size = seq(5,9,2)) {
  lapply(pooling.method, \(method) {
    lapply(initial.size, \(size) {
      max.layers <- min(max(n.layers), log2(size) - 1)
      lapply(1:max.layers, \(layers) {
        lapply(ker.size, \(ker){
          folds <- make.folds(data$CN, K)
          pblapply(
            1:K,
            \(k) {
              # split the data
              train <- list(CN = data$CN[-folds[[k]], ], 
                            segment = data$segment[-folds[[k]], ],
                            purity = data$purity[-folds[[k]]])
              val <- list(CN = data$CN[folds[[k]], ], 
                          segment = data$segment[folds[[k]], ],
                          purity = data$purity[folds[[k]]])
              # train the model on training data
              model.obj <- run.pipeline(train, 
                                        feature.size = size, 
                                        penalty = penalty, 
                                        n.layers = layers, 
                                        pooling.method = method, 
                                        ker.size = ker %>% as.integer,
                                        save.path = NULL)
              # run model prediction on validation data
              denoised <- run.prediction(model.obj$model, val$CN)[,,1]
              val.loss <- tv.mse(val$segment, denoised, penalty)
              data.frame(penalty = penalty, 
                         val.loss = val.loss,
                         train.val.loss = min(model.obj$history$val_loss),
                         pooling.method = method,
                         feature.size = size,
                         layers = layers,
                         kernel.size = ker,
                         iteration = k)
            }) %>% do.call(rbind, .)
        }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

# ----- Helper function for analyze results  -----------------------------------

# Get the best models from the tuning results dataframe 
get.best.models <- \(df, n.per.lambda = "all") {
  lambdas <- unique(df$penalty)
  lapply(lambdas, \(l) {
    # extract current lambda values and compute the mean validation loss
    df.curr <- df[df$penalty == l, ] %>% 
      group_by(layers, feature.size, kernel.size, pooling.method, .add = T) %>%
      summarize(val.loss.sd = sd(val.loss), val.loss = mean(val.loss)) %>%
      as.data.frame
    df.curr$penalty <- rep(l, nrow(df.curr))
    df.curr <- df.curr[order(df.curr$val.loss, decreasing = F), ]
    if (n.per.lambda != "all") {
      # select the top n models (in terms of minimal validation loss)
      df.curr %<>% .[1:n.per.lambda, ] 
    }
    df.curr
  })
}

# Loads the n best files from disk
load.best.models <- \(model.path = "models", n.best = 1) {
  # recall that the last position in the name is the rank
  all.models <- list.files(model.path, full.names = T)
  model.names <- 
    sapply(all.models, \(file) strsplit(file, ".", fixed = T)[[1]][1:2]) %>% .[1]
  ranks <- sapply(model.names, \(file) strsplit(file, "_")[[1]][3]) %>%
    as.numeric
  print(sapply(model.names, \(file) strsplit(file, "_")[[1]][3]))
  if (max(ranks) < n.best) {
    warning(sprintf(
      "Too many models requested. Loading %i best models", max(ranks)))
  }
  lapply(all.models[ranks <= n.best], load.denoiser)
}

# ---- Refit the models --------------------------------------------------------

# Retrains the n best models and saves them for further analysis
fit.best.models <- \(df, data, n.models = 1, save.models = T, 
                     save.path = "models") {
  # concatenate the models into one df to make it easier to work with
  df.models <- get.best.models(df, n.per.lambda = n.models) %>% 
    do.call(rbind, .)
  df.models$rank <- rep(1:n.models, nrow(df.models) / n.models)
  # train the models on the new data, and return as list. the models are saved
  # in the models folder on disk
  pblapply(1:nrow(df.models), \(row) {
    setup <- df.models[row, ]
    # naming convention : model_penalty_rank (use print.summary for params)
    save.path <- setup %$% 
      paste0(save.path, "/model_", penalty, "_", rank, ".keras")
    print(sprintf("Model saved on path %s", save.path))
    # use the into in the df row to construct the model
    # the model will be saved in the models folder
    setup %$% 
      run.pipeline(
        dataset = data %>% shuffle.dataset, # make sure to shuffle dataset!!
        feature.size = feature.size,
        n.layers = layers, 
        ker.size = kernel.size, 
        penalty = penalty, # ugh
        pooling.method = pooling.method,
        save.path = ifelse(save.models, save.path, NULL) # NULL if not save
      )$model
    # we discard the history in this case
  })
}

# ----- Pipeline for computing MSE of dataset ----------------------------------

# Computes the MSE of predicitons of entire dataset
compute.mse.pipeline <- \(dataset, method = "rollmed", model = NULL, 
                          compute.segmentation = FALSE, k = 11,
                          ..., show.plots = F) {
  nrow <- nrow(dataset$CN)
  if(length(unique(rownames(dataset$CN))) != nrow(dataset$CN)) {
    rownames(dataset$CN) <- 1:nrow %>% as.factor
  }
  print("[INFO] Denoising...")
  if(method == "conv1d") {
    if (is.null(model)) {
      stop("Method conv1d selected, but not no model provided. Stopping...")
    }
    denoised <- pbapply(dataset$CN, 1, run.prediction, scale = T,
                      model = model) %>% 
      t %>% rescale.denoised.df(., dataset$CN)
  } else if(method == "rollmed") {
    denoised <- pbapply(dataset$CN, 1, roll, k = k) %>% as.data.frame %>% t
  } else if(method == "bcp") {
    denoised <- pblapply(segs, \(s) s$segment) %>% do.call(cbind, .) %>%
      as.data.frame %>% t
    if (compute.segmentation) {
      warning("You have selected bcp and compute.segmentation. 
              Skipping segmentation step...")
    }
  } else {
    stop("Invalid method. Choose from 'rollmed', conv1d' and 'bcp'")
  }
  rownames(denoised) <- rownames(dataset$CN)
  # check if any prediction is constant
  is.constant <- apply(denoised, 1, \(row) sum(diff(row) != 0) == 0)
  n.constant <- sum(is.constant)
  na.df <- data.frame(
    purity = rep(NA, n.constant), 
    mse = rep(NA, n.constant),
    est.purity = rep(NA, n.constant)
  )
  if (all(is.constant)) {
    warning("No non constant predictions. Skipping..") 
    na.df
  } else {
    if (n.constant > 0){
      warning(sprintf("Found %i constant predictions.", n.constant))
    }
    # remove constant predictions
    denoised <- denoised[!is.constant, ]
    print("[INFO] Estimating purities...")
    # purity interval and adjustments
    p <- seq(0.005, 1, length.out = 200) 
    # compute estimated purities 
    pure <- pblapply(1:nrow, \(i) compute.purity(denoised[i,], 
                                                 name = rownames(denoised[i,]), 
                                                 p = p, show.plots = F)) %>%
      do.call(cbind, .) %>% t
    rownames(pure) <- rownames(denoised)
    colnames(pure) <- c("mean", "median")
    pure.df <- data.frame(pure, true = dataset$purity %>% as.numeric)
    if(show.plots) {
      print(plot.purity(pure.df))
    }
    print("[INFO] Computing MSE...")
    # if we have chosen conv1d or rolling median, we should compute segments
    if ((method != "bcp") & compute.segmentation) {
      denoised %<>% as.matrix
      denoised <- pbapply(denoised, 1, \(sample) {
        compute.segmentation(sample, show.plot = FALSE) %>% get.profile
      }) %>% t
    }
    # correct profiles 
    corr.profs <- pblapply(1:nrow, \(i) correct.profile(denoised[i,], 
                                                        method = "none",
                                                        pure.df$median[i]))
                                                        #dataset$purity[i]))
    corr.segs <- pblapply(1:nrow, \(i) correct.profile(dataset$segment[i,], 
                                                       method = "none",
                                                       dataset$purity[i]))
    mses <- lapply(1:nrow, \(i) mean((corr.profs[[i]]-corr.segs[[i]])^2)) %>%
      unlist
    df <- data.frame(purity = dataset$purity %>% as.numeric, mse = mses, 
                     mean = pure.df$mean, median = pure.df$median)
    if(sum(is.na(df$mse)) > 0) {
      warning("Purity estimation produced NAs.", call. = FALSE)
    }
    if(show.plots) {
      df.plot <- aggregate(mse ~ purity, df, \(x) 
                           c(mean = mean(x), sd = sd(x))) %>% 
        do.call(data.frame, .)
      df.plot %<>% na.omit
      print(plot.mse(df.plot))
    }
    rbind(df, na.df)
  }
}

# ----- Pipeline for segmentation ----------------------------------------------

# compare breakpoints between original and denoised data
# returns a list with the number of true positives, false positives and false
# negatives
# orig.bp: true breakpoints
# den.bp: predicted breakpoints
# margin: allowed margin of error
.compare.bp <- \(orig.bp, den.bp, margin = 5) {
  tp <- 0; fp <- 0; fn <- 0
  for (i in 1:length(orig.bp)) {
    if (any(abs(den.bp - orig.bp[i]) <= margin)) {
      tp <- tp + 1
    } else {
      fn <- fn + 1
    }
  }
  fp <- length(den.bp) - tp
  return(list(tp = tp, fp = fp, fn = fn))
}

# from original and estimated breakpoints, compute precision and recall
.get.metrics <- \(orig.bp, den.bp, margin = 5) {
  cm <- .compare.bp(orig.bp, den.bp, margin = margin)
  precision <- cm$tp / (cm$tp + cm$fp)
  recall <- cm$tp / (cm$tp + cm$fn)
  return(list(precision = precision, recall = recall))
}

# Pipeline for computing segmentation quality
# dataset: dataset to denoise
# denoised: pre-denoised data
# method: the method to use
# model: already trained denoising autoencoder
run.tuning.segmentation <- \(dataset, denoised = NULL, method = "conv1d",
                             model = NULL) {
  # if not null, we have already performed segmentation etc.
  if (is.null(denoised)) {
    if(method == "conv1d"){
      if (is.null(model)) {
        model <- run.pipeline(train, model = run.conv1d.denoiser, 
                              penalty = lambda)  
        model <- model$model
      }
      denoised <- run.prediction(model, dataset$CN)[,,1]
    } else if(method == "bcp") {
      segs <- pbapply(dataset$CN %>% t, 2, compute.segmentation.bcp)
      denoised <- pblapply(segs, post.process.bcp, show.plots = F) %>% 
        lapply(., get.profile) %>% do.call(rbind, .) %>% as.data.frame
    } else {
      stop("Method not recognized. Valid methods are 'conv1d', 'bcp'. Stopping..")
    }
  }
  # check if any prediction is constant
  is.constant <- apply(denoised, 1, \(row) sum(diff(row) != 0) == 0)
  n.constant <- sum(is.constant)
  na.df <- data.frame(
    value = rep(0, 2*n.constant), 
    metric = rep(c("precision", "recall"), each = n.constant),
    purity = rep(dataset$purity, times = 2*n.constant)
  )
  if (all(is.constant)) {
    warning("No non constant predictions. Skipping..") 
    na.df
  } else {
    if (n.constant > 0){
    warning(sprintf("Found %i constant predictions.", n.constant))
    }
    # remove constant predictions
    denoised <- denoised[!is.constant, ]
    # perform segmentation and compare with the original
    if(method == "conv1d"){
      segments <- apply(denoised, 1, \(s) 
                        compute.segmentation(s, show.plot = F) %>%
                          postprocess.segment)
    } else if(method == "bcp") {
      segments <- apply(denoised, 1, segment.df)
    }
    n <- nrow(dataset$segment)
    metrics <- lapply(1:n, \(i) 
                      .get.metrics(dataset$segment[i,] %>% get.breakpoints %>% 
                                     which, 
                                  segments[[i]]$start[-1])) %>% transpose
    precision <- metrics$precision %>% unlist
    recall <- metrics$recall %>% unlist 
    df.metrics <- data.frame(
      value = c(precision, recall), 
      metric = rep(c("precision", "recall"), 
                   times = c(length(precision), length(recall))),
      purity = rep(dataset$purity, times = 2)
    )
    rbind(df.metrics, na.df)
  }
}

# ----- Function for running the tuning pipeline -------------------------------

# Run pipeline for several lambdas
# models: list of models to use
# data: dataset to tune on
# lambdas: vector of penalty values used in the tuning
# compute.seg: TRUE if segmentation should be performed, otherwise just denoising
# pipeline.fun: the pipeline to run, either mse or seg
run.penalty.tuning <- \(models, data, lambdas = seq(0.0, 0.025, 0.05, 0.075, 1),
                        compute.seg = F, show.plot = FALSE,
                        pipeline.fun = "mse") {
  # make sure that valid option is provided
  if (!(pipeline.fun %in% c("mse", "seg"))) {
    stop(sprintf("Option %s not recognized. Choose from 
                 mse and seg. Stopping.."))
  }
  if ((length(models) != length(lambdas))) {
    stop("The number of models and lambdas must be the same. Stopping..")
  }
  if (pipeline.fun == "mse") {
    results <- 
      pblapply(models, \(m) 
               compute.mse.pipeline(data, method = "conv1d", model = m,
                                    compute.segmentation = compute.seg, 
                                    show.plots = F)) %>% 
      do.call(rbind, .)
    results$model.id <- rep(1:length(models), each = nrow(data$CN))
    results$penalty <- rep(lambdas, each = nrow(data$CN))
    if (show.plot) {
      print(plot.penalty.tuning(results))
    }
  }
  else { # at this stage, only the option "seg" is feasible
    results <- 
      pblapply(models, \(m) 
               run.tuning.segmentation(data, model = m)) %>% do.call(rbind, .)
    results$model.id <- rep(1:length(models), each = 2*nrow(data$CN))
    results$penalty <- rep(lambdas, each = 2*nrow(data$CN))
    if (show.plot) {
      print(plot.penalty.seg(results %>% remove.missing))
    }
  }
  results
}

# ----- Plot functions ---------------------------------------------------------

# this only makes sense if we use a small number of purity values
plot.penalty.tuning <- \(df, n.groups = 4, log = F) {
  df %<>% na.omit
  df$penalty %<>% as.factor
  df$model.id %<>% as.factor
  df$purity.group <- df$purity %>% cut(n.groups)
  if(log) df$mse %<>% log
  df %>% ggplot(aes(x = model.id, y = mse, group = model.id,
                    color = penalty)) + geom_boxplot() + 
    facet_grid(purity.group ~ ., scales = "free_y") +
    theme(legend.position = "top") + 
    labs(y = "MSE of purity corrected reconstruction", x = "Model id")
}


plot.penalty.seg <- \(df, n.groups = 5) {
  df$penalty %<>% as.factor
  df$model.id %<>% as.factor
  df$purity.group <- df$purity %>% cut(n.groups)
  df %>% ggplot(aes(x = model.id, y = value, group = model.id,
                    color = penalty)) + geom_boxplot() + 
    facet_grid(purity.group ~ metric, scales = "free_y") +
    theme(legend.position = "top") + 
    labs(y = "Metric", x = "Model id")
}

# ----- Plots of autoencoder tuning --------------------------------------------

# custom plot function for the tuning results
plot.models.custom <- \(df.tuned, ylim = NULL, remove.outlier = T, 
                        lambda = 0.1) {
  # subset the dataframe to only include on penalty value
  df.tuned %<>% {.[10 * .$penalty == 10 * lambda, ]}
  kernels <- unique(df.tuned$kernel.size)
  feature.size <- unique(df.tuned$feature.size)
  plot.list <- lapply(
    kernels,
    \(kernel) {
      df.temp <- df.tuned[df.tuned$kernel.size == kernel, ]
      df.temp.large <- df.temp[df.temp$feature.size != max(feature.size), ]
      df.temp.small <- df.temp[df.temp$feature.size == max(feature.size), ]
      pl.large <- df.temp.large %>% 
        ggplot(aes(x = layers, y = val.loss, group = layers)) + 
        geom_boxplot(position = position_dodge(.9)) + 
        labs(color = NULL, y = NULL, x = NULL) +
        stat_summary(fun = "mean", geom = "point", shape = 8, 
                     size = 3, color = "red",
                     position = position_dodge(.9)) + 
        facet_wrap(~ feature.size, scales = "free_y") + 
        scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
      if (kernel != min(kernels)) {
        pl.large <- pl.large + theme(
          strip.background = element_blank(),
          strip.text.x = element_blank()
        )
      }
      if (!is.null(ylim)) {
        pl.large <- pl.large + coord_cartesian(ylim = ylim)
      }
      pl.small <- df.temp.small %>% 
        ggplot(aes(x = layers, y = val.loss, group = layers)) + 
        geom_boxplot(position = position_dodge(.9))  + 
        stat_summary(fun = "mean", geom = "point", shape = 8, 
                     size = 3, color = "red",
                     position = position_dodge(.9)) +
        labs(x = NULL, y = NULL) + 
        scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
      if (kernel == min(kernels)) {
        pl.small <- pl.small + facet_grid(kernel.size ~ feature.size)
      }
      else {
        pl.small <- pl.small + facet_grid(kernel.size ~ .)
      }
      if (!is.null(ylim)) {
        pl.small <- pl.small + coord_cartesian(ylim = ylim)
      }
      plot_grid(pl.large, pl.small, rel_widths = c(0.65, 0.35))
    }
  )
  helper.plot <- ggplot() + geom_line() + 
    labs(x = "# layers", y = "validation loss")
  plot <- plot_grid(plotlist = plot.list, ncol = 1,
                    rel_heights = c(1.18, 1, 1, 1))
  y.grob <- textGrob("validation loss", gp = gpar(fontsize = 14), rot = 90)
  x.grob <- textGrob("# layers", gp = gpar(fontsize = 14))
  title <- textGrob(paste0("Penalty value: ", lambda), gp = gpar(fontsize = 16))
  grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob), top = title)
}

# ----- Plots for rolling median tuning ----------------------------------------

plot.roll.tuning <- \(df, n.groups = 5, log = T) {
  df$k %<>% as.factor
  if(log) df$mse %<>% log
  df$purity.group <- df$purity %>% cut(n.groups)
  ggplot(df, aes(x = k, y = mse, group = k, color = k)) + geom_boxplot() +
    theme(legend.position = "top") + 
    labs(y = "MSE of purity adjusted profiles") +
    facet_grid(noise ~ purity.group, scales = "free")
}

