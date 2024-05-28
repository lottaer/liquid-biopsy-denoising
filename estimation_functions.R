#     Estimation of purity and subclonal ratio from liquid biosy sequencing. 
#     Inspired by Lakatos et al., <https://github.com/elakatos/liquidCNA>.
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
library(gtools) # needed for permutations

# ----- Purity estimation ------------------------------------------------------
## ----- Help functions --------------------------------------------------------
# plot the peaks of a density estimate
.plot.peaks <- \(density.est, modes, name = NULL) {
  data.frame(
    x = density.est$x,
    y = density.est$y
  ) %>%
    ggplot(aes(x = x, y = y)) + 
    geom_area(fill = "lightgrey") + geom_line() +
    geom_vline(xintercept = modes, linetype = "dashed", color = "firebrick") +
    labs(y = NULL, x = NULL) + 
    ggtitle(name)
}

# plot the peaks of a density estimate in a grid
.plot.peaks.grid <- \(fits.list, name) {
  pl.a <- plot_grid(plotlist = fits.list$plots, ncol = 2, align = "vh")
  title <- ggdraw() + draw_label(name)
  pl.a <- plot_grid(title, pl.a, ncol = 1,
                    rel_heights=c(0.1, 1))
  pl.a
}

# plot the probabilities of purities for different adjust values
.plot.purity.probs <- \(p.fits, p, mins, adj, name) {
  p.fits$p <- p
  p.fits <- melt(p.fits,id='p')
  p.fits$adj <- rep(adj, each = length(p))
  pl.f <- ggplot(p.fits, aes(x=p, y=value, color=factor(adj)))
  if(!is.null(name)) {
    pl.f <- pl-f + geom_point(size=2, alpha=0.3, label = sprintf("%.2f"))
  } else {
    pl.f <- pl.f + geom_point(size=2, alpha=0.3)
  }
  pl.f <- pl.f + theme_bw() +
    geom_vline(xintercept = mins, linetype=c('dashed','dotted'), 
               colour=c('red','blue')) +
    labs(title = name, color = "adjust") +
    theme(legend.position = "top")
  pl.f
}

# compute the modes of a density estimate
.compute.modes <- \(density, epsilon = 0.1) {
  cutoff <- density$y %>% max * epsilon
  density$x %<>% .[which(density$y>(0 + cutoff))]
  density$y %<>% .[which(density$y>(0 + cutoff))]
  density$x[which((diff(sign(diff(density$y))) == -2))]
}

# adjust a density s.t. the highest peak is 2
.density.adjusted <- \(denoised, adjust, lower, upper, epsilon = 0) {
  density.est <- density(denoised, from = lower, to = upper, adj = adjust)
  modes <- .compute.modes(density.est, epsilon = epsilon)
  modes.max <- which.max(density.est$y) # find index of highest peak
  max.x <- density.est$x[modes.max]
  mode.2 <- which.min(abs(modes-max.x))
  modes.y <- density.est$y[which(density.est$x %in% modes)]
  modes.adj <- modes / modes[mode.2] * 2
  density.est$x <- density.est$x / modes[mode.2] * 2
  return(list(density.est = density.est, modes.adj = modes.adj, 
              modes.y = modes.y))
}

# evaluate the purity density
.eval.purity.density <- \(p, modes, modes.y) {
  if((length(modes) == 1) | (length(modes) >= 6)) {
    return(NA)
  } 
  mode.2 <- which(modes == 2)
  nb.modes <- ((2-mode.2+1):(length(modes)+2-mode.2))[-mode.2]
  modes <- modes[-mode.2]
  modes.y <- modes.y[-mode.2]
  sapply(1:(length(modes)), \(i) 
         min(((2+p*(nb.modes[i]-2))-modes[i])^2) * modes.y[i]^2) %>% 
    sum %>% sqrt
}

## ----- Estimation of purity --------------------------------------------------
# compute the purity of a sample
# denoised: the denoised data
# r: the range of purity values to be tested
# n.adj: the number of adjustments to be tested
# adj.min, adj.max: the minimum and maximum adjustment value
# adj: the adjustment values to be tested, if not provided, will be computed
#      according to n.adj, adj.min and adj.max
# show.plots: whether to show the plots or not
# grid.plots: how many plots to show in the grid
# name: the name of the plot
compute.purity <-
  \(denoised, p = seq(0.005, 1, length.out = 200),
    n.adj = 10, adj.min = 0.8, adj.max = 1.6, 
    adj = seq(adj.min, adj.max, length.out = n.adj), 
    show.plots = FALSE, grid.plots = 1:(min(4,n.adj)), name = NULL, 
    epsilon = 0.1, n.sd = 2.5) {
    denoised %<>% as.numeric
    # define the region to keep (otherwise too many modes introduced)
    lower <- denoised %>% {mean(.) - n.sd * sd(.)}
    upper <- denoised %>% {mean(.) + n.sd * sd(.)}
    fits.list <- lapply(adj , \(a) {
      d.adj <- .density.adjusted(denoised, a, lower, upper, epsilon = epsilon)
      pl <- .plot.peaks(d.adj$density.est, d.adj$modes.adj) + 
        ggtitle(sprintf("Adjust = %.1f", a)) + 
        theme(legend.position = "none", plot.title = element_text(size = 10))
      p.fits <- sapply(p, \(p) .eval.purity.density(p, d.adj$modes.adj, 
                                                    d.adj$modes.y))
      list(p.fits = p.fits, plots = pl)
    }) %>% purrr::transpose(.)
    p.fits <- fits.list$p.fits %>% as.data.frame
    p.cols <- p.fits[ ,1:length(adj)] %>% as.matrix
    which.mins <- apply(p.cols, 2, which.min) %>% unlist
    mins <- c(mean(p[which.mins], na.rm = TRUE),
              median(p[which.mins], na.rm = TRUE))
    # mins <- c(mean(p[apply(p.cols, 2, which.min)], na.rm=T),
    #           median(p[apply(p.cols, 2, which.min)], na.rm=T))
    names(mins) <- c("mean", "median")
    if(show.plots){
      .plot.peaks.grid(fits.list[grid.plots], name) %>% print
      .plot.purity.probs(p.fits, p, mins, adj, name) %>% print
    }
    if(length(mins) == 0) {
      mins <- c(NA, NA)
    }
    return(mins %>% as.data.frame)
  }

# compute the purity of a set of insilico samples
# data: the data matrix, rows are samples, columns are bins
# metadata: the metadata data frame, containing the purity values and the number
#           of reads
# method: the denoising method to be used
#         bcp: Bayesian change point detection
#         nn: neural network, denoising autoencoder
#         rollmed: rolling median
#         none: no denoising
# model: the denoising model to be used, if method is "nn" 
# k: the window size for the rolling median, if method is "rollmed"
run.purity.insilico <- \(data, metadata, method = "none", model = NULL, k = 21) {
  if (!(method %in% c("nn", "rollmed", "bcp", "none"))) {
    stop("Choose from methods {nn, rollmed, bcp, none}")
  }
  print(sprintf("[INFO] Running denoising for method %s", method))
  if (method == "bcp") {
    denoised <- pbapply(data, 1, \(row) bcp(row)$posterior.mean) %>% t
  }
  else if (method == "nn") {
    if (is.null(model)) {
      model <- load.denoiser("denoising_model.keras")
    }
    denoised <- run.prediction(model, data)[,,1]
  }
  else if (method == "rollmed") {
    denoised <- apply(data, 1, \(row) roll(row, k = k)) %>% t
  }
  else {
    denoised <- data # do nothing
  }
  print(sprintf("[INFO] Estimating purities %s", method))
  pblapply(1:nrow(denoised), \(i) {
    data.frame(
      purity = compute.purity(denoised[i, ])[2, ], # median
      true.purity = metadata$purity[i],
      reads = metadata$n.reads[i],
      method = method
    )
  }) %>% do.call(rbind, .)
}

## ----- Plotting --------------------------------------------------------------
# plot the results of the purity estimation
# purity.df: the data frame containing the results from compute.purity
# replicates: the number of replicates (if NULL, does nothing)
# mean: whether to plot the mean or the median
plot.purity <- \(purity.df, replicates = NULL, mean = TRUE) {
  if(is.numeric(replicates)) {
    if(length(replicates) == 1) {
      replicates <- rep(rep(1:5, each = replicates), 
                        length.out = nrow(purity.df))
    }
  } else if(is.null(replicates)) {
    replicates <- numeric(nrow(purity.df))
  } else {
    stop("'replicates' has to be numeric")
  }
  purity.df$replicate <- factor(replicates)
  if(mean) {
    p <- ggplot(purity.df, aes(x = true, y = mean, group = replicate))
  } else {
    p <- ggplot(purity.df, aes(x = true, y = median, group = replicate))
  }
  p <- p +
    geom_point(aes(
      shape = replicate,
      color = replicate
    ), 
    size = 3, alpha = 0.7) + 
    geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", 
                # color = "black", 
    ) +  
    geom_smooth(aes(group = NULL), method = "lm", se = FALSE, 
                color = "grey",  
                formula = "y~x") +
    theme_light() +
    labs(#title = "Estimated vs. True purity",
      x = "True purity",
      y = "Estimated purity") +
    # theme(legend.position = "none", axis.title.x=element_blank(), 
    #       axis.title.y=element_blank()) +
    scale_shape_manual(values = c(16, 15, 17:19, 8))
  if(length(unique(replicates)) == 1) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

# plot the results of the purity estimation for multiple sigmas
plot.multiple.purities <- \(purity.df, mean = TRUE) {
  purity.df$sigma %<>% as.factor
  if(mean){
    p <- ggplot(purity.df, aes(x = true, y = mean, color = sigma))
  } else {
    p <- ggplot(purity.df, aes(x = true, y = median, color = sigma))
  }
  p <- p + geom_point() + facet_wrap(. ~ sigma, ncol = 2) + 
    geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
    theme(legend.position = "top") + 
    labs(x = "True purity", y = "Estimated purity")
  p
}

# plot the results of the purity estimation for in silico data
plot.purity.insilico <- \(df) {
  df$reads  %<>% as.factor
  ggplot(df, aes(true.purity, purity, color = reads, shape = reads)) + 
    geom_smooth(aes(shape = NULL, color = NULL), method = "lm", 
                formula = "y~x", se = F, color = "darkgray", linewidth = 0.5) +
    geom_point() +
    scale_shape_manual(values = c(19, 17, 15, 18)) + 
    theme(legend.position = "top") +
    labs(color = "# of reads [million]", shape = "# of reads [million]",
         x = "True purity", y = "Estimaded purity") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "darkgray") +
    facet_wrap(. ~ method)
}

# ----- Subclonal ratio estimation ---------------------------------------------
## ----- Segment functions -----------------------------------------------------

# computes the purity corrected CNA profiles using estimated purity
# profile: the CNA profile to be corrected, as a numeric vector
# purity: the estimated purity of the profile
# method: median   -- uses the median as cn 2
#         mode     -- uses the maximum mode as cn 2
#         none     -- assumes cn 2 is at level 1 in the profile
#         constant -- uses a given constant value as cn 2. if const = 1,
#                     'constant' is the same as 'none'
correct.profile <- \(profile, purity, method = "none", base = 2, const = 1) {
  profile %<>% as.numeric
  if(method == "median") {
    base / purity * (profile / median(profile) - 1) + base
  } else if(method == "mode"){
    dens <- density(profile)
    mode <- which.max(dens$y)
    base / purity * (profile / dens$x[mode] - 1) + base
  } else if(method == "none") {
    base / purity * (profile - 1) + base
  } else if(method == "constant"){
    base / purity * (profile / const - 1) + base
  } else {
    stop("Unknown method. Valid methods are 'median', 'mode', 'none', 'constant'.")
  }
}

# compute the union of the breakpoints from longitudinal data and remove the 
# ones that are shorter than min.bins
# segment.df.list: a list of dataframes with start and end columns
# min.bins: the minimum number of bins for a segment to be included
get.long.bp <- \(segment.df.list, min.bins = 12) {
  starts <- lapply(segment.df.list, \(x) x$start) %>% unlist %>% unique %>% sort
  ends <- lapply(segment.df.list, \(x) x$end) %>% unlist %>% unique %>% sort
  breakpoints <- data.frame(start = starts, end = ends)
  breakpoints %<>% .[((.$end - .$start) >= min.bins),]
  return(breakpoints)
}

## ----- Compute subclonal segments --------------------------------------------
# identify subclonal segments from set of segments
.identify.segments <- \(segments, cutoff = 0.05) {
  # find the segments that deviates enough from 0
  which(sapply(1:nrow(segments), \(i) {
    # might be better to just look at the largest difference and not all of them
    # all(abs(segments[i, -c(1)]) > cutoff)
    any(abs(segments[i, -c(1)]) > cutoff)
  }))
}

# identify unstable segments from set of subclonal segments
.identify.unstable <- \(segments, delta = 0.005) {
  # find monotonous segments
  which(sapply(1:nrow(segments), \(i){
    !(all(diff(segments[i, ]) > -delta) | 
        all(diff(segments[i, ]) < delta))
  }))
}

# plot the segments from compute.subcl.segs with color corresponding to their
# type -- clonal, subclonal or unstable
.plot.delta.segment <- \(delta.cn, subclonal, unstable = NULL) {
  delta.cn %<>% as.data.frame
  delta.cn$id <- 1:nrow(delta.cn)
  delta.cn.melt <- melt(delta.cn, "id")
  delta.cn.melt$sample.id <- rep(1:(ncol(delta.cn) - 1), each = nrow(delta.cn))
  sample.types <- rep("clonal", nrow(delta.cn))
  sample.types[subclonal] <- "subclonal"
  sample.types[unstable] <- "unstable"
  delta.cn.melt$segment.type <- rep(sample.types, times = ncol(delta.cn) - 1)
  ggplot(delta.cn.melt, 
         aes(x = sample.id, y = value, group = id, color = segment.type)) + 
    geom_point(aes(color = segment.type)) + geom_line() +
    theme(legend.position = "top") + labs(color = NULL) +
    scale_colour_manual(values=c('lightgray','firebrick3','dodgerblue4')) 
}

# computes which segments to include in the subclonal estimation
# data.corr: purity-corrected data, as a matrix or dataframe, 
#            with samples in rows and bins in columns
# bp.list: a list of breakpoints, with start and end columns
# cutoff: the cutoff value for the identification of subclonal segments
# baseline: the baseline value for the CNA profile (where is CN 2)
# show.plots: whether to show the plots or not
# delta: the delta value for the identification of unstable segments
# optim: whether to return the optimal cutoff value or not
compute.subcl.segs <- \(data.corr, bp.list, cutoff = 0.05, 
                        baseline = 1, show.plots = FALSE, delta = 0.075, 
                        optim = FALSE) {
  bps <- get.long.bp(bp.list)
  # compute the segment length (will be used to sort segments)
  bps$length <- bps %$% (end - start + 1)
  cn.segs <- apply(data.corr, 1, \(s) seg.from.bps(s, bps)$value)
  delta.cn.segs <- apply(cn.segs, 2, \(s) s - cn.segs[,1])
  # get subclonal (stable and unstable) segment indices
  subclonal <- .identify.segments(delta.cn.segs, cutoff = cutoff)
  if(length(subclonal) <= 1) {
    if(optim) {
      return(data.frame(n.stable = 0, n.unstable = 0, ratio = 0))
    } else {
      return(matrix(NA, ncol = 1, nrow = nrow(delta.cn.segs)))
    }
  }
  n.samples <- ncol(delta.cn.segs)
  perms <- permutations(n.samples-1, n.samples-1, 2:n.samples)
  # browser()
  unstable.list <- apply(perms, 1, \(p){
    order <- c(1, p)
    u <- .identify.unstable(delta.cn.segs[subclonal, order], delta = delta)
  })
  if(length(unstable.list) == 0) {
    unstable.list <- 0
  }
  unstable <- unstable.list %>% .[[which.min(lengths(.))]]
  subcl.unst <- subclonal[unstable]
  subcl.st <- setdiff(subclonal, subcl.unst)
  if(show.plots) {
    print(.plot.delta.segment(delta.cn.segs, subcl.st, subcl.unst))
  }
  if(optim){
    # returns a list used to identify the optimal cutoff value
    data.frame(n.stable = length(subcl.st), n.unstable = length(subcl.unst),
               ratio = length(subcl.st) / length(subclonal))
  } else {
    # return the values of the subclonal segments, including multiplicity
    subclonal.segments <- delta.cn.segs[subcl.st, ]
    # plot.profiles(apply(cn.segs, 2, \(s) rep(s, times=bps$length)) %>% t) %>% plot
    # apply(subclonal.segments, 2, \(s) rep(s, times=bps$length[subcl.st]))
    subclonal.segments
  }
}

# plots the optimal cutoff value
.plot.optim.cutoff <- \(df, opt) {
  p1 <- ggplot(df, aes(x = cutoff, y = ratio)) + geom_line() + 
    labs(y = "Ratio of stable segments") + 
    geom_vline(xintercept = opt, linetype = "dashed")
  p2 <- ggplot(df) +
    geom_line(aes(x = cutoff, y = n.stable, color = "Subclonal")) +
    geom_line(aes(x = cutoff, y = n.unstable, color = "Unstable")) +
    geom_vline(xintercept = opt, linetype = "dashed") +
    # theme(legend.position = "top") + 
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) +
    labs(color = NULL, y = "Number of segments")
  plot_grid(p1, p2)
}

# get the optimal cutoff value for the choice of subclonal segments
# data.corr: the purity corrected data, as a matrix or dataframe, with samples in 
#            rows and bins in columns
# bp.list: a list of breakpoints, with start and end columns
# cutoff.vec: a vector of cutoff values to be tested
# min.nseg, max.nseg: the minimum/maximum number of segments to be included
# baseline: the baseline value for the CNA profile (where is CN 2)
# show.plots: whether to show the plots or not
# delta: the delta value for the identification of unstable segments
optim.cutoff <- \(data.corr, bp.list, cutoff.vec, min.nseg = 5, max.nseg = 20, 
                  baseline = 1, show.plots = FALSE, delta = 0.075) {
  df <- lapply(cutoff.vec, \(cutoff) {
    compute.subcl.segs(data.corr, bp.list, cutoff = cutoff, baseline = baseline, 
                       show.plots = F, delta = delta, optim = T)
  }) %>% do.call(rbind, .)
  ind <- which((df$n.stable >= min.nseg) & (df$n.stable <= max.nseg))
  max.ind <- df$ratio[ind] %>% which.max
  opt <- cutoff.vec[ind][max.ind]
  if(show.plots) {
    df$cutoff <- cutoff.vec
    rbind(df) %>% .plot.optim.cutoff(opt = opt) %>% plot
  }
  opt
}

## ----- Compute subclonal ratio -----------------------------------------------
# evaluate the sybclonal ratio
.eval.subcl.ratio <- \(r, modes, which.modes = NULL) {
  if(length(modes) <= 1 | length(modes) >= 4) {
    return(NA)
  }
  if(is.null(which.modes)){
    pos.ix <- which.max(modes > 0)
    if(!any(modes > 0)) {
      pos.ix <- length(modes) + 1
    }
    which.modes <- (1:(length(modes)+1)-pos.ix)[-pos.ix]
  }
  if(max(which.modes) >=4 | min(which.modes) <= -3) {
    return(NA)
  }
  sapply(1:(length(which.modes)), \(i) {
    min(((r*(which.modes[i]))-modes[i])^2)
  } ) %>% 
    sum %>% sqrt
}

# computes the subclonal ratio based on a vector of subclonal segments
# subclonal: the subclonal segments, as a numeric vector
# which.modes: the indices of the modes to be used, if NULL, it is assumed the 
#              modes are evenly spaced
# r: the range of subclonal ratios to be tested
# n.adj: the number of adjustments to be tested
# adj.min, adj.max: the minimum and maximum adjustment value
# adj: the adjustment values to be tested, if not provided, will be computed
#      according to n.adj, adj.min and adj.max
# epsilon: the epsilon value for the identification of modes
# n.sd: the number of standard deviations to be used for the range of the
#       subclonal segments
# show.plots: whether to show the plots or not
# return.median: whether to return the median or the mean of the optimal ratios
compute.subcl.ratio <- \(subclonal, which.modes = NULL, 
                         r = seq(0.01, 0.6, length.out = 100), 
                         n.adj = 10, adj.min = 1, adj.max = 3, 
                         adj = seq(adj.min, adj.max, length.out = n.adj),
                         epsilon = 0.1, n.sd = 4, show.plots = FALSE, 
                         return.median = TRUE) {
  mean <- mean(subclonal); sd <- sd(subclonal)
  from <- (subclonal[subclonal >= (mean - n.sd * sd)] %>% min) - 0.1 * sd
  to <- (subclonal[subclonal <= (mean + n.sd * sd)] %>% max) + 0.1 * sd
  adj.list <- lapply(adj , \(a) {
    density <- density(subclonal, adjust = a, from = from, to = to)
    modes <- .compute.modes(density, epsilon = epsilon)
    r.fits <- sapply(r, \(r) .eval.subcl.ratio(r, modes, which.modes)) %>%
      as.numeric
    # plot the different smoothed densities for diagnostic reasons
    if (show.plots) {
      pl <- .plot.peaks(density, modes) + 
        ggtitle(sprintf("Adjust = %.1f", a)) + 
        theme(legend.position = "none", plot.title = element_text(size = 10))
      list(optimal.ratio = r[which.min(r.fits)], plots = pl)
    }
    else {
      list(optimal.ratio = r[which.min(r.fits)])
    }
  }) %>% transpose
  if (show.plots) {
    print(.plot.peaks.grid(adj.list, NULL))
  }
  ifelse(
    return.median, 
    median(adj.list$optimal.ratio %>% unlist, na.rm = T), 
    mean(adj.list$optimal.ratio %>% unlist, na.rm = T)
  )
}

# ----- Pipeline functions -----------------------------------------------------
# help function for preprocess.subcl
.pre.subcl <- \(data.cn, prob.cutoff = 0.1) {
  n <- nrow(data.cn)
  segs <- compute.segmentation.bcp(data.cn %>% t)
  denoised <- segs$segment %>% as.data.frame %>% t
  # compute estimated purities
  purities <- lapply(1:n, \(i) compute.purity(denoised[i,] %>% as.numeric, 
                                              show.plots = F)) %>%
    do.call(cbind, .) %>% t %>% as.data.frame
  bps <- which(segs$probs > prob.cutoff)
  bp.list <- data.frame(start = c(1, bps+1), end = c(bps, nrow(segs$segment)),
                        value = rep(0, length(bps) + 1))
  # bp.list <- postprocess.segment(bp.list)[,-3]
  # bp.list <- .filter.bps(bp.list)[,-3]
  bp.list <- bp.list[,-3]
  list(bp.list = bp.list, purity = purities$median)
}

# preprocess the data for subclonal tracking
# data: the data to be preprocessed, as a list with the following elements:
#       CN: the CNA profiles, as a matrix or dataframe, with samples in rows and
#           bins in columns
#       purity: the estimated purities of the samples
# prob.cutoff: the probability cutoff value for the breakpoints
# insilico: whether the data is in silico or not
# n: the number of samples to be denoised together in the insilico case. has to
#    be a divisor of the number of samples
preprocess.subcl <- \(data, prob.cutoff = 0.1, insilico = F, n = 2) {
  N <- length(data$purity)
  if(insilico) {
    if(N %% n != 0) {
      stop("The number of samples must be divisible by n")
    }
    pre <- lapply(1:(N/n), \(i){
      samples <- (i*n-1):(i*n) #    and here 
      .pre.subcl(data$CN[samples,], prob.cutoff) 
    }) %>% transpose
    pre$purity %<>% unlist
  } else {
    pre <- .pre.subcl(data$CN, prob.cutoff)
    pre$bp.list <- list(pre$bp.list)
  }
  return(pre)
}

# function for computing the subclonal ratio
# data.cn: the CNA profiles, as a matrix or dataframe, with samples in rows and
#          bins in columns
# bp.list: a list of breakpoints, with start and end columns
# r: the range of subclonal ratios to be tested
# cutoff.vec: a vector of cutoff values to be tested
# min.nseg: the minimum number of segments to be included
# adj: the adjustment values to be tested
# show.plots: whether to show the plots or not
# epsilon: the epsilon value for the identification of modes
# baseline: the baseline value for the CNA profile
subcl.tracking <- \(data.cn, bp.list, r = seq(0.01, 1, 0.01), cutoff.vec = r,
                    min.nseg = 5, adj = c(1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3), 
                    show.plots = 0, epsilon = 0.1, baseline = 1) {
  # note that the dataset has to be longitudinal, denoised and purity corrected
  cutoff <- optim.cutoff(data.cn, bp.list, cutoff.vec, min.nseg = min.nseg, 
                         baseline = baseline, show.plots = show.plots)
  subcl.segs <- compute.subcl.segs(data.cn, bp.list, cutoff = cutoff, 
                                   baseline = baseline, 
                                   show.plots = show.plots)
  if(length(subcl.segs) == 0 | length(cutoff) == 0) {
    return(rep(NA, nrow(data.cn)-1))
  } else {
    ifelse(show.plots == 2, show.plots <- T, show.plots <- F)
    est <- apply(subcl.segs[,-1] %>% as.matrix, 2, compute.subcl.ratio, r = r, 
                 adj = adj, epsilon = epsilon, show.plots = show.plots) %>% 
      as.data.frame %>% t
    return(est)
  }
}

# pipeline for subclonal tracking with bcp, inluding denoising and 
# purity estimation. can be combined with preprocess.subcl to save time
# data.list: a list of data to be analyzed, with the following elements:
#            CN: the CNA profiles, as a matrix or dataframe, with samples in rows
#                and bins in columns
#            purity: the estimated purities of the samples
#            ratio: the true subclonal ratios of the samples
# pre.list: a list of preprocessed data, as returned by preprocess.subcl.
#           if NULL, the data will be preprocessed automatically
# min.purity: the minimum purity value for the inclusion of samples
# cutoff.opt: whether to optimize the cutoff value or not
# show.plots: whether to show the plots or not
# insilico: whether the data is in silico or not
# adj: the adjustment values to be tested
subcl.tracking.pipeline <- \(data.list, pre.list = NULL, min.purity = 0.1, 
                             cutoff.opt = F, show.plots = F, insilico = F,
                             adj = c(seq(0.5, 2.5, 0.25))) {
  pblapply(1:length(data.list), \(i) {
    data <- data.list[[i]]
    N <- nrow(data$CN)
    if(is.null(pre.list)){
      pre <- preprocess.subcl(data, insilico = insilico)
    } else {
      pre <- pre.list[[i]]
    }
    purity <- pre$purity
    ifelse(insilico, pur.corr <- purity, pur.corr <- data$purity)
    range <- c(1, which(pur.corr[-1] > min.purity)+1)
    bp.list <- pre$bp.list
    ratios <- c(0, rep(NA, N-1))
    data.corr <- lapply(range, \(i) {
      corr.profs <- correct.profile(data$CN[i,], pur.corr[i], method = "none")
    }) %>% do.call(rbind, .)
    if(cutoff.opt){
      cutoff.vec <- seq(max(0.01, max(data$ratio[range]) - 0.20), 
                        min(1, max(data$ratio[range]) + 0.10), 0.01)
      ratios[range[-1]] <- subcl.tracking(data.corr, bp.list, 
                                          cutoff.vec = cutoff.vec, min.nseg = 5, 
                                          show.plots = show.plots, adj = adj)
    } else {
      ratios[range[-1]] <- subcl.tracking(data.corr, bp.list, min.nseg = 7, 
                                          show.plots = show.plots, adj = adj)
    }
    if(length(range) < 4) {
      ratios <- rep(NA, N)
    }
    data.frame(purity = purity,
               true.purity = data$purity,
               ratio = ratios,
               true.ratio = data$ratio,
               time = paste0("X", 1:N),
               id = rep(i, N))
  }) %>% do.call(rbind, .)
}

# ----- Plots ------------------------------------------------------------------
# plots the results of the subclonal tracking pipeline
# plot.df: the output of the subclonal tracking pipeline
plot.ratio <- \(plot.df) {
  ggplot(plot.df, aes(x = true.ratio, y = ratio)) + 
    geom_point(alpha = 0.3, size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(title = NULL,
         x = "True ratio",
         y = "Estimated ratio")
}

# plot the results of the subclonal tracking pipeline for insilico data
# plot.df: the output of the subclonal tracking pipeline
plot.ratio.insilico <- \(plot.df) {
  pl <- ggplot(plot.df, aes(resistant, resistant.est, color = reads, shape = reads)) + 
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "darkgray") +
    facet_wrap(. ~ method)
  if (!(sum(df$reads == "NA")) == nrow(plot.df)) {
    pl <- pl + scale_shape_manual(values = c(19, 17, 15, 18)) +
      theme(legend.position = "top") +
      labs(color = "\\# of reads [million]", shape = "\\# of reads [million]",
           x = "True subclonal-ratio", y = "Estimated subclonal-ratio")
  }
  else {
    pl <- pl + labs(x = "True subclonal-ratio", 
                    y = "Estimated subclonal-ratio") +
      theme(legend.position = "none")
  }
  pl
}