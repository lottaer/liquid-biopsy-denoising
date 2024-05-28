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

#  ----- Load packages ---------------------------------------------------------
library(purrr)
library(sn)

# Define default alteration and subclonal probabilities
.cna.probs <- c(2,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,4,4,5) %>% 
  {table(.) / length(.)} 
.subcl.probs <- c(rep(0, 18),-1,-1,1,1,1,2) %>%
  {table(.) / length(.)}

# ----- Helper functions -------------------------------------------------------

# Simulates a "true" sample (a copy number profile)
.simulate.cnv <- \(n.segments, min = 1, max = 5, 
                  cna.probs = .cna.probs,
                  min.length = 12, max.length = 130, genome = NULL,
                  step = NULL, normalize.sample = \(s) s / 2) {
  if (is.null(genome)) {
    genome <- rep(1, n.segments) # the baseline is 2 in diploid organism
  }
  position <- 1
  while (position < n.segments) {
    state <- sample(min:max, prob = cna.probs, size = 1)
    alteration.length <- sample(min.length:max.length, size = 1)
    end.position <- min(position + alteration.length, n.segments)
    genome[position:end.position] <- state
    position <- end.position
  }
  normalize.sample(genome)
}

# Simulate segment data and add random noise 
# N.samples: number of samples to simulate of each provided purity
# N.segments: number of genomic bins
# min: minimim copy number
# max: maximum copy number
# cancer.prop: purity of sample
# noise.mean: mean of Gaussian noise
# noise.sd: standard deviation of Gaussian noise
# simulated.data: data to alter
# cna.dep.noise: makes noise dependent on copy numbers
# min.length: minimal length of alterations
# max.length: maximal length of alterations
# downsample: simulates binning 500kb bins into 50kb bins
.simulate.noisy <- \(N.samples, N.segments, min = 1, max = 5,
                    cancer.prop = 1, noise.mean = 0.2, 
                    noise.sd = 0.005, simulated.data = NULL,
                    cna.dep.noise = FALSE,
                    min.length = 12, max.length = 130, downsample = FALSE) {
  if (is.null(simulated.data)) {
    simulated.data <- .simulate.cnv(N.segments, max = max, min = min,
                                    min.length = min.length, 
                                    max.length = max.length) 
  }
  # reduce the signal (normal contamination)
  segment.data <- sapply(1:N.samples, \(i) simulated.data) * cancer.prop  + 
                  (1 - cancer.prop) * rep(1, N.samples) 
  if(downsample) {
    segment.data.down <- segment.data
    segment.data %<>% apply(2, rep, each = 10)
  }
  noisy <- matrix(0, nrow(segment.data), ncol(segment.data))
  if (cna.dep.noise) {
    for (i in 1:N.samples) {
      sigma0 <- max(0.01, rnorm(1, noise.mean, noise.sd)) # use a minimal noise
      sigma <- sigma0 + 0.075*sigma0*segment.data[,i]
      noisy[, i] <- segment.data[, i] + rnorm(nrow(segment.data), 0, sigma)
    }
  }
  else {
    # create noisy measurements
    for (i in 1:N.samples) {
      sigma <- max(0.05, rnorm(1, noise.mean, noise.sd)) # use a minimal noise
      noisy[, i] <- segment.data[, i] + rnorm(nrow(segment.data), 0, sigma)
    }  
  }
  if(downsample) {
    noisy <- sapply(1:N.samples, \(s) sapply(1:N.segments, \(i) 
                                             mean(noisy[(10*i-9):(10*i),s])))
    segment.data <- segment.data.down
  }
  list(segment = segment.data, CN = noisy)
}

# Helper function, simulates dataset of one purity 
.simulate.noisy.dataset <- \(N.samples, N.segments, N.replicates = 2, 
                            max = 5, min = 1, cancer.prop = 1,
                            noise.mean = 0.2, noise.sd = 0.005,
                            simulated.data = NULL, cna.dep.noise = FALSE,
                            min.length = 12, max.length = 130, 
                            downsample = FALSE) {
  lapply(1:N.samples, \(i) .simulate.noisy(N.samples = N.replicates, 
                                           N.segments = N.segments, 
                                           max = max, min = min, 
                                           cancer.prop = cancer.prop, 
                                           noise.mean = noise.mean, 
                                           noise.sd = noise.sd, 
                                           simulated.data = simulated.data,
                                           cna.dep.noise = cna.dep.noise,
                                           min.length = min.length, 
                                           max.length = max.length, 
                                           downsample = downsample)) %>%
    {list(
      CN = lapply(., \(list) data.frame(list$CN) %>% t) %>% 
        do.call(rbind, .) %>% t,
      segment = lapply(., \(list) data.frame(list$segment) %>% t) %>% 
        do.call(rbind, .) %>% t
    )}
}

# ----- Simulate synthetic data -----------------------------------------------

# Simulates dataset of multiple purities and different noise levels
# N.samples: number of samples to simulate of each provided purity an noise level
# N.segments: number of genomic bins
# min: minimim copy number
# max: maximum copy number
# cancer.prop: vector of purities
# noise.mean: vector of means of Gaussian noise
# noise.sd: vector of standard deviation of Gaussian noise
# simulated.data: data to alter
# cna.dep.noise: makes noise dependent on copy numbers
# min.length: minimal length of alterations
# max.length: maximal length of alterations
# downsample: simulates binning 500kb bins into 50kb bins
simulate.dataset <- \(N.each, N.replicates, cancer.prop,
                      N.segments = 1000, min.CN = 1, max.CN = 5,
                      simulated.data = NULL, noise.mean = 0.2,
                      noise.sd = 0.005, cna.dep.noise = FALSE,
                      min.length = 12, max.length = 130, downsample = FALSE) {
  # the dataset will be of size N.each * length(props)
  if (length(cancer.prop) == 1) {
    cancer.prop %<>% c # add to vector to be able to apply over it
  }
  if (length(noise.mean) != length(noise.sd)) {
    stop("noise.mean and noise.sd must be of equal length")
  }
  if (length(noise.mean) == 1) {
    noise.mean %<>% c
    noise.sd %<>% c
  }
  if (length(noise.mean) != length(noise.sd)) {
    stop()
  }
  lapply(1:length(noise.mean), \(i) {
    lapply(cancer.prop,
           \(prop) {
             dataset <- .simulate.noisy.dataset(N.each, N.segments, N.replicates,
                                                cancer.prop = prop,
                                                noise.mean = noise.mean[i],
                                                noise.sd = noise.sd[i],
                                                cna.dep.noise = cna.dep.noise,
                                                min.length = min.length,
                                                max.length = max.length, 
                                                min = min.CN, max = max.CN, 
                                                downsample = downsample)
             list(CN = dataset$CN, segment = dataset$segment)
           }
    ) %>% purrr::transpose(.) %>%
      {list(
        CN = .$CN %>% do.call(cbind, .) %>% t,
        segment = .$segment %>% do.call(cbind, .) %>% t,
        purity = rep(cancer.prop, each = N.each * N.replicates),
        noise.mean = rep(
          noise.mean[i], each = N.each * N.replicates * length(cancer.prop)),
        noise.sd = rep(
          noise.sd[i], each = N.each * N.replicates * length(cancer.prop))
      )}
  }) %>% purrr::transpose(.) %>%
    {list(
      CN = .$CN %>% do.call(rbind, .),
      segment = .$segment %>% do.call(rbind, .),
      purity = .$purity %>% unlist, # concatenate the purities
      noise.mean = .$noise.mean %>% unlist,
      noise.sd = .$noise.sd %>% unlist
    )}
}

# a method to sample parameters for longitudinal data, for testing purposes
sample.params.long <- \(N, p = seq(12, 46, by = 2), p0 = seq(15, 46, by = 2), 
                        r = seq(5, 80, by = 5)) {
  data.frame(
    time = paste0("X", 1:N),
    purity = c(sample(p0, 1), sample(p, N-1, replace = T)) / 100,
    ratio = c(0, sample(seq(5, 80, by = 5), N-1, replace = T)) / 100
  )
}

# Simulate longitudinal data where the CN-profile is the same but purity
# and subclonal ratios are different
# N.segments: number of genomic bins
# purities: vector of purities to use
# ratios: vector of rations to use
# noise.mean: mean of Gaussian noise
# noise.df: standard deviation of Gaussian noise
# cna.dep.noise: makes noise dependent on copy numbers
# downsample: simulates binning 500kb bins into 50kb bins
simulate.long <- \(N.segments, purities, ratios, noise.mean = 0.2, 
                   noise.sd = 0.005, min.length = 12, max.length = 130,
                   cna.dep.noise = FALSE, downsample = F) {
  if(length(purities) != length(ratios)) {
    stop("The length of purities and ratios must be the same")
  }
  if(length(purities) == 1) {
    purities <- c(purities)
    ratios <- c(ratios)
  }
  norm.fun <- \(s) s / 2
  cn.sen <- .simulate.cnv(N.segments, cna.probs = .cna.probs, 
                          min.length = min.length, max.length = max.length, 
                          normalize.sample = norm.fun)
  seg.sen <- segment.df(cn.sen)
  change <- sample((-1):2, prob = .subcl.probs, size = nrow(seg.sen), 
                   replace = TRUE)
  seg.res <- seg.sen %>% mutate(value = pmax(1, value + norm.fun(change)))
  cn.res <- seg.res %>% get.profile
  lapply(1:length(ratios), \(i) {
    cn.data <- (1-ratios[i]) * cn.sen + ratios[i] * cn.res
    .simulate.noisy(N.samples = 1, N.segments = N.segments, max = 5, min = 1, 
                    cancer.prop = purities[i], noise.mean = noise.mean, 
                    noise.sd = noise.sd, simulated.data = cn.data,
                    cna.dep.noise = cna.dep.noise, downsample = downsample)
  }) %>% transpose %>%
    {list(
      CN = .$CN %>% do.call(cbind, .) %>% t,
      segment = .$segment %>% do.call(cbind, .) %>% t,
      purity = purities,
      ratio = ratios
    )}
}

sample.dataset.insilico <- \(data, parameters, n.datasets = 10, noise.id = 1) {
  # right now we discard the noise level
  CN <- lapply(1:n.datasets, \(i) {
    # make sure that we are not able to sample the baseline again
    data[c('S0', sample(rownames(data)[1:(noise.id*30)], size = 5)), ]
  }) %>% lapply(\(d) {
    CN <- d
    samples <- rownames(d)
    purity <- sapply(samples, \(s) {
      parameters$purity[parameters$sample == s]
    })
    ratio <- sapply(samples, \(s) {
      parameters$ratio[parameters$sample == s]
    })
    list(CN = CN, purity = purity, ratio = ratio, sample = samples)
  })
}

# ----- Other functions --------------------------------------------------------

# Shuffles the order of the rows in the dataset
shuffle.dataset <- \(dataset) {
  CN <- dataset$CN 
  segments <- dataset$segment 
  N <- nrow(CN) # each column is a sample
  order <- sample(1:N, size = N)
  list(CN = CN[order, ], segment = segments[order, ],
       purity = dataset$purity[order],
       noise.mean = dataset$noise.mean[order],
       noise.sd = dataset$noise.sd[order]
  )
}

