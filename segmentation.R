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
library(cumSeg)

# ----- Post-processing of the segmentation (optional) -------------------------

# the value assigned to the segment is the weighted avg. of the merged segments
.weighted.average <- \(segment.df, first, second) {
  values <- segment.df$value
  lengths <- segment.df$end - segment.df$start + 1
  total.length <- lengths[first] + lengths[second]
  (lengths[first] * values[first] + lengths[second] * values[second]) / 
    total.length
}

# in the liquidCNA-code we remove bins that are shorter than 12*500kbp
postprocess.segment <- \(segment.df, min.bins = 12) {
  width <- segment.df$end - segment.df$start + 1
  to.merge <- which(width < min.bins)
  while (length(to.merge) > 0) {
    segment <- to.merge[1]
    curr <- segment.df$value[segment]
    # if last segment, merge into the last
    if (segment == nrow(segment.df)) {
      segment.df$end[segment - 1] <- segment.df$end[segment]
      segment.df$value[segment - 1] <- 
          .weighted.average(segment.df, segment - 1, segment)
    }
    # if first segment, merge into the next segment
    else if (segment == 1) {
      segment.df$value[segment + 1] <- 
        .weighted.average(segment.df, segment, segment + 1)
      segment.df$start[segment + 1] <- segment.df$start[segment]
    } else {
      diff.before <- abs(segment.df$value[segment - 1] - curr)
      diff.after <- abs(segment.df$value[segment + 1] - curr)
      # the new mean value of the segment will be weighted avg. of the segments
      if (diff.before < diff.after) {
        segment.df$value[segment - 1] <-
          .weighted.average(segment.df, segment - 1, segment)
        segment.df$end[segment - 1] <- segment.df$end[segment]
      } else {
          segment.df$value[segment + 1] <- 
            .weighted.average(segment.df, segment, segment + 1)
          segment.df$start[segment + 1] <- segment.df$start[segment]
      }
    }
    # remove merged segment
    segment.df <- segment.df[-c(segment), ]
    # recompute segment lengths and segments to merge
    width <- segment.df$end - segment.df$start + 1
    to.merge <- which(width < min.bins)
  }
  segment.df
}

# remove unwanted segments instead of postprocessing them
filter.segments <- \(data, min.bins = 12) {
  # all the segments that we would merge, we can just throw away
  width <- data$end - data$start + 1
  to.remove <- which(width < min.bins)
  data$value[to.remove] <- NA # to be able to just remove NA
  data 
}

# compute the union of the breakpoints from longitudinal data and remove the 
# ones that are shorter than min.bins
get.long.bp <- \(segment.df.list, min.bins = 12) {
  starts <- sapply(segment.df.list, \(x) x$start) %>% unlist %>% unique %>% sort
  ends <- sapply(segment.df.list, \(x) x$end) %>% unlist %>% unique %>% sort
  breakpoints <- data.frame(start = starts, end = ends)
  breakpoints %<>% .[((.$end - .$start) >= min.bins),]
  return(breakpoints)
}

# ----- Extract the signal from breakpoints ------------------------------------

# construct the segment from the segment.df
get.profile <- \(seg.df) {
  profile <- numeric(seg.df$end[nrow(seg.df)])
  for (i in 1:nrow(seg.df)) {
    profile[(seg.df$start[i]):(seg.df$end[i])] <- seg.df$value[i]
  }
  profile
}

# finds the breakpoints from CNA data, returned as a boolean array
get.breakpoints <- \(cna.profile) {
  c(F, sapply(2:length(cna.profile), \(i) cna.profile[i] != cna.profile[i-1])) 
}

# construct a segment dataframe from the entire copy number profile
# contains the start and end of breakpoints and the segment values
segment.df <- \(cn.profile) {
  is.breakpoint <- get.breakpoints(cn.profile)
  starts <- which(is.breakpoint)
  data.frame(
    start = c(1,starts),
    end = c(starts-1, length(cn.profile)),
    value = cn.profile[c(1,is.breakpoint %>% which)]
  )
}

# construct a segment data frame from set breakpoints
seg.from.bps <- \(cn.sample, breakpoints) {
  values <- sapply(1:nrow(breakpoints), \(i) {
    cn.sample[breakpoints$start[i]:breakpoints$end[i]] %>% median
  })
  data.frame(
    start = breakpoints$start,
    end = breakpoints$end,
    value = values
  )
}

# ----- Do model selection for selecting the number of breakpoints -------------

# plot the criterion, e.g. BIC curve to validate the selection rule
.plot.criterion <- \(values, y.lab = "BIC", selected = NULL) {
  gg <- data.frame(
    n.segments = 0:(length(values) - 1),
    criterion = values
  ) %>% ggplot(aes(x = n.segments, y = criterion)) + 
    geom_line(col = "#999999") +
    geom_point(col = "#999999") + labs(x = "number of breakpoints", y = y.lab)
  if (!is.null(selected)) {
    gg <- gg + geom_vline(xintercept = selected, 
                          col = "red", linetype = "dotted")
  }
  gg
}

# Defines the selection criteria
.compute.BIC <- \(object, n, Cn = \(n) log(log(n))) {
  RSS <- object$RSS
  edf <- 1 + 2 * (object$df - 1)
  log(RSS / n) + edf * log(n) / n * Cn(n)
}

# Function for selecting the number of breakpoints
.select.n.segments <- \(values, tol = 0.01) {
  values <- (values - min(values)) / (max(values - min(values)))
  d <- diff(values, lag = 2)
  above <- which(d < -tol)
  above[length(above)] + 1
}

# NOTE : the code below is adapted from the authors of cumSeg-package
# see their Github for the entire code: https://github.com/cran/cumSeg
# the function is used to re-compute the segment means after selecting the 
# number of breakspoints
seg.lm.fit0<- \(y, Z, PSI, control, round = FALSE, ...){
  it.max <- old.it.max <- control$it.max
  toll <- control$toll
  visual <- control$visual
  last <- control$last
  stop.if.error<-control$stop.if.error
  h <- min(abs(control$h), 1)
  if (h < 1)
    it.max <- it.max + round(it.max/2)
  it <- 1
  epsilon<-10
  k <- ncol(PSI)
  psi.values <- NULL
  H <- 1
  psi<-PSI[1,]
  XREG<-cbind(Z[,1])
  obj<-list(residuals=rep(10,3))
  while (abs(epsilon) > toll) {
    U <- pmax((Z - PSI), 0)
    V <- ifelse((Z > PSI), -1, 0)
    X<-cbind(XREG, U, V)
    dev.old <- sum(obj$residuals^2)
    rownames(X) <- NULL
    if (ncol(V) == 1) {
      colnames(X)[ (ncol(XREG) ):ncol(X)] <- c("firstSlope","U", "V")
    } else {
      colnames(X)<-rev(c(paste("V", ncol(V):1, sep = ""),
                         paste("U",ncol(U):1, sep = ""),
                         rep("firstSlope",ncol(X)-ncol(U)-ncol(V))))
    }
    obj <- lm.fit(x = X, y = y) #drop(solve(crossprod(X),crossprod(X,y)))
    dev.new <- sum(obj$residuals^2)
    if (visual) {
      if (it == 1) cat(0, " ", formatC(dev.old, 3, format = "f"),"",
                       "(No breakpoint(s))", "\n")
      spp <- if (it < 10) "" else NULL
      cat(it, spp, "", formatC(dev.new, 3, format = "f"), "---",ncol(V),
          "breakpoints","\n")
    }
    epsilon <- (dev.new - dev.old)/dev.old
    obj$epsilon <- epsilon
    it <- it + 1
    obj$it <- it
    class(obj) <- c("segmented", class(obj))
    if (k == 1) {
      beta.c <- coef(obj)["U"]
      gamma.c <- coef(obj)["V"]
    } else {
      #se ci sono contrasti i beta.c quali sono?
      beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
      gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]
    }
    if (it > it.max) break
    psi.values[[length(psi.values) + 1]] <- psi.old <- psi
    if (it >= old.it.max && h < 1) H <- h
    psi <- round(psi.old + H * gamma.c/beta.c,0)
    PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
    #check if psi is admissible..
    a <- apply((Z <= PSI), 2, all)
    b <- apply((Z >= PSI), 2, all)
    if(stop.if.error) {
      if(sum(a + b) != 0 || is.na(sum(a + b))) stop("(Some) estimated psi out of its range")
    } else {
      id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
      Z <- Z[,id.psi.ok,drop=FALSE]
      psi <- psi[id.psi.ok]
      PSI <- PSI[,id.psi.ok,drop=FALSE]
    }
    if(ncol(PSI)<=0) {
      warning("No breakpoint estimated", call. = FALSE)
      obj<-lm.fit(x = XREG, y = y)
      obj$fitted.values<-rep(obj$coef,length(y))
      obj$est.means<-obj$coef
      return(obj)
    }
  } # end while
  if(round) {
    psi<-round(psi,0)
    PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
    V <- ifelse((Z > PSI), -1, 0)
  }
  obj$psi <- sort(psi)
  obj$beta.c <- beta.c[order(psi)]
  obj$gamma.c <- gamma.c[order(psi)]
  obj$epsilon <- epsilon
  obj$V<- V[,order(psi)]
  obj$psi<-obj$psi[!is.na(obj$beta.c)]
  obj$V<-as.matrix(as.matrix(obj$V)[,!is.na(obj$beta.c)])
  obj$beta.c<-obj$beta.c[!is.na(obj$beta.c)]
  obj
}


# NOTE : the code below is adapted from the authors of cumSeg-package
# see their Github for the entire code: https://github.com/cran/cumSeg
compute.segmentation <- \(denoised, K = 200, show.plot = FALSE, 
                          criterion = .compute.BIC) {
  # overestimate the number of change points
  init.splits <- jumpoints(denoised, output = "1", k = K)
  # the length of psi in init.splits is typically lower than K as some values
  # are disgarded in the iterative process.
  max.breaks <- init.splits$n.psi
  stepwise.fit <- lars(abs(init.splits$V), denoised, type = "stepwise",
                       normalize = FALSE, intercept = TRUE)
  crit <- criterion(stepwise.fit, n = length(denoised))
  # return the model selected
  n.splits <- .select.n.segments(crit)
  if (length(n.splits) == 0) {
    n.splits <- 1
  }
  if (show.plot) {
    plot(.plot.criterion(crit, selected = n.splits - 1))
    pl <- plot(.plot.criterion(crit, selected = n.splits - 1))
  }
  id.var.entry <- (1:ncol(init.splits$V))[order(stepwise.fit$entry)]
  id <- sort(c(0, id.var.entry)[1:n.splits])
  id %<>% {.[-1]} # exclude intercept
  psi <- init.splits$psi[id]
  V.keep <- cbind(1, abs(init.splits$V[, id]))
  k <- length(psi)
  n <- length(denoised) 
  x <- 1:n
  Z <- matrix(rep(x, k), nrow = n)
  PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
  obj <- seg.lm.fit0(y = cumsum(denoised), Z = Z, PSI = PSI, 
                   control = fit.control(toll = 1e-04, it.max = 20),
                   round = F)
  fitted.v<-drop(abs(obj$V)%*%obj$beta.c)
  if ("firstSlope" %in% names(coef(obj))) {
    fitted.v <- fitted.v + obj$coef["firstSlope"]
  }
  return(fitted.v %>% segment.df)
}

# ----- Plot denoised sample and segmentation ----------------------------------

# Plots the segmentation together with the sample ot the ground truth
# sample: the sample
# true.segment: the ground truth
# segmentation: already segmented data
# post.process: should the segment be postprocessed? TRUE if yes
# data: original data
plot.segments <- \(sample, true.segment = NULL, segmentation = NULL, 
                   post.process = FALSE, data = NULL, ...) {
  if (is.null(segmentation)) {
    segmentation <- compute.segmentation(sample) 
  }
  if (post.process) {
    segmentation %<>% post.process.old(segmentation, data)
  }
  p <- data.frame(
    index = 1:length(sample),
    value = sample %>% as.numeric
  ) %>% ggplot(aes(x = index, y = value)) + 
    geom_line(color = "lightgray") +
    geom_segment(data = segmentation, 
                 aes(x = start, y = value, xend = end, 
                     yend = value), linewidth = 1.2, color = "red")
  if (!is.null(true.segment)) {
    df.true <- data.frame(index = 1:length(true.segment), 
                          value = true.segment %>% as.numeric)
    p <- p + geom_line(data = df.true, aes(x = index, y = value), 
                       linetype = "dashed", color = "darkgrey")
  }
  p
}
