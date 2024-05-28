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

#  ----- Load packages ---------------------------------------------------------
library(QDNAseq)

# ----- Load the data -----------------------------------------------------------

# loads all QDNAseq data-files in the folder path.to.files into .Global.Env
load.data <- \(path.to.files = "data/cell-line-50k") {
  .data.files <- list.files(path = path.to.files, pattern = "\\.RData$",
                            full.names = TRUE)
  .N.loaded <- 0
  for (.file in .data.files) {
    .file.name <- tools::file_path_sans_ext(basename(.file))  
    if (!exists(.file.name)) {
      .N.loaded <- .N.loaded + 1
      assign(.file.name, load(.file, .GlobalEnv) %>% mget(., parent.frame()),
             envir = .GlobalEnv)
    }
    else {
      print(paste0(.file, " already loaded. Continuing..."))
    }
  }
  all.variables <- ls(envir = .GlobalEnv)
  if ("copyNumbersSegmented" %in% all.variables) {
    rm("copyNumbersSegmented", envir = .GlobalEnv)
  }
  if ("readCountsFiltered" %in% all.variables) {
    rm("readCountsFiltered", envir = .GlobalEnv)
  }
}

# ----- Functions for handling QDNAseq-objects ---------------------------------

# Gets the raw count matrix from a QDNAseq object
get.counts <- \(read.counts) {
  read.counts$readCountsFiltered@assayData$counts 
}

# Gets the CN-matrix from a QDNAseq object
get.CNVs <- \(copy.numbers) {
  copy.numbers$copyNumbersSegmented@assayData %>% 
    {list(CN = .$copynumber, segment = .$segmented)}
}

# Get positions of chromosomes, using the genome-file in QDNAseq object
get.chromosome.position <- \(data, read.counts = T) {
  if (read.counts) {
    data <- data$readCountsFiltered
  }
  else {
    data <- data$copyNumbersSegmented
  }
  data@featureData@data$chromosome %>% 
    {data.frame(
      chromosome = unique(.),
      start.position = sapply(unique(.), \(chrom) which(. == chrom)[1])
    )}
}

# Get information about gc-content and mappability for normalization from a
# QDNAseq object
get.segment.info <- \(data, read.counts = T) {
  if (read.counts) {
    data$readCountsFiltered@featureData@data 
  }
  else {
    data$copyNumbersSegmented@featureData@data
  }
}

# Removes NA:s from a dataset
remove.missing <- \(df) {
  if(is.null(dim(df))){
    df %<>% as.matrix %>% t
  } else if(dim(df)[1] == 1) {
    df %<>% as.matrix
  }
  to.remove <- apply(df, 2, \(col) any(is.na(col)))
  df[, !to.remove]
}

# Blacklists problematic features in a genome data frame
# genome.df: genome file from QDNAseq object
# blacklist: indices to blacklist
perform.blacklisting <- \(genome.df, blacklist, by.col = TRUE) {
  if(by.col){
    genome.df[, blacklist] <- NA
  } else {
    genome.df[blacklist, ] <- NA
  }
  return(genome.df)
}

# ----- Rename 500k bin files --------------------------------------------------

# Rename to filenames to a readable format
rename.files <- \(path.to.files, copy.numbers = T) {
  file.names <- list.files(path.to.files) %>% get.file.names(copy.numbers)
  for (file in file.names) {
    split <- strsplit(file, "[.]")[[1]]
    # if we have a . in the name, we will rename the file
    if (length(split) > 2) {
      print(paste0("Re-naming: ", file))
      new.name <- paste0(path.to.files, "/", split[1], "_", split[2], 
                         collapse = "", ".RData")
      file.rename(paste0(path.to.files, "/", file), new.name)
    }
  }
}

# ----- Functions reading objects from folder ----------------------------------

# Get the name of either all copy-number files, or readcount files
get.file.names <- \(files, copy.numbers = T) {
  file.start <- ifelse(copy.numbers, "Copynumber", "Readcounts")
  to.keep <- sapply(files, \(file) strsplit(file, "_")[[1]][1] == file.start)
  files[to.keep]
}

# Get filenames of samples corresponding to a certain patient
get.patient.filenames <- \(path.to.files, patient.id = "4968") {
  all.files <- list.files(path = path.to.files)
  patient <- sapply(all.files, \(file) strsplit(file, "_")[[1]][4])
  all.files[patient == patient.id]
}

# Read all files in folder into one dataframe
# copy.numbers: TRUE if CN should be imported, otherwise readcounts
# segments: TRUE if segmented data (by QDNAseq) should be read, CN otherwise
# remove.blacklist: TRUE if we should remove problematic regions
# blacklist.path: path to file containing indices to remove
# all.file.names: list of files to read
files.to.dataset <- 
  \(path.to.files, copy.numbers = T, segments = F, remove.blacklist = TRUE,
    blacklist.path = "data/RData/blacklist_500k.RData",
    all.file.names = NULL) {
  if (is.null(all.file.names)) {
    all.file.names <- list.files(path.to.files)
  }
  file.names <- all.file.names %>% get.file.names(copy.numbers)
  object.names <- sapply(file.names, \(file) strsplit(file, "[.]")[[1]][1])
  loaded.objects <- ls(envir = .GlobalEnv)
  is.loaded <- sapply(object.names, \(obj) obj %in% loaded.objects)
  if (!all(is.loaded)) {
    print("Files not in .GlobalEnv. Loading files...")
    load.data(path.to.files)
  }
  if (copy.numbers) {
     dataset <- object.names %>%
      lapply(\(file) {
        data <- get.CNVs(get(file)) #$CN %>% t
        if (segments) {
          data$segment %>% t
        }
        else  {
          data$CN %>% t
        }
        }) %>% 
      do.call("rbind", .) %>% data.frame   
  } else {
    dataset <- object.names %>%
      lapply(\(file) get.counts(get(file)) %>% t) %>% 
      do.call("rbind", .) %>% data.frame   
  }
  # remove the files that we did not have prior to running the function
  if (!all(is.loaded)) {
    print("Cleaning loaded files...")
  }
  to.remove <- setdiff(ls(envir = .GlobalEnv), loaded.objects)
  rm(list = to.remove, envir = .GlobalEnv)
  if (path.to.files == "data/cell-line-50k" | 
      path.to.files == "data/cell-line-500k") {
    to.remove <- which(rownames(dataset) == "E5_2")
    dataset <- dataset[-to.remove, ]
  }
  if (remove.blacklist) {
    if (!file.exists(blacklist.path)) {
      print("Could not find extended blacklist file. Remove bins manually")
    }
    else {
      to.remove <- readRDS(blacklist.path)
      dataset <- dataset[, -to.remove]
    }
    dataset %<>% remove.missing
  }
  dataset
}

# Get sample from rownames
get.sample <- \(cn.samples, sample.name = "A1_1") {
 cn.samples[which(rownames(cn.samples) %in% sample.name), ] %>% 
    as.matrix
}

# Read the list of indices to blacklist
get.extended.blacklist <- \(path.to.blacklist = NULL) {
  if(is.null(path.to.blacklist)) {
    path.to.blacklist <- "data/RData/blacklist_500k.RData"
  }
  readRDS(path.to.blacklist)
}

# ----- Outlier detection ------------------------------------------------------

# finds outliers in the data, either for the full region or for only the
# lower or upper half
find.outliers <- \(data, n.sd = 3, region = "full") {
  if(region == "full"){
    which(abs(data - median(data)) > (n.sd * sd(data))) %>% return
  } else if(region == "lower") {
    which((data - median(data)) < - (n.sd * sd(data))) %>% return
  } else if( region == "upper") {
    which((data - median(data)) > (n.sd * sd(data))) %>% return
  } else {
    stop("invalid region")
  }
}

# help function for find.outliers.df
.outlier.df <- \(outlier.data, name) {
  data.frame(dataset = rep(name, length(outlier.data)), x  = outlier.data)
}

# identifies outliers for a data set
find.outliers.df <- \(imputed.df, n.sd = 3, region = "full", names = NULL) {
  outliers <- apply(imputed.df, 1, find.outliers, n.sd = n.sd, region = region)
  if(is.null(names)) {
    names <- names(outliers)
  }
  outliers.df <- lapply(1:length(outliers), \(i) 
                        .outlier.df(outliers[[i]], names[i])
  ) %>% do.call(rbind, .)
  return(outliers.df)
}

# plot outlier data for different data sets
plot.outliers <- \(df, range = NULL) {
  if(!is.null(range)) {
    df <- df[which(df$x %in% range),]
  }
  df %>% ggplot(aes(x = x, y = dataset)) + 
    geom_point(alpha = 0.2, color = "#56B4E9") +
    theme(axis.title = element_blank())
}
