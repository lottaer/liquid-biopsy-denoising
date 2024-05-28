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

# ------------------------------------------------------------------------------
# computes the MSE for the entire dataset
compute.mse.matrix <- \(prediction, true) {
  (abs(prediction - true) %>% rowSums) / ncol(prediction)
}

# corrects the denoised data by the estimated purities
correct.eval <- \(denoised.mat, purities = NULL, post.process = F) {
  if (is.null(purities)) {
    purities <- 
      pbapply(denoised.mat, 1, \(s) compute.purity(s)[2, ])
  }
  if (post.process) {
    sapply(1:nrow(denoised.mat), 
           \(i) correct.profile(denoised.mat[i, ], purities[i]) %>%
             segment.df %>% postprocess.segment %>% get.profile) %>% t
  }
  else {
    sapply(1:nrow(denoised.mat), 
           \(i) correct.profile(denoised.mat[i, ], purities[i])) %>% t
  }
}