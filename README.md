# Deconvolution methods for quantifying copy number variations in liquid biopsy sequencing

This repository contains the code used in the thesis _Deconvolution methods for quantifying copy number variations in liquid biopsy sequencing_, that evaluates methods for devoncolving copy number variations (CNV) profiles from noisy liquid biopsy sequencing data. We explore two methods for this purpose: a denoising autoencoder and Bayesian change point detection. The thesis workflow is presented in `main.R`.

We denoise the signals, and then segment them, i.e. divide them into regions of homologous copy number state.


The methods are developed for liquid biopsy sequencing data, where the sequenced DNA is a mixture of both healthy DNA (without genomic alterations) and tumor DNA (containing genomic alterations). We implement a variation of the liquidCNA-algorithm for estimating sample purity and removing normal contamination, i.e. contribution of healthy cells to the overall copy number profile. 

For more details, see the manuscript [TBA]

## Contents

The files are organized as follows: 

- `autotuning.R`: methods for tuning model parameters in denoising autoencoder and rolling median tuning.
- `denoising.R`: methods for denoising: including rolling median denoising, Bayesian change point detection, constructing denoising autoencoder and making predictions.
- `estimation_functions.R`: functions for estimating sample purity and subclonal ratio.       The estimation of the sample purity and subclonal ratio is based on `liquidCNA` created by [Lakatos et al](https://github.com/elakatos/liquidCNA).
- `evaluation.R`: helper functions for method evaluation.
- `handle_data.R`: methods for reading data from disk.
- `main.R`: contains the thesis workflow.
- `models.py`: Python implementation of denoising autoencoder.
- `segmentation.R`: method for segmenting denoised signals.
- `simulate_data.R`: simulates CNV data inspired by high grade ovarian cancer data, at different purities and noise levels. It is also possible to simulate longitudinal datasets used for subclonal ratio estimation.

## Example data

The folder `example_data` contains in silico mixed samples of blood and ovarian cancer, see the [Lakatos et. al. (2021)](https://www.sciencedirect.com/science/article/pii/S2589004221008579) for further reference.

## License
Estimation of purity and subclonal ratio from liquid biosy sequencing. Inspired by Lakatos et. al., <https://github.com/elakatos/liquidCNA>.

Copyright (C) 2024  Lotta Eriksson lottaer@chalmers.se & Linnea Hallin hallinl@chalmers.se

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.