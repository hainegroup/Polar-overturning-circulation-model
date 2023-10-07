Polar Overturning Circulation model
==============================
[![Build Status](https://travis-ci.com/ThomasHaine/polar_overturning_circulation_model.svg?branch=master)](https://travis-ci.com/ThomasHaine/polar_overturning_circulation_model)
[![codecov](https://codecov.io/gh/ThomasHaine/polar_overturning_circulation_model/branch/master/graph/badge.svg)](https://codecov.io/gh/ThomasHaine/polar_overturning_circulation_model)
[![License:MIT](https://img.shields.io/badge/License-MIT-lightgray.svg?style=flt-square)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/257373994.svg)](https://zenodo.org/badge/latestdoi/257373994)


Polar Overturning Circulation (POC) model MATLAB code and `sympy` theory for H21, plus subsequent development. 

  * Run `run_POC_model.mlx` to run experiments and build figures.
  * Run `run_glacial_POC_model.mlx` to run LGM experiments and build figures (under development).
  * Run `run_H21_POC_model.mlx` to run the suite of experiments and build figures from H21 paper.
  * Run `interactive_POC_model.mlapp` app.
  * Run `symbolic_matrix_inverse.ipynb` notebook for POC theory used in the H21 paper.
  * Run `POB_entrainment_parametrization.mlx` for a comparison with the Price & O'Neil Baringer (1994) entrainment model parameters.
  * Run `make_FSBSO_TS_figure.m` to make the H21 figure of data at Fram Strait and the Barents Sea Opening.
  
 The MATLAB code uses [Gibbs-Seawater (GSW) Oceanographic Toolbox functions](http://www.teos-10.org/software.htm#1) (it needs to be installed). It was written in MATLAB version 2018b and has been tested in MATLAB 2023a.
 
A paper (H21) entitled [*A Conceptual Model of Polar Overturning Circulations*](https://journals.ametsoc.org/view/journals/phoc/51/3/JPO-D-20-0139.1.xml) has been published by JPO.
