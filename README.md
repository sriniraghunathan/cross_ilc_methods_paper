# Cross-ILC methods paper

## Paper:
* ....

## Overview:
* Cross-ILC technique:
  * For kSZ: Compute cross-spectrum between tSZ-free and CIB-free ILC maps to measure the kSZ power spectrum.
  * For CMB lensing: Pass tSZ-free and CIB-free ILC maps in the two legs of lensing quadratic estimator to mitigate both tSZ and CIB-induced lensing biases.
  
## Notebooks:
* [make_paper_plots.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/make_paper_plots.ipynb): Reproduce plots in the paper.
* [perform_ksz_fisher.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/perform_ksz_fisher.ipynb): Forecasts total SNR of the kSZ power spectrum for current and future CMB surveys.
* [get_radio_residuals.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/get_radio_residuals.ipynb): Model power spectrum of radio sources and also get the ILCed radio residuals analytically.

## [Data products](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/data):
* [lagache_2019_ns150_radio.dat](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/lagache_2019_ns150_radio.dat): Radio \$dN/dS$ = Lagache et al. 2020 model [(arXiv: 1911.09466)](https://arxiv.org/abs/1911.09466) used to model radio point source power.
* [ilc_weights_residuals_agora_fg_model.npy](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/ilc_weights_residuals_agora_fg_model.npy): ILC residuals for multiple ILC combinations and all experiments considered in this work.
* [cmb_cl_planck_2015_lensedCls.dat](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/cmb_cl_planck_2015_lensedCls.dat): {\it Planck} 2015 CMB power spectra.
* [CAMB](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/data/CAMB): CMB power spectra and derivatives of the power spectra as a function of cosmological parameters:
  * These are used for [Fisher](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/perform_ksz_fisher.ipynb) forecasting.
  * Both [lensed](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/cmb_spectra_lensed.txt) and [unlensed](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/cmb_spectra_unlensed.txt) versions are available.
  * Fiducial: Lensed spectra.

## Dependencies:
* Standard python packages: jupyter notebook, numpy, scipy, matplotlib.
