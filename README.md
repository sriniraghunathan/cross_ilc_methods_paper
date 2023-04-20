# Cross-ILC methods paper

## Paper:
* Raghunathan & Omori 2023; arXiv:[2304.09166](https://arxiv.org/abs/2304.09166); Submitted to ApJ.

## Overview:
* Cross-ILC technique:
  * For kSZ: Compute cross-spectrum between tSZ-free and CIB-free ILC maps to measure the kSZ power spectrum.
  * For CMB lensing: Pass tSZ-free and CIB-free ILC maps in the two legs of lensing quadratic estimator to mitigate both tSZ and CIB-induced lensing biases.
  
## Notebooks:
* [make_paper_plots.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/make_paper_plots.ipynb): Reproduce plots in the paper.
* Fisher formalism:
  * [s1_get_fisher_matrix.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/s1_get_fisher_matrix.ipynb): Forecasts total SNR of the kSZ power spectrum and also cosmological / reionisation constraints for current and future CMB surveys.
  * [s2_analyse_fisher_matrix.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/s2_analyse_fisher_matrix.ipynb): Analyse the results from [s1_get_fisher_matrix.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/s1_get_fisher_matrix.ipynb). This notebook will help in making tables and figures.
* [get_bandpower_errors.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/get_bandpower_errors.ipynb): Obtain the bandpower errors for different ILC techniques / CMB experiments analytically using Knox formula.
* [estimate_bias_using_fisher.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/estimate_bias_using_fisher.ipynb): Estimate biases in kSZ and cosmological parameters due to unmodelled CIB/tSZ residuals. This used a modifed Fisher formalism based on [Amara & Refregier 2007](https://arxiv.org/abs/0710.5171).
   * There are also couple of more example notebooks to understand Fisher-based bias estimation in [fisher_bias](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/b6eef6e608ef574324ce874112d2db7f638efc29/fisher_bias) folder.
     * [example1_bias_in_line_fitting.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/b6eef6e608ef574324ce874112d2db7f638efc29/fisher_bias/example1_bias_in_line_fitting.ipynb) - Bias estimation in simple line fitting and a comparison of the result using a lileihood approach.
     * [example2_bias_in_cosmo.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/b6eef6e608ef574324ce874112d2db7f638efc29/fisher_bias/example2_bias_in_cosmo.ipynb) - Inject biases as xx $\sigma$ of cosmological parameters and recover them using the modified Fisher appraoch. 
* [get_radio_residuals.ipynb](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/get_radio_residuals.ipynb): Modelling the power spectrum of radio sources and also get the ILCed radio residuals analytically.

## Data products:
* [Data](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/data):
  * [lagache_2019_ns150_radio.dat](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/lagache_2019_ns150_radio.dat): Radio \$dN/dS$ = Lagache et al. 2020 model [(arXiv: 1911.09466)](https://arxiv.org/abs/1911.09466) used to model radio point source power.
  * [ilc_weights_residuals_agora_fg_model.npy](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/ilc/ilc_weights_residuals_agora_fg_model.npy): ILC residuals for multiple ILC combinations and all experiments considered in this work.
  * [CAMB](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/data/CAMB): CMB power spectra and derivatives of the power spectra as a function of cosmological parameters:
      * Both [*Planck* 2018](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/planck_2018/) and [*Planck* 2015](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/planck_2015/) are available.
      * [*Planck* 2018](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/planck_2018/) is the default.
    * These are used for [Fisher](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/s1_get_fisher_matrix.ipynb) forecasting.
    * Both [lensed](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/planck_2018/cmb_spectra_lensed.txt) and [unlensed](https://github.com/sriniraghunathan/cross_ilc_methods_paper/blob/main/publish/data/CAMB/planck_2018/cmb_spectra_unlensed.txt) versions are available.
    * Fiducial: Lensed spectra.
* [ILC results](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/ilc):

* [Fisher results](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/fisher):

* [Lensing results](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/lensing):
   * [Lensing N0](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/lensing/lensing_n0): Results available for standard QE, [MH18](https://arxiv.org/abs/1802.08230), and the cross-ILC.
     * Simply do: `ls *qewithcmbmv*` or `ls *mh*` or `ls *xilc*` to see the files in the linked [Lensing N0](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/lensing/lensing_n0) folder.
   * [Lensing cross-correlations](https://github.com/sriniraghunathan/cross_ilc_methods_paper/tree/main/publish/lensing/lensing_xls):
 

## Dependencies:
* Standard python packages: jupyter notebook, numpy, scipy, matplotlib.
* For Fisher bias calculations:
  * CAMB.
