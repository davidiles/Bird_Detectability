The scripts in this .zip file contain the necessary R code for replicating the analysis of golden-crowned sparrow (GCSP) analysis found in:

Martin-Schwarze, A., Niemi, J. & Dixon, P. Joint Modeling of Distances and Times in Point-Count Surveys. JABES (2021). https://doi.org/10.1007/s13253-021-00437-3  

The analysis is a model of  animal abundance from point-count surveys utilizing both distance sampling and removal sampling.  Importantly, the model provides a flexible density function for the joint distribution of times and distances to first detection.  The proposed model features Bayesian MCMC sampling based on a hierarchical N-mixture model featuring:

* a Poisson-normal model for abundance
* a constant-rate model for availability
* a half-normal distance function for IJ (animal-based) perceptibility, and
* a half-normal distance function for PP (cue-based) perceptibility.

The code requires (and the .zip file includes) data downloaded from the USGS Alaska Science Center:

Amundson, C. L., Handel, C. M., Ruthrauff, D. R., Tibbitts, T. L. and Gill, R. E., Jr , 2018, Data for Montane-breeding Bird Distribution and Abundance across National Parks of Southwestern Alaska, 2004-2008: U.S. Geological Survey data release, https://doi.org/10.5066/F7MW2GD3

The two data files are:

"montaneBird_distribution_SW_obsdata_amundson_2004_2008.csv"
and
"montaneBird_distribution_SW_sitecovs_amundson_2004_2008.csv"

The proper order to run code is presented in the umbrella script 'Code_Roadmap.R', which should successfully reproduce the entire analysis.
