README for UTM/UFIT aux subdirectory.

This directory contains additional programs with varying degress of polishedness that may be useful for running UTM/UFIT, or for the interpretation of UFIT results. These programs may need modification to adapt to specific needs; see program headers.

- ufitrandomrep.pro :  runs UFIT repeatedly with random offsets of initial values. Typically used for AMOEBA fits. 

- ufitloop.pro : example of a wrapper that steps UFIT through a set of initial values; (in this case, the parameter 1radi).  Needs modification for any use.

- utmparplot.pro :  reads a list of .utm setupfiles and displays distributions and correlations of selected parameters.  Intended to interpret sequences of  output (.fout ) setup-files from ufitrandomrep.pro  or ufitloop.pro. 

- rdutmpar.pro : reads a list of .utm setupfiles and saves selected paramters to text file (.csv format), for reading into Excel or similar. Intended for output (.fout setup-files) from ufitrandomrep.pro 


- ufit_mcshow.pro : displays .sav output file from MCMC runs of UFIT. Will show temporal evolution of a selected parameter and generates plots of distributions and correlations of the free parameters. This program calls the following ones in sequence:

    - parevol.pro : Shows temporal evoltion of parameters from MCMC plot.

    - plot_par1_vs_params: scatterplots of values of one parameter (on Y axis, intended to be chisq or stddev) versus all other parameters values. 

    - ufit_exofast_plotdist.pro : Plots distributions of parameters; is modified version of exofast_plotdist.pro. Can be used standalone; see main-level example program at end. Is called by ufit_mcshow.pro.

    - scattermatrix.pro :  Creates a matrix of scatter plots or correlations. By default data are shown as dots, but contour lines indicating confidences are also possible. Can be used standalone; see main-level example program at end. Is called by ufit_mcshow.pro. 

    - ufit_exofast_errell : Routine needed by scattermatrix.pro (gives the path of the error ellipses at constant probability given two parameter distributions.)

- ufit_limitshow: reads a setupfile and the correspoding .mc.sav output from an MCMC run and plots an overlay with all models with a chi-square within 1 (=1sigma) of the best fit. 


- ufit_mccomb.pro : Combines mc.sav output files from several MCMC runs of UFIT. 

- mcextract.pro  ;template code to extract some specific param-values from mc.sav file

- rdtab.pro : A flexible reader for alphanumeric tables extracting up to 30 columns

