This directory contains all code needed to run UTM/UFIT and associated programs.
The IDL_PATH of your IDL installation should include this directory and its subdirectories.

Its content is:
(The UTM/UFIT core-files are in the root of this directory)
utm.pro	      	source for UTM (*)
ufit.pro      	source for UFIT
utm_ancil.pro    collection of anciliary routines used by UFIT, UTM
setupfile.pro       required for dealing with setup files	

/exofast     original and modified EXOFAST routines that are required for UTM/UFIT

/aux         Several additional programs, some on prototype level, that are useful for running UTM/UFIT or for the interpretation of UFIT results. See _Readme.txt in there.

/modellers    Modeller routines that are callable via setup-file by UFIT and which are alternatives to UTM; e.g. for polynome fitting. They are also useful as templates to incorporate further fitting functions. See also /examples/alt_modellers. 



