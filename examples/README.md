This directory contains example setups and (in some cases) the required data and preprocessor-codes to perform a variety of UTM or UFIT simulations or fits.

*.utm files are example setup files. The names of most are self-explanatory; see the description in their headers. Just try them out with UTM. Some are mentioned in the documentation.

*_prep.pro  are preprocessor files  (for use, see kepler_earth_sun.utm  or K446spot.utm)

mktestlc.utm and fitest.utm: Basic fit with UFIT, see the fitting example in the UFIT documentation.

EBfit : example of a model fit to an eccentric binary from TESS data

prior_examples : demos for the use of priors in UFIT

HD144548triple : setup (triple_HD144548.utm ) to generate the model used in the Alonso et al. 2015 paper on the eclipsing triple system HD144548.  HD144548eclips.lc resp. HD144548full.lc are the corresponding processed K2 lightcurves.

extern_body : examples to generate transits from objects with arbitrary shapes

alt_modellers : example setups for using UFIT with alternative modellers to UTM (polynomes, sinusoidals)

