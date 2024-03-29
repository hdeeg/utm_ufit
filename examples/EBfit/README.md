This folder contains an example of a binary model resp. a fit to the eccentric EB HD38600. Its lightcurve (hd38600ph.lc) is a phase-folded lightcurve based on data obtained by the TESS mission. The lightcurve was trimmed to sections around the two eclipses and normalized to 1 in the off-eclipse secton. It was also rebinned to phase increments that correspond to intervals of ~10mins.

To plot a model lightcurve (which corresponds very closely to the observed lightcurve) use:\
IDL> utm,'hd38600fit.utm',3  \
Set the parameter 'wflag 1'  to save a text-version of the model-curve

To generate some statistics about the agreement between the input model and the lightcurve, use:

IDL> ufit,'hd38600fit.utm',/nofit

To run an MCMC with 10 chains and 5000 steps, use:\
IDL> ufit,'hd38600fit.utm',2

(for good parameter distributions, set 'mc_maxsteps  50000'. The Gelman-Rubin mixing criterium should be fulfilled at 12000-15000 steps, taking 1-2hours)

The output chain of parameters can be visualized with:\
IDL> ufit_mcshow,'hd38600fit.mc0.sav'   (see /aux folder for the program) 

To run amoeba fits, change parameter 'fitalgo' to 'am' and modify any of the starting parameters (the current excellent fit will not improve from further amoeba runs). The scale parameters in hd38600fit.utm have been set for use with amoeba for poorly known input parameters. Also, amoeba will likely not find a model of a quality similar to the original model; this will rather succeed from one or several MCMC runs. 

--
Explanations on the setup:

Setup file hd38600fit.utm  is an edited version of a setup that was  autogenerated from previous fits, and is hence uncommented. This setupfile makes use of the preprocessor bin_prep_k2.pro, which relates several parameters between the components (namely period, eccentricy, omega, inclination, eclipse-epoch of the primary component with transit-epoch of the secondary) and derives the component's further parameters from a parameterization consisting of:

- krad :        secondary /primary radius ratio,  r2/r1
- (r1+r2)/a    : sum of component radii over semimajor axis between components
- qlbin :        secondary /primary luminosity ratio L2/L1
- contlum:        fraction of light from a third star, with  contlum = L3/(L1+L2+L3)

The preprocessor converts these parameters into the component's individual parameters, using normalizations of: 
- L1+L2+L3=1    (corresponding to the observed lightcurve's flux normalization)
and 
- a=1=a1+a2  
The ratio a2/a1 doesn't matter in a binary fit and a1=a2=0.5 is used. For models of triple systems however, the barycenter needs to be defined. Outcommented code is therefore present that takes a mass ratio 'qmass' as input parameter, from which a1 and a2 will be calculated.

Note that the period is 1, since the curve is in units of orbital phase. Above normalization can however also be used with timeseries in units of day, requiring then the real period for '0period', and 'tunit' should be set to 'day'.




