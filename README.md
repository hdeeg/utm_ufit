# UTM/UFIT Universal Transit Modeller and Universal Fitter

UTM (Universal Transit Modeller)) is a highly flexible light-curve simulator for all kinds of mutual transits or eclipses between arbitrary numbers or stars, planets, moons, planetary rings, on Keplerian orbits that may be hierarchical and may include osculating elements. Common limb darkening laws are provided, and oblique object shapes or star-spots may be modelled as well. Resampling for the analysis of slowly sampled light curves and the modelling of different types of out-of-transit brightness variations is also possible. All input parameters are provided in setup-files. Short codelets permit the adaptation to any parameter set used in the modelling of transits or eclipses. Transit fluxes are calculated with analytical routines from EXOFAST (by Jason Eastman) or through pixelized object representations; by default, UTM selects automatically the appropriate calculation mode. UTM provides configurable textual and graphical output, with optional animations of transit events. 

UFIT (Universal Fitter) is a fitter for UTM with an interface for alternative modellers; some examples are provided. UFIT can fit any parameter that is specifiable to the modeller, and uniform or gaussian priors may be given. Fitting may be done with AMOEBA or with an MCMC routine based on EXOFAST. UFIT's best-fit parameters are written to setup-files with an identical format as the input to UTM or UFIT. These can be used in further invocations of UFIT, with a versioning that facilitates the documentation of the work performed. Several levels of graphical and printed output can be specified. For the MCMC output, analysis routines are provided which generate corner and distribution plots and several statistical indicators.

UTM/UFIT is written in IDL, and is released under the GNU General Public License. A thorough set of documentation is included, together with numerous applications examples. 

## Package Content
The subdirectories are: 
/code   All required IDL code  
/documentation      Documentation of UTM/UFIT 
/examples       Files requires for a variety of example runs of UTM or UFIT

For more details, see README files in the subdirectories


## Installation:
IDL Version 8.0 or higher is required. The code does not work with GDL.

As the root for the UTM/UFIT install, create a directory that is either:
- A subdirectory of a directory that is already pointed to by the IDLPATH. To see the current path, in the IDL cmd-line type: `print,!path`   
- If you want to install UTM/UFIT in some other directory:
    Create the directory, e.g. `yourpath/utm_ufit'`
    If using the the bash-shell, add the follwing line to the .bashrc file:
    ``export IDL_PATH="<IDL_DEFAULT>:+yourpath/utm_ufit/"``

Either download the ZIP file (see green 'Code' button) or clone this github repositary in your newly created directory. If you download the ZIP, it gets expanded into a directory named `utm_ufit-master`. You should then move the content of the `code` subdirectory into the previously created directory.

## Getting started:
`cd` to the `examples` directory, start the IDL command line, and run UTM with any of the `.utm` setup-files provided, e.g.
`IDL> utm,'ring.utm',3 `         where `3` controls the level of display.

Switching from one to another example, you get usually an error of the kind:
`% Conflicting data structures: <BYTE      Array[385, 385]>,STAR.
% Execution halted at: UTM   1295`
This means that the next example uses different-sized object arrays. This requires to give a `.reset` command to IDL.

Demo of a fit of a transit model:
Generate a noisy lightcurve (named `testtrans.lc` ) of a transit by typing
`IDL> utm,'mktestlc.utm',3`
Then fit a model to it, using an AMOEBA fit:
`IDL> ufit,'fitestlc.utm',3`
The fitted paramters are displayed and also written to an output setup-file named `fitestlc.fout0` (with a duplicate named `fitestlc.f_in1`). The fitted model curve is saved as `fitestlc.fout0.utm`
You may compare the fitted paramters (`1radi, 1trepoch, 1inclin` which are the planet size, epoch and inclination) in `fitestlc.fout0` with the fit's (poor) input parameters in `fitestlc.utm` and with the original parameters in `mktestlc.utm`

To run an MCMC using the parameters from the previous AMOEBA fit as a starting point:
On the duplicate, `fitestlc.f_in1.utm`, use a text-editor to change the line with `fitalgo am` to  `fitalgo mc`.  
Then run ufit using this setup:
`ufit,'fitestlc.f_in1.utm, 2`    (a display level of 1 or 2 is recommended for speed) Execution will need several minutes

Best fit parameters are shown again and saved to file `fitestlc.fout1.lc` (with a further duplicate named `fitestlc.f_in2.utm`). The MC chain is saved as `fitestlc.mc1.sav`.
As the output text indicates, you may analyze the output (printing of paramter uncertainties and generation of corner plot and others) using:
`ufit_mcshow,'fitestlc.mc1.sav'`
and generate a plot with the model lightcurve's 1-sigma uncertainty relative to the best fit with
`ufit_limitshow,'fitestlc.fout1.utm','fitestlc.mc1.sav'`

The further study of some of the setup-files in conjunction with the UTM and UFIT manuals (`documentation/utm_doc.txt` and `ufit_doc.rtf`) is then recommended.

## Citing UTM or UFIT
The preferred way is through its [entry in the Astrophysics source code library](https://ascl.net/1412.003):
 `Deeg, H. 2014, Astrophysics Source Code Library, record ascl:1412.003`
  
Entry in ADS: [ `2014ascl.soft12003D` ](https://ui.adsabs.harvard.edu/abs/2014ascl.soft12003D) 

### Acknowledgement
The authors thanks Prof. Jean Schneider, who proposed a transit modeller and fitter for multiple objects as a project for a stay of the author at Paris Observatory in 1999, where the initial version of UTM was written.  

### Repository of earlier UTM/UFIT versions
UTM/UFIT versions until 20210327 (2021 March 27) are available at
[ftp://tep:fu9dufa5@ftp.iac.es/idl_hans_lib/utm/](ftp://tep:fu9dufa5@ftp.iac.es/idl_hans_lib/utm/)
The UTM/UFIT distribution on that ftp server will not longer be updated and UTM versions posterior to 2021 March 27 will only be available here on github.




