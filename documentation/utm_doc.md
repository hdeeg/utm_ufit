# UTM
A universal transit modeller

H.J. Deeg, 27 Mar 2020  (version history at end) 
UTM Version described:  utm.pro of 18 Mar 2021 


## 1. General Description
———————————

UTM is a modeller for all kinds of configurations between any number and kind of objects, such as stars, planets, moons,..  
It calculates the resulting brightness, and displays graphical representations of the transit-configurations. For an overview see the paper by Deeg 2009 that is included in the UTM distribution (file Deeg_UTM_proc_Boston08_2.pdf ). For citing UTM, see Appendix at end.

Depending on the settings of the tflag parameter (see 'Setup file-, Sect. 2.1), UTM can be used to:

- generate model lightcurves from scratch (tflag=0),with or without noise (parameter 'onoise') 
- generate model lightcurves for given time-points (tflag=1), with or without noises
- add models of planetary transits to given lightcurves (tflag=2). This option is useful to generate test-data based on light-curves of non-transiting stars.
- subtract models of planetary transits from given lightcurves (tflag=3). This option returns the residual data - model.
- in combination with UFIT (see ufit_doc.txt), UTM can be used to fit models to an observed lightcurve.

UTM allows the definition of three general kinds of objects:
- dark objects (planet, moon, ring)
- luminous objects (star)
- objects that are neither and do not cause eclipses (point). They are intended mainly to define positons of barycenters in hierachical systems.


### 1.1 Functional description of UTM:

UTM will read all required input parameters from a setup-file, which is a text-file with parameter - value pairs.
The setup contains the parameters of the bodies that are simulated and also all generic parameters that deviate from defaults.
For each body (star, planet, ring..), a computational object is then generated, that describes its kind, size in the sky-plane, orbital parameters, actual position, together with type-specific parameters, such as limb-darkening coefficients for stars, or the inclination for rings. In the default operational mode of UTM, a pixelized 2D representation of the transmittance of each object is also generated, and bright objects (emitters) get a 2D representation of their luminance towards the observer. 
The objects' orbits are prescribed by Keplerian parameters that are supplied from the setup-file, with defaults that simplify their use. Keplerian parameters are used since they are independent for each orbiting object. More interdependent parameterisations can be obtained through user-supplied 'preprocessor' codelets, which -for example- may replace a planet's orbital inclination by its impact parameter on the central star. UTM is not dependent on the use of any specific units for length, time, luminosity, but units may be declared as input parameters. This is useful for the labelling of graphical output, but it permits also the use of relative units. For example, we may replace a planet's absolute orbital semimajor half-axis 'a' by its relative one, a/R*, by setting R* (the central star's radius) as the unit for length. With the appropriate unit choices and preprocessor codelets (examples are given in the /example directory of the UTM/UFIT distro), UTM can be adapted to any of the common parameterizations used in the simulation of planetary transits or stellar eclipses.

In UTM's main loop, the program will step through a series of time values or epochs. These epochs can be red from an input file (usually a light-curve for which a simulation is to be produced), but equidistant epochs can also be generated internally from parameters in the set-up file. At each epoch, the objects' positions in a 3D coordinate system are calculated from their Keplerian orbital parameters. If  2D representations have been generated, they lie in the sky plane. UTM checks then if any of the emitters is overlapped by any other object. If occultations occur, they are broken down into a hierarchy of mutual occultations between two objects, where the one in the background is always an emitter. UTM starts with the mutual occultation in which both the emitter and the occulter are closest to the observer, and calculates the remaining brightness. If the same emitter is occulted by more objects, remaining brightnesses after each mutual occultation are calculated, until arriving at the occulter that is closest to the emitter. If there are further emitters in the system, UTM calculates the occulted fluxes for the next emitter at increasing distance from the observer, starting again with the occulter that is closest to the observer. Following this procedure, multi-body occultations of arbitrary complexity can be resolved. 

The calculations of the mutual occultations' brightnesses are performed either in a pixelized or in an analytical mode, depending on the circumstances (described in detail in Sect. 2.4):

- Pixelized mode: The overlapping areas between the background and foreground object are extracted as rectangular sections from the 2D arrays of their luminance respectively their transmittance. These sub-arrays are scaled to a common pixel-resolution, with the coarser (smaller) one of the two being re-binned to the finer one, so that both will be 2D arrays of identical size. The remaining flux in the overlapping area is then calculated, and from that the remaining flux of the entire background object. If the background object suffers a further occultation (or if a display of its 2D image is requested), the objects entire 2D luminosity distribution is also updated after each mutual occultation.

- Analytical mode:  Given a set of conditions (a bright object is occulted by only one object; both must be of circular shape; no 2D graphical display; not all limb-darkening laws admitted; see Sect 2.4 for details), the resultant brightness from a mutual occultation is calculated using Mandel & Agol (2002) algorithm, as implemented in the 'exofast_occultquad.pro' routine from the EXOFAST package by J. Eastman(2013). The code 'exofast_occultquad_cel.pro' from EXOFASTv2 (Eastmann et al. 2019, and described in more detail in Agol et al. 220, AJ 1259,123) can be used alternatively, but is slightly slower, see 'occultcelflag'. The analytical mode has generally a significantly improved speed over pixel-mode and a precision that is not limited by the coarseness of the pixel-arrays.

By default, UTM will use a hybrid mode for the calculation of transit fluxes. The analytical mode is then used whenever permitted by the circumstances, and else, the pixelized mode is used. This behaviour can be overriden, see 'anlycalcflag' below.


UTM implements several limb-darkening 'laws': the linear, the quadratic and the root-square one. For the quadratic one, several specifications of the coefficients are possible, see entry to 'limblaw'

The orbits for all objects can be freely positioned in a three-dimensional coordinate system, and a full set of circular or elliptic Keplerian orbital parameters can be given, as described in Sect. 2.2. Arbitrary sets of input parameters can be ingested with preprocessor conversion routines, which are explained in Sect. 2.3. This permits, for example, the use of the stellar mass and planet period as input parameters, with the preprocessor calculating the semi-major axes using Kepler's law.

Plots of the system luminosity and animated displays of the occulted objects and of orbital positions can be shown in a customisable way.

### 1.2 System requirements

IDL 8.0 version 8 is required for the current version.

The minimum source files needed to run UTM are:  
  utm.pro  
  utm_ancil.pro   

The .pro files have to be in the current directory, or be in directories that are discoverable for IDL through the  IDL_PATH.


## 2. Use of UTM
 
Start UTM from the IDL command line with:

`IDL> utm,'setupfile',dflag`
example: `utm,'6bodies.utm',3`

where setupfile is the name of the setup file (this file is described
below).  All initial parameters of UTM are red in from the setup file, except for `dflag`. Some more command line-parameters are available to UTM which are intended for calls to UTM within wrapper programs (e.g. ufit.pro) ; they are described in the Appendix. 

`dflag` may take one of the values 0-4 and controls the amount of information that is displayed:

0: *None* (only the resulting lightcurve may be written to a file)

1: *Printout of the lightcurve:*  
The time and brightness values of the modeled lightcurve are 
	written to the terminal. The brightness may be displayed in three 
	different ways, depending on 'oflag' (see description of setup 
	files)

2: *Plot of the lightcurve* 
A plot of the last 100 time-points during UTM exection is shown. Its symbols are similar to the final light-curve plot described below. 
	Some warning messages are also shown only if dflag ≥2.

3: *Displays of bodies and orbits of the simulated system*  
	The object's orbits are shown in three panels, with the objects in the x-y plane (the view as seen from the observer), the y-z plane (see from the 'side' with the observer to the right) and in the x-z plane (the view from 'above' onto the system). In the y-z and the x-z plane, the observer is indicated by a triangle. Also, for each luminous object, its current surface brightness distribution and eventual occultation is shown. The object's size is always scaled to fill the panel (Hence two stars of different size appear equal). 

4: *Single step mode*   
	After each time-step, the program stops (type .cont to continue) and the current positions 
	(in X,Y,Z coordinates) of the objects are printed, together with
	a list of the objects which are transiting.

Note that there are setup-parameters (such as `plotorb`) which may override `dflag`, in order to turn on or off some specific type of displays. 

Independently of the setting of `dflag `, a final lightcurve is displayed in an interactive plot-window (this can be suppresed with  `plcflag -1` in a setupfile). 

The symbols this plot are  (these depend on the `tflag`parameter, see below):  
		- solid-line: transit-model (UTM-output for tflag =0,1)  
		- crosses: input light-curve (only for tflag=2,3)  
		- diamonds: residuals from 'data minus model' (for tflag=3) or the sum of 'data plus model' (tflag=2). If oflag=0 (output in system-luminosity) , the residuals are shifted to the off-eclipse luminosity in order to overlap with the data outside of eclipses. 





2.1 Setup files
---------------

The setup file defines completely the initial state of UTM. In it, all
the objects and their parameters are defined.  Examples of setup files
can be found in the UTM distribution, see files ending with *.utm. (Note: the older .set extension works as well. *.utm files may however not work with older versions of UTM). Each line in a setup
file consists of two tokens which are separated by one or more white spaces:

   keyword  value  

If three or more tokens are present, they are ignored (to be used for comments)
For the creation of a setup-file, the following rules apply:

- the order of the parameters is without relevance 
- each keyword has to be unique (in case of duplication, the value of
the first occurrence will be used). An exception to both rules are the 'fit' and 'scale'
keywords, which are only used by UFIT (see its documentation).

- between the keyword and the value has to be at least one white space
  (thus, keywords cannot contain spaces. Use underscore)
- alignments and numerical formats can be varied freely 
- unrecogized keywords are ignored by UTM
- keywords should not be repeated. If done, the first value is used .

- the third and further words in a line is ignored
- everything behind a # sign is ignored
- lines beginning with # are ignored
- very short lines (less then 3 characters) are ignored

- Parameters for the individual objects have to be prepended by a one-digit
object number N, following these rules:
- There has to be an object with number 0 (zero)
- The  object numbers must to be consecutive ( no 'holes' in the numbering). 
(If an intermediate object with number N is to be removed, it can be replaced by an 'Ntype point' statement, which will generate a zero-sized point instead.  )
- The object-number of a body whose position is pending on a main-object needs to be higher then 
the main-object's number


Below is a list of the parameters that may be defined in a
setup-file for UTM. Given is the name of the parameter, and the
variable-type of the associated value, followed by the default.  For the generation of new setup-files, consultation of examples in the UTM distribution is recommended


Notes:
 - bin is a binary with possible values 0 or 1.
- strings may not contain white-space or need to be in quotes.
- {A,B,C} or {0-1}means that values have to be from the set or range in {}. 
- the value after the comma is the default, which is used if the parameter is not given.
- only the major paramters (radi, period, size) give a warning if they are missing; others are set to default without warning.


#General parameters  
#unit settings. Are only used for display purposes
sunit  string,Rsun    	#units for size
tunit  string,day     	#unit for time
lunit  string,Lsun    	#unit for luminosity (see also oflag)


tflag	{0,1,2,3},0      #0: generate time-points from tinit,tinc,tfin
			 #1: use time-points given in file 'indatafile
			 #2: add model to data given in file 'indatafile'. 
			 #3: subtract model from data given in file 'indatafile' (generate residuals). 

			 #tinit,tfin,tinc are only used if tflag=0
tinit  double,0           #if tflag=0: time at start of UTM (in units of 'tunit')
tfin   double,10      	#if tflag=0: time at end of UTM
tinc   double,0.1           #if tflag=0: step-size of time

resample int,1  #if larger then 1: Over times given by ‘sample_span’, the light curve is calculated 
			#on a resampled time scale that is finer by the factor ‘resample’, and that is then
			#rebinned (unless keep_resample=1, see below) to the original time-points.

sample_span float,0.02043  #the time-span over which the lightcurve is being resampled, in units
			# of ‘tunit’ . Ignored unless resample > 1. The default is the integration time of
			# Kepler long-cadence data, valid if  tunit = day (else derive from 29.4244min)

keep_resample bin,0  #if 1, the UTM output will be the resampled (higher time-res.) lightcurve. 
			#If 0, output will be rebinned (averaged) to the input time-points. Ignored 
			#unless resample > 1. 
			
indatafile string  #name of ascii table with time or time and flux values. Required if tflag>0. 	
			#NOTE: alternative parameters 'inlcname' or 'infilename' are also accepted
			 #if tflag =1: file needs to contain time-points in first row. The rest is ignored.
			 #if tflag =2 or 3, a lightcurve needs to be given. Time-points in first row, brightness 
			 #values (in same unit as given by oflag and/or lunit) in second row, rest ignored.

inlcdigit  int       #if  present, then leading digits of time-values in indatafile are stripped. 
			#inlcdigit indicates the number of remaining digits before the decimal
			#point. This is useful if light-curves with full Julian Dates are given. For example,
			#with ‘inlcdigit 1’ , a time-stamp like 2456789.123 will be
 			#stripped to 9.123. Any parameters indicating epoch values need to be
	 		#compatible with the stripped time-stamps.


anlycalcflag  {-1,0,1},0  #This flag  determines the mode in which fluxes are calculated during transits. 
			 #-1 : fluxes are always calculated from pixelated representations of transiting objects. 
			#0 (default): fluxes are calculated analytically whenever possible,  and through pixel representations elsewise.  
			#1: fluxes are always calculated analytically. UTM will stop if conditions for analytical calculation are not met (see below)

occultcelflag  bin, 0   #0 (default) For the analytical calculation of transit fluxes, exofast_occultquad.pro from EXOFAST is used.
			#1 exofast_occultquad_cel.pro from  EXOFASTv2 is used. This should have a calculation of higher precision near limbs, but is approx. 15% slower.
			
gsize   integer,24  #size of 2D arrays that are used to represent objects. 
			#Minimum recommened size is 12. gsize cannot be changed 
			#without restarting IDL! 

gsizemultstar  float,4   #multiplier to gsize for stars  
gsizemultplan  float,1   #multiplier to gsize for planets  
gsizemultrung  float,2   #multiplier to gsize for rings

finestep   bin,1       #finestepping: if set, overlapping regions in arrays of occulter or star are being shifted by 
			#subpixel lengths, using interpolation, to increase precision. 


preproc  filename,0   #program to prepocess setupfile (without .pro); 
			  #outcomment or 0 for none (See Sect.2.3)
saveprepset     bin,0    #1: save preprocessed setupfile, 0: don't

limblaw {lin,quad,pmquad,usquad,root},lin  # limbdarkening law of luminous objects. The options are as follows, with
				#limbd and limbd2 being input parameters in setup-file (see below):
  				#lin:  linear law, I(μ) = I(1)(1 − u(1 − mu))  ; limbd :=u  
				#quad: quadratic law, I(mu) = I(1)(1 − u1(1 − mu) − u2(1 − mu)^2);  limbd :=u1 , limbd2 :=u2
				#pmquad: quad. law with 'plus-minus' coefficients (Gimenez 2006, A&A 450, 1231): limbd := u+ = u1 + u2; limbd2 := u- = u1-u2
				#usquad: quad. law  with uniformly sampling coefficients that avoid unphysical configurations (Kipping 2013, MNRAS 435, 2152):
				#        limbd := q1 =(u1+u2)^2;    limbd2 := q2 = 0.5u1 / (u1 + u2) ; inputs of q1 and q2 outside allowed range [0..1] are clamped into it.
				#root: sq.-root law (Diaz-Cordoves & Gimenez 1992, A&A, 259, 227): I(mu)=I(1)(1−u1(1−mu)−u2(1−sqrt(mu))); limbd:=u1 , limbd2:=u2, 
                                #where I(1) is the intensity at the disk-center; mu = cos(gamma), with gamma being the angle between line 
				#of sight and normal to stellar surface; u1 and u2 are the LD coefficients; 
				#For the lin. and quad. law: LD coefficients leading to unphysical conditions (see Kipping 2013) will cause the star to take the shape of an 'X', and 
				#a warning is given. The motivation for this behaviour is as follows:  For unphysical LD values encountered during a fitting or MCMC run (with UFIT),
				#the X-shape will generate a poor transit model, and the fit will be directed away from these bad values.

oflag  {0,1,2},0         #0: output/input is given as system-luminosity, which is the summed luminosity of all objects 
			 #1: output/input is given as relative flux variation dF/F, normalized to dF/F=0 outside of eclipses.  dF/F = (L - L0)/L0 , where L0 is the total 
		               #luminosity of the system in the absence of eclipses, and L the actual luminosity. Transits have negative dF/F values. 
			       #dF/F +1 gives the normalized luminosity L/L0.
			 #2: output/input is given as magnitude variations, normalized to delta-mag=0 outside of eclipses. This is dF/F converted to delta-magnitudes. 
				#Positive delta-mag values mean less flux.
 
			 #oflag affects only the graphic display and the input/output-lightcurve. Values for the contlum and Nlum 
			 #parameters (see below under ‘stars’) always have to be given in the flux-unit defined by ‘lunit’.

onoise float,0    #Gaussian noise is added to output lightcurve, with a standard deviation given by onoise, in units like oflag/lunit 	

contlum float,0  #adds additional (contaminating) non-transiting 'third-light'
			#to system luminosity, in units of lunit






wflag  bin,0      #1: write output model-lightcurve to file, 0: don't do it
outlcname string,model.lc  	#name of model-lightcurve 
                               #(NOTE: the depreciated parameter 'lcfilename' is also accepted)
lchead bin,1      #write setupfile as header of output curve (only if wflag =1)

#The following  parameters add constant, linear or quadratic offsets to the output data, using as unit whatever setting of lunit and oflag (system-luminosity, dF/F, delta-mag) is given. These offsets are added as a final step, after all occultations have been calculated in UTM, and after the conversion to the final type of output (given on oflag) has been done.  If y(t) is an output-value at a time ’t’ before offsets are applied, the output  with the offsets will be:
y’(t) = y(t) + ooff + oslope x (t - tozero) + oquad x (t - t0zero)^2
These offsets are intended to permit f. ex. the simulation of slopes from airmass effects. If a slope is being fitted using UFIT, either tozero or ooff should be free parameters, but NOT both (else the fit is overdetermined).

ooff   float,0     #adds an offset to output after sim, in units like oflag/lunit
oslope float,0 #adds a linear slope to output, units are lunit per time-unit
oquad float,0  #adds a quadratic offset to output, units are lunit per (time-unit)^2
tozero double,0  #reference time against which oslope and oquad are calculated. If not given, the first time-point is used.

#Display flags to show or suppress graphic output.
#Flags with three options {-1/0/1} behave as follows: -1 always supressfes; 0 uses default pending on the value of dflag; 1 always shows the corresponding output. Giving a value of 0 is similar to not setting the flag.
 
plcflag         {-1,0,1},0   #suppresses/defaults/shows final interactive plot of model.
plcsimflag      {-1,0,1},0   #suppresses/defaults/shows last 100 pts of model during simulation
plotorbflag     {-1,0,1},0   #suppresses/defaults/shows plots with objects' orbital positions
orbtrflag          bin,0     #1: leave traces in orbital position plots., 0: plot only the objetcs current positons
dispstarflag    {-1,0,1},0   #suppresses/defaults/shows figures of the stars
xpflag  bin,0    #if 1: makes final plot of x-positions with double-lines that indicate the size of the objects.  (this special plot is not generated by any setting of dflag)

#t for ALL object-types (point, planet,star, ring) the following parameters are recognized:
(N is the object number)

Ntype  {star,planet,ring, point, extocc (,moon)}  #type of object
	                              #point: point-like object (of size 0) that does not cause any eclipses. Use for reference positions in hierachical systems  
				      #extocc: external occulter (see Appendix F)
				      #moon: depreciated, is identical to planet.
		 
  
Nperiod  double,0  		#period; if 0.0 then is fixed. In units of 'tunit'
Ntrepoch    double   	#epoch of transit (see 'coordinate system') 
Necepoch  double   	       #epoch of eclipse (see 'coordinate system') 
Nphase    double            #phase (for defin. see 'coordinate system') 
Nepoch    double,0          #epoch of periapsis (see 'coordinate system') 

Nrdist   double,0     		#orbital semimajor axis, in units of 'sunit' 
Ninclin  double,90       #inclination of orbit in degree; 0=face-on, 90=edge-on
				    #(for defin. see 'coordinate system')
Necc   {0-1}float,0             #orbital eccentricty
Nomega  float,270    #arg. of periapsis or periastron angle in degree (see 'coordinate system') 
Nsqrtecsin  float    #sqrt(ecc)*sin(omega), alternative overriding params for ecc, omega 
Nsqrteccos  float    #sqrt(ecc)*cos(omega)

Nomegadot  float,0    #change of omega in degree per time-unit (see 'coordinate system') 

Npa      float,270      #long. of ascending node resp. position angle in degree (see 'coordinate system')
Nmainobj  {0 to (N-1)},x    #number of associated main object. If set, the body's positon is calculated relative to the main object's position, whose number N must be lower. If Nmainobj is absent or set to 'x', the body's position is calculated relative to the origin of the coordinate system. 

#for object-types  'planet', 'star' , the following parameters are recognized:

Nradi  float,1         		#(equatorial) radius of object , in units of 'sunit'
Nmass float,1         		#mass of object  (not used in current UTM version)
Noblate float,0        #oblateness = (radi - Rpolar)/radi, were radi is the equatorial radius
Nspinpa float,0        #pos. angle of the spin-axis in degree, from 0=North (+y) to East (-x). 
				#Only used if oblate <> 0.0
	  
#for object type 'star' additionally:
Nlum   float         		#luminosity of object, in units of 'lunit'
Nlimbd {0-1},0.3	      #first limb-darkening coefficient, for stars only
Nlimbd2 {0-1},0.3      	#second limb-darkening coefficient,  not needed for linear law.

# for object type 'ring' some parameters are additional, or have different meaning :
Ntransp   {0-1},0.5      #transmittance (transparency) 0=opaque,1=transparent
Nradi      float       		#is the outer radius of the ring
Nrdist     float                #is the inner radius of the ring  
Ninclin   double,90         used for ring-plane in same sense as orbital planes for stars,planets
Nspinpa double,0          used for ring-plane in the same sense as spin-axes for stars,planets

Notes: 
- After changing the gsize or the gsize-multiplier parameters in the setupfile, IDL needs
to reset (use .reset).  Otherwise, an error 'Conflicting data structures' will occur.
This problem is due to IDL keeping structure definitions between
program invocations.
- Starting with UTM version 5may2015, the object type 'moon' has become obsolete (since all body-types permit now the definition of a main-body) and is now identical to 'planet'. It is kept for compatibility with older setup-files.


2.2 Coordinate system
----------------------
The following coordinate system is used internally in UTM:

        y	
        ^
        |
        |
        |_____>x
        /
       /
      z

This is the standard for 3D image processing. The observer is at a
positive z distance. The x-y plane corresponds to the sky-plane and North is the top(+y) and East at the left (-x). Transits in front of the coordinate center occur at (0,0,z). 

The orbits are parametrized as follows:

By default, orbits are calculated relative to the coordinate center at (0,0,0). If a main-object is defined, then the orbits are relative to that object. UTM performs no checks if an orbit is physically viable.

'rdist' is the semimajor axis  of an object's orbit. It is the radius on circular orbits. 

'ecc' is the eccentricty. If not specified, a value of zero (circular orbit) is taken. Negative input values (which might occur during MCMC chains) will be converted to their absolute value. This avoids 'edge-effects' in MCMC runs (see also appendix to Eastman+ 2013, Exofast paper)

All the following angles are to be specified in degrees:

'omega' is the argument of periastron in degree. It rotates the orbit in the direction of the body's motion. If not specified, a value of 270deg is taken. (see below for comment on choice for this default)

'sqrtecsin' and 'sqrteccos' provide an alternative parameterization of eccentricy and arg. of peristron by specifying sqrt(e)*sin(omega) and sqrt(e)*cos(omega). This parameterization provides better coverage of phase-space in fitting routines (see e.g. Ciszmadia et al, 2020, The Transit and Light Curve Modeller). Also note that the eccentricity cannot becomce negative in this case.  
To convert sqrtecsin and sqrteccos back to ecc and omega: With utm compiled, use:
IDL> print,e_omega(sqrtecsin, sqrteccos) 
which prints ecc and omega.

'omegadot' indicates the change of omega in degree per time-unit. Its default is 0. This parameter is permits the simulation of apsidial motion. For each time-point, UTM calculates a value omega(t) = omega + omegadot * (t-t0), where t0 is the reference epoch of periapsis. Note that this equation also applies if 'phase', 'trepoch' or 'ecepoch' are given for inital conditions (see below), since internally UTM converts them into an epoch of periapsis (which can be displayed running UTM with dflag=4)

'inclination' is the angle between the sky-plane (x-y) plane and the orbital plane in direction of the orbital motion, measured at the ascending node on the side of the reference direction (North or +y, but rotation by position angle may change this). If not specified, a value of 90deg is taken (edge-on orbital plane). 

'position angle' is the longitude of the ascending node, defined as an angle in sky plane between the +y direction (North), and the nodal line (towards the ascending node) . From +y (North), it goes counterclockwise to -x (or to East). For pa=0, the ascending node lies therfore on the +y axis and the nodal line from descending to ascending node points from S to N (-y to +y). If not specified, a value of 270deg is taken (see below for comment on choice for this default). Note: for rings, pa indicates the counterclockwise rotation of the normal vector of the ring-plane. As before, from +y (North), it goes counterclockwise to -x (or to East)  

The ascending node is the point in an orbit were the object crosses the sky-plane moving away from the observer. This corresponds to the point where an object's z coordinate changes from a positive to a negative value.  
 
Intial conditions (the locaction of the body on the orbit at the begin of simulation) need to be set with either of 'trepoch','ecepoch','epoch' or 'phase' parameters, defined as follows:

'epoch' indicates a moment in time when the body passes through its periapsis. The body's initial position is then calculated for the first time-point that is to be simulated (e.g. for tflag = 0, this is 'tinit'). Internallly, UTM uses the epoch as the fundamental value for the initial location. 

'phase' is a value from 0…1 that indicates for the first time point the fraction of an orbital period (in time) that has passed since the previous epoch.  Phase is identical to the Mean anomaly divided by 2 pi. If phase is given, a corresponding value for the epoch will be calculated.

'ecepoch' is the moment when a body passes through the 'eclipse position' behind the coordinate center (e.g.  if i~90deg, it gets eclipsed by the coordinate center), crossing the xz plane at y=0 with negative z-coordinates. Given that this happens at nue=1/2 pi - omega, where nue is the true anomaly, UTM converts this to the epoch of periapsis.

'trepoch' is similar to 'ecepoch', indicating the moment when a body passes through the 'transit position' in front of the coordiante center (e.g. it transits in front of the coordiante center if i~90deg), crossing the xz plane at y=0 with positive z-coordinates. Given that this happens at nue=1.5 pi - omega, where nue is the true anomaly, UTM converts this to the epoch of periapsis.

For all initial condition parameters, the internally used periapsis epoch is being displayed if UTM is executed with dflag=4.
If several of the initial condition parameters are given, they override from highest to lowest preference:  trepech, ecepoch, phase, epoch. 

Default values: The reason that both omega and the positon angle have default values of 270 is as follows: For inclinations near 90deg,  the body will then be in the 'transit' position (near 0,0,z) at times that correspond to the 'epoch'  or to phase =0, independently if being on a circular or eccentric orbit. Furhtermore, near phase=0, it will then move towards the right (to +x) to the ascending node (away from observer). With these defaults and for inclinations <90 deg, transits will be across  the lower half of a body placed  at the coordinate center.

Old setupfiles (prior to 8May2015) on version UTM8May15 and posterior: 
Setup files created for previous versions that contain position angles or the 'eccplanets' type will give different results. Setupfiles without these should run with identical results, except for a reversal of y-axis positions (transits for i < 90deg will now be below the coordinate-center).

For positon angles, this is due to the different definitions of the positon angle and the 'epoch'. To correct old setupfiles, add 180 to all position angles (if they indicated pa =90, which was the old default value, it is better to delete them, as the new UTM will then default to pa=270deg).

For eccplanets: this class has been removed and should be replaced by 'planet'. Values for omega, epoch, phases will however need a deeper revision, since the orbital orientation of the old 'eccplanet' class was not well defined. Also, for the 'planet' class, the rdist parameter needs to be defined (for the old 'eccplanet', this was calculated from mass and period). A prelimary fix to change eccplanet parameters to those for 'planet' is to: 
-change the sign of the inclination
-change the sign of the positon angle
-subtract 90 from the positon angle

 
2.3. Preprocessor programs for arbitary sets of paramters
-------------------
Sets of parameters different to the UTM standard ones (e.g. different to the orbit-paramters of Sect. 2.2) can be used in setup-files, if a corresponding 'preprocessor' routine is provided.
This permits, for example, the use of the stellar mass and planet period as input paramters, with the preprocessor calculating the planet-distance 'rdist' using Kepler's law.  The example-setupfile 'keplertest.utm' invokes kepler_prep.pro to perform exactly this operation. Another preprocessor is used in example-setup 'K446sp.utm', where the input parameters Rpl/Rst (planet/star size) and a/Rst (semi-maj. axis /star-size) are used to calculate the absolute planet size and distance used by UTM, by running 'rpl_apl_prep.pro' .The preprocessor is invoked before running the UTM main routine, which then uses the temporary setup file generated by the prepocessor.

The standard scheme for the use of a prepocessor is as follows:
In setupfile, these lines should be defined
preproc  filename   #program to prepocess setupfile (without .pro)
saveprepset     1    #1: save preprocessed setupfile, 0: don't

The preprocessor, named 'filename.pro', needs to be in the directory of the setupfile, in the directory where utm.pro resides, or anywhere else in the IDL-PATH.  The use of 'saveprepset' is optional, setting it to 1 permits revision of the intermediate setup-files (file-names ending with *preprocd.set).   

The preprocessor is written in usual IDL code, but its layout should follow the following simple scheme:

pro example_prep    ;preprocessing code for UTM
	;read parameter's values from setup-file
	inval1= double(parasign('1inpar'))  ;1inpar and 2inpar are parameters in setupfile
	inval2= double (parasign('2inpar'))  
	;perform some operation on these values
	1outval= inpar1 / inpar2   
	;save the resulting value into parameter '1outpar'
	paradd,'1outpar',1outval
end

parasign reads a setupfile parameter and paradd adds one.
For more detail on these commands, see code in setupfile.pro




2.4 Hybrid versus analytical versus pixel-based calculation mode
--------------------------------------------------------
The default hybrid mode ﻿(anlycalcflag= 0) combines the generally better speed and precision of the analytical calculation with the flexibility of the pixelated calculation, and therefore is able to resolve any transit-configuration. It can therefore always be used. The analytical-only mode ﻿(anlycalcflag= 1) is recommended if it is certain that it can be used throughout the entire simulation. It offers usually a small gain in speed over hybrid mode (this gain can be large in the case of large gsize values and short lightcurve), since the initialization of the pixel-objects' 2D arrays is suppressed. The pixel-only mode ﻿(anlycalcflag= -1) serves to reproduce results from older UTM versions and for tests against analytical calculations. If an entire simulation has to use pixel-calculations (e.g. for never meeting the conditions for analytical calculation), there is no difference in speed between hybrid and pixel-only mode.
In analytical mode, currently exofast_occultquad.pro from EXOFAST and exofast_occultquad_cel.pro from  EXOFASTv2 are included. The original routine from EXOFAST, based on the Mandel&Agol(2002) algorithm, has been maintained as default (see 'occultcelflag'), due its 10-15% higher speed. The improvement in precision of the EXOFASTv2 version is on very small scales (see Agol+2019, AJ 159,123, Fig. 14). In practice, the flux-differences between these algorithms will be of much smaller scales than the errors from a 2-parameter Limb-Darkening law to describe the true fluxes of an occultation event.

For the use of the analytical calculation, all the following conditions have to be met:
- The 2D display of the occulted stars needs to be turned off (see dispstarflag). This is usually the case if dflag 0,1 or 2 is used. If the use of dflag >= 3 is desired ( e.g for printing of detailed orbital info), the 2D display needs to be turned off explicitly by setting dispstarflag -1.
- The objects need to be stars or planets with circular shape (obliquity =0) or 'points'. Rings and external objects are not supported.
- The the square-root limb-darkening law is inadmissible. All others are fine.
- Only mutual occultations between two objects are supported. However, systems with more than two objects may be simulated purely analytically, as long as only two of them transit each other at any time. In hybrid mode, if simultaneous transits among several objects occur, pixelized and analytical calculations may occur concurrently.



3. Precision of UTM
-------------------

Below, results are shown from tests using setup file 'setcalib', where
a star is transiting in front of a uniformly bright disk. The tests
were performed for two sizes of 2D arrays (setup parameter gsize).

gsize=12

Rplan/Rstar	theor.loss	loss in UTM	rel. error
0.001		1e-6		0.000001	<0.5 (limited by display)
0.01			1e-4		0.000101	<0.015
0.1			0.01		0.010050	0.005
0.2  			0.04		0.040200	0.005
0.25			0.0625	0.062810	0.005
0.5			0.25		0.251258	0.005
0.9			0.81		0.824226	0.017

gsize=24

0.001		1e-6		0.000001	<0.5 (limited by display)
0.01			1e-4		0.000100	<0.005
0.1			0.01		0.010000        <0.00005
0.2			0.4		0.039998        0.00005 
0.25			0.625	0.062496	0.00006
0.5			0.25		0.249989	0.00004 
0.9			0.81		0.810408	0.0005


TIP:

If modelled lightcurves display significant steps, the gsize parameter
should be increased. Steps result from a time-increment that
is so small that the rasterized representations of the objects remain
in the same relative position to each other.


4. Scope of current version of UTM and possible future extensions
-----------------------------------------------------------------

UTM is a modelling program only. 
For the fitting of transit-models to lightcurves, use the program 'UFIT'

UTM  implements circular and elliptic orbits in aritrary orientations. The orbits are indicated by purely geometricl parameters and UTM performs no check on their physical viability. To obtain geometrical paramters (e.g. the semi-major axis 'a') from other parameters like mass and period, a preprocessor-program should be written (see e.g. kepler_prep.pro)

UTM calculates the lightcurve only in a single color (wavelength). If
multi-color models need to be generated, one could run UTM repeatedly
with different values of limb-darkening. A wrapper-program for
multi-color lightcurves was specifically written for simulations
related to COROT and is available on request from the author.

Reflection of stellar light on planets is not taken into account (this
would best be done by creation of a new 'luminous' object-type
'reflecting_planet')

=======
Appendices

A1:  Description of the code for future changes
---------------------------------------------------

Care has been taken to comment the code line-by-line, and here only a
general overview is given. UTM has been written in IDL, using the
object-oriented extensions that are documented in the IDL manual
'Objects and Object Graphics'. This manual (as are all IDL manuals) is
available in the on-line IDL help.

The general layout of the source file (utm.pro) is:
-sub-procedures and functions
-object-methods (functions and procedures)
-main program

A1.1 Main program
----------------
The main program is relatively short, and at the end. It starts at the line:
pro utm,setfname,dflag

An overview about it is given here:
-A few general parameters are read in (mostly related to size of the geometric
 arrays)
-The different types of object structures are defined. 
    tform is a 2D array with image of the transmittance (0=transparent, 255=opaque)
    lform contains an image of the intrinsic luminosity (luminous objects only)
    lformoc is an image of the current luminosity (this one is displayed in 
	the graphic panels when dflag >= 3 )
-The objects are initialized (=filled with data)
-some more parameters are read in to be used in the time-increment loop (the 
 main loop), and graphics is initialized
-the main loop starts (for t=tinit,tfin...):
-the objects are sorted in z-distance (most distant, at -z, is first)
-starting from the most distant object, they are tested if they are a luminous
 object (star), with:  
	if obj_isa(obj[izs],'star') then ...
-if it is a star, a list is made o& objects that are transiting in front of it,
 with:
   	frind= obj[izs] -> findinfront(i,zind,pos,rad)
-if this list is empty, the  star's brightness is added to the total 
 brightness with:
	syslum=syslum+obj[izs] ->getlum() ;no transits
-if this list contains any other object(s), the star's occulted brightness is 
 added to the total brightness with:
         syslum=syslum+ obj[izs] ->gettrlump(obj[frind])
-the remainder of the loop is again concerned with graphical display

The subroutines and methods called by the main program should be
fairly self-explanatory, only the gettrlump() function is
significantly more complicated. 

A1.2 gettrlump
---------------
(For list of arrays used, see end of this section.)

gettrlump is the main routine for pixelized brightness calculation, and 
is the most complicated routine in UTM. 
It returns the brightness of the calling object in the
presence of one or more occulters. The idea is, that the overlapping
zones of the occulter and of the luminous object are scaled to the
same size. The luminosity behind the occulter is then multiplied by
the occulters transmittance.

The operations are done on the self.lformoc array (which is a 'working
copy' of the static self.lform array), and on the occulters tform
array.

The overlapping zones are calculated in the lines defining the
variables lindmin,lindmax,tindmin,tindmax, which contain the indices
of the overlapping subarry from the self.lformoc and the tform array.
These corresponding subarrys are lformbh and tformoc

Also calculated are the scales of both arrays (lscale, tscale) which
are in units of 'sunit' per pixel (sunit is the arbitrary size unit
defined in setup).  The coarser of the two subarrys (lformbh and
tformoc) is then regridded (calls using 'congrid') to the size of the
finer one.

The normal case is that the tformoc array of the occulter is finer
then the lformbh array (since the occulter is normally much smaller).

The program first considers the other case in line
 if tscale*tsize ge lscale*lsize … ,
That should only occur for extreme size differences in eclipsing
binary stars. Here, the tformoc array is enlarged to correspond to
lformbh, giving the tformocproj array. lformbh is then scaled by the
transmittance, and glued back into self.lformoc.

The normal case (tscale*tsize lt lscale*tsize) follows next in line 
 endif else begin  ;the lform array is scaled to tform...

The basic structure is the same, but more care is involved to
obtain a precise calculation. A simple reversal of the procedure above
results in brightness losses that have errors of up to 20% due to the
imperfect size-match between lformbh and tformoc. This match is not
perfect, since these arrays can only increment their sizes by full
pixels, and the scale of the projected array lformbhproj deviates from
the occulters scale. The brightness-loss that is obtained in the
lformbhproj and (after re-sizing) in the lformbh array is therefore
scaled by a correction factor derived from the ratio between the scale
of lformbhproj and tformoc. Errors in the loss calculation are now
orders of magnitude smaller. 

Overview of arrays of gettrlump

self.lform  : luminous object being occulted; original byte-array
self.lformoc  :   same, but array that is being operated on; is also output 

tform:  occulter shape; original size
tformoc: subarray of  tform that is occulting lform

lformbh:  subarray of lform that is behind occulting object (before occultin)
	lumbhorg:  luminosity of lfrombh before transit
lformbhproj:   lformbh scaled to size of tformoc. Occultation is operated on this.
	lumbhprojorg:   luminosity of lformbhproj before transit
lformbh:  subarray of lform that is behind occulting object (after occulting, gets recreated)
	lumbhafter:  luminosity of lfrombh after transit
	desiredlumbh : luminosity that lumbh should have to be correct


Pendning Improvements /revisions:
Better precision of ingress/egress can be obtained, if the oversizing
of the arrays by 0.1 is reduced, to a constant of 1 pixels.  Precision
for binary eclipsing lightcurves can be improved if loss-correction for the
case tscale ge lscale is implemented.


---------------------------------------------------------------------------
Appendix B: Periastron angle

Regarding the periastron angle omega: In the literature, the quoted value is custumarily that of the primary object in a system (e.g. central star in a planet sys) even if that value is meant to apply to a planet. In UTM, omega needs to be defined however individually for each object whose orbit is to be calculated. For a planet in a star-planet sys, this means that omega needs to be increased by 180deg over literature values.

 —————————————————————————————————
Appendix C: Setting up eclipsing 2-body systems.

2 body systems orbiting a common center should follow these rules: 
- inclination, position angle and eccentricity need to be the same
- their omega values must differ by 180 deg.
- i) Either their epochs of periapsis (parameter 'epoch') have to be the same, or ii) the eclipse-epoch 'ecepoch' of one body (usually the primary) should be identical to the transit-epoch 'trepoch' of the other. This, because one body's eclipse is the other's transit. With the epochs defined that way, their eccentricity or omega values can be varied (maintaining them 180 deg apart) without changing the epoch of the primary eclipse (the phase of the secondary eclipse will vary of course). 

The prepocessor template 'bin_prep.pro' takes care of these rules; so that above parameters need to be defined only for one of the bodies. (This is especially relevant if performing fits with UFIT)

For hierachical systems, barycenters of the sub-systems should be defined as the 'point' class and their positional paramters be given. The stellar components of the sub-system  should then point to that  barycenter as 'mainobject'.

 —————————————————————————————————
Appendix D: Further command-line parameters to UTM

The full definition of UTM is given by:
pro utm,setfname,dflag,tarr,lumarr,intrarr=intrarr,modeloffvalue=modeloffvalue

were setfname and dflag have been described in main-text above. The others, intended for use of UTM within wrapper-routines, are:
tarr  : input array with time-points for which UTM calculates it model. If present, it overrides any input-file (indatafile) defined in the setup file. Requires tflag ≥ 1 in setup.
lumarr: input array with data values (used only for tfag =2,3). Output array with model (tflag=1) or data±model (tflag =2,3)
intrarr: output array that is 1 if there is an eclipsing event at a corresponding time of tarr, else 0.
modeloffvalue: output scalar with the off-eclipse value of the model generated by UTM
 
====
Appendix E: Citing UTM

The preferred way is through its entry in the Astrophysics source code library. Its bibcode in ADS is 2014ascl.soft12003D  ; below the corresponding bibtex entry:

@MISC{2014ascl.soft12003D,
   author = {{Deeg}, H.~J.},
    title = "{UTM: Universal Transit Modeller}",
howpublished = {Astrophysics Source Code Library},
     year = 2014,
    month = dec,
archivePrefix = "ascl",
   eprint = {1412.003},
   adsurl = {http://cdsads.u-strasbg.fr/abs/2014ascl.soft12003D},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

===
Appendix F: External occulters

The use objects of arbitrary shapes or 'external occulters' (type 'extocc', as they are generated outside of the UTM core program) is currently in a prototype implementation.
Such objects are long-integer arrays of 100x100 pixel side-length that describe the transmissivity (where 0 is transparent and 255 opaque), and which are saved as a variable named 'transarr' into a .sav file. These .sav files have to be generated before UTM is invoked. For examples of programs to generate such .sav files, and for setups to run UTM with them, see /examples/extern_body in the UTM/UFIT distribution.

Specific setup parameters for the 'extocc' type are

Ntype    extocc     specifies the extocc object type
Nextname filename  #full name of the .sav file with the 2D transmissivity array
Nsizex  float      #length of object in X direction,  in length-units from 'lunit'
Nsizey  float      #lenght in Y
Furthermore, the usual orbital parameters for period, epoch, distance, etc (see Sect. 2.2) can be specified. 
The gsize parameter does not affect the 'extocc' type


----------------------------------
Version history
28/1/1999 first version for UTM 28/1/1999
5/5/1999 version for UTM 30/4/1999. Updated use of dflag, incl. additional params tflag, oflag etc., some rewrite
3/10/2000 inserted gsize into description of setup file
24/8/2005 added description of limb-darkeing laws, some changes in 'scope of UTM', example setupfile
19/5/2006 added definition of position angle to 2.2 (coord sys)
19/6/2011 (a few intermediate versions not included above; ecplanet was added to appendix.). Added tflag=3 description. fixed error with plcflag  (was plcpflag previously)
7jun2013 added description of tflag=3. A few minor edits.
25jun2013 added desc of contlum paramter
4aug2014 minor fixes to eccplanet description
16jan2015 added reference to UTM in ASCL and the overview paper from 2009.
3feb2015 added list of array-variables to  description of gettrlump
4feb2015 added Sect.2.3 description of preprocessor programs for setupfile.
5-12May2015 many changes for actual UTM version, princially due to integration of elliptical orbits
19May2015  further cleanup; removed appendizes with example-setups
21May2015 added omegadot
17Jun 2015 actualized of description graphics display and oflag
9Oct2015 added description of linear and quadratic offsets (oslope, ocuad, tozero)
23Oct2015 minor fixes on preprocessing 
12Nov2015 added dispstarflag  
27Nov2015 added q1,q2 and u+, u- limb-darkening parameters
3Feb2017 replaced transparency by the more correct term transmittance
20Feb2019 removed mention of mass as setupfile parameter.
11Apr2019 modified description of display flags to agree with UTM version of today
12Apr2019 added 1sqrtecsin and 1sqrteccos parameters for UTM version of today
27Jun2019 added finestep param for UTM version of today
3jul2020  removed duplicated description of oflag
26mar2021 adding of analytical calculation and external objects
27mar2021 improved functional description of UTM
 

