;pro utm.pro
;version history
;26/1/99 first working version
;27/1/99 losses are now calculated precisely, writing of lc to output file
;28/1/99 in gettrlump, the case 'if tscale ge lscale' is disabled. This  
        ;results in much better
	;precision for big occulters but ultimately it would be better to
	;re-enable this case and include loss-correction
;29/1/99 re-enabled above case and included correction calc for it.
;27/4/99 opt. reading of lightcurve (tflag); setting output-format (oflag)
;29/4/99 removed paramfile subroutines into separate files, changes to
         ;plotlcutm
;30/4/99 added input param 'tarr' and output param 'lumarr'
;15/11/99 changed case statement for tflag
;7/4/00  added input param 'epoch'
;11/4/00 added one digit precision to output lc
;10/7/00 fixed error that tflag=0 even if tarr is given in command-line
;10/10/00 inserted obj_destroy at end to avoid leakage of obj vars
;5/2/01 inserted position angle 
;21/10/02 lchead param, allows prepending setupfile (header) to output lcs
;24/10/02 ooff param, adds offset to output lc
;9/5/2003 made unused mass parameter optional in setupfile
;22/8/2005 limblaw param selects linear or quadratic limbdarkening
;23/8/2005 tflag=2 now ADDS model (if oflag=2, adds in mags; needs input in
;         mags, if oflag=1, needs input in some flux unit, subtracts 
;         relative flux loss )
;24/8/2005 command line parameter 'lumarr' is now input (if tflag=1)
;         as well as input/output (if tflag=2)        
;26/8/2005 added copyright and license notes
;2/9/2005 fixed bug if period=0 and only epoch was given
;7/9/2005   changed loss correction algorithm for large occulters (tscale > lscale)
;8/9/2005   changed output to exp numbers for dflag 1, oflag 1
;9/9/2005   changed criteria when lform is scaled or tform is scaled
                   ;in occultation calc. Right now, tform is only scaled if total
                   ;eclipses are possible
;19/12/2005 (JMAV )limblaw param selects linear, quadratic or root-square limbdarkening
;22/3/2006 (HD) made hook: if a negative radius is given it is set to
;               1e-6, and a warning is printed.
;22/4/2008 (HD) first adapatation to eccentric orbits using the rKeplersolve routine by Achille Nucita. New object type eccplanet; still needs testing
;24/6/2008 (HD) minor changes regarding the xy and xz plots
;18/7/2008A  option 3 for tflag which subtracts model
;18/7/2008B   The next two 'B' changes were once inserted into UTM (18 and 
;             19/9/2006) but then development was made based on 22/3/2006 vers.
;             Now they have been reinserted again.
;18/7/2008B (HD) inserted intrarr output keyword, which returns an
;                array indicating the points that are on or
;                off-transit. Note that any point were the luminosity
;                sum varies is considered on-transit; which happens
;                once the box of the planet touches the star.
;18/7/2008B (HD) made transparency of a ring independent of its inclination
;21/7/2008 (HD) fixed error in relative flux addition to input lc (if oflag=1, tflag=2 or 3)
;19/3/2009 (HD) in text-output of lightcurves increased precision of brightness
;15/11/2009 (HD) clarify eccplanet input parameters: rename tau to omega.
                                ;older setupfiles need to reflect this cahnge to work!
;16/11/2009 (HD) possible input in eccplanet is now period OR rdist (=semi-major axes); rewrote rkepler_solve to work with simpler input parameters.
;30/4/2010 (HD) added  xpflag at program end insert optional plot of elongation (x-position); similar added plcflag to plot entire lc at end
;7/5/2010 (HD) added first version with optional preprocessing macros for setup-files
;26apr2011 (HD) fixed bug in  body::init, inserted in post 20080721 version, that generated NAN positons for fixed bodies in e.g. setsimple
;18jun2011 (HD) added some comments into code, no functional change.
;7jun2013 (HD) Returns now graphics sys variables !p,!x,!z to intial status after execution. 
;25jun2013 (HD) enabled tflag=2,3 for oflag=0 (relative flux loss) by multiplying/dividing input data with model. Added 'contlum' for additional (contaminating) light-source, in units of lunit. Fixed error on saveprepset parameter being ignored.
;3feb2015(HD)  revised corrective scaling in gettrlump. The entire output array lcformoc (not only the occulted part) is now corrected. Brightness values are identical to before, but visible squares have gone away.
;5may15(HD) added mainobj parameter to all object types (supressed internal subbody type).
;6may15(HD) added point object class of zero-sized reference points
;7may15(HD) added eccpoint class with defintions of orbital elements omega, inclination, pos angle (= longitude of asc. node) following official conventions. As consequence, transits with incl <90 are now below the star. Added orbtrflag and array posarr with all positions during sim. Added z-y panel panel to plotxy. Default positon angle  was changed to 270 deg, so that default transiting systems keep moving left ot right along +x. To correctly interpret circular orbits  using this default, reversed direction of orbital motion in point::incrempos.
;8may15(HD) all body classes are now based on eccpoint and permit eccentric orbits. Numerous minor changes
;9may15 renamed plotlc to plotlcutm to avoid naming conflicts. Added plot of data vs model for tflag 2,3 
;10may15 added ecepoch parameter. Deleted self.phase and obslete eccplanet-related code. Merged eccpoint class into point.
;13may15  modified paramter-names inlcname and outlcname for input and output lightcurves. Old param-names still work.
;17may15 added syslum0 output-keyword for off-eclipse system-luminosity (needed in ufit to plot correctly residuals)
;21may15 added omegdot input parameter
;29may15 fixed an error that had broken the ring-type. Added oblate and spinpa parameters.
;2 june15 changed default for oflag to 0 (=sys lum) and gsize to 36. Fix of scaling in plotlcutm. 
;8 june 15 a serious precison issue for small planets was identified in the feb15 revision of gettrlump. Inserted switch feb15flag to turn on/off the feb15 version of the correction. Is set to 0, so that the prior correction is again active. 
;11jun 15 added resample and related keywords for resampling and posterior rebinning of data.
;17jun15 revision of plotlcutm. oflag 1 (flux variation) produces now negative values during eclipses
;12sep15 added several display flags (plotorbflag etc. ) that override dflag settings if present. Changed syslum0 output kw to modeloffvalue, indicating the off-eclipse value of the data.
;9oct15 added optional linear and quadratic offset to output (lumfin)
;15oct added inlcdigi param. Replaced int function calls by floor. Moved the adding of offsets, slopes (from 9oct15) to be before the combination of model and input data (for tflag=2,3). Removed adding of ooff from modeloffvalue output param
;17oct15 modified geolcalc to use array ops; is much faster
;19oct15 also geotcalc with array ops. Fixed bug from 17oct
;25oct15: inserted hook that avoids NAN output if lumbhafter = 0 in star::gettrlump (as can happen in grazing transits)
;4Nov15:  added self.scale to all body types; removed dependency on the fixed extensions within arrays (e.g the radfrac=0.45).  
;6Nov15: removed all assumptions about symetrical occulters/emmitters or square-arrays from  star::gettrlump. Sizes of one array (t or lform) will now be constant in the other's pix-coord system, avoiding position-dependent rounding effects. Changed self.lum and .lconv from float to double
;8nov15: introduced setting of self.scale during init process which gives an exact surface area in pixel-units
;9nov15 in gettrlump, removed the case of regridding tformoc array to lformbh, since the congrid introduces uncorrectable errors.
;10nov15 increment size of lformbh to next larger odd one; then adding offset to objects' baryctr, for symmetry. Changed radfac to 0.48.
;12nov15 fix in gettrlump to treat correctly point objects.
;27nov15 added input for u+,u- and q1,q2  coefficients for the quad lim-dark law, and supression of non-physical LD for linear and quad laws. 
;1dic2015 implemented  rapid calculation (fastflag)  in gettrlump for simple 2-body occultations without display.
;24jul16 inserted totecflag for total eclipses in gettrlump to avoid N/A's or very small values in occulted flux.
;25jul16 fixed another condition that may lead to N/A fluxes in gettrlump, in cases of total occultations.
;11apr19  changed display flags to -1,0,1 for a more predictable behaviour. 
;12apr19  added sqrtecsin, sqrteccos as alternative params for ecc, omega
;13apr19 fixed seed to onoise that may give noisy lcs that are identical
;26-29apr19 added mechanism to raise oobflag and jump the main loop if out-of-bound conditions are encountered
;9may19 added gsize multipliers. Fixed bug in flux-calc that became apparent on total occultations behind large occulters with fastflag=1, in line rellosproj=... 
;9may19 fixed error of returing different sized tarr, lumarr if resampling was used and oobflag got raised 
;16may19  only some improvements in commenting in gettrlump
;17may19 test version for finestepping in time
;4jun19  modified radfac to leave 1-pix wide empty frame in object arrays tform, lform
;27jun19  fixed bug causing crashes in eclipses in up-down direction due to some arrays becoming 1D
;         finestepping is implemented. Works well if pixels in tform are larger than pixels in lform. Otherwise, some precision issues remain; see comments in code
;3jul19: raised temporalily the limit for tlscale to 1.5 for the sub-pix shifting of tformoc. Outcommented again (search for '3jul19' )
;15mar20 removed oobflag from command-line, is now handled through pnarr,pvarr
;22mar20  added small offset in calculation of intrarr, to avoid spurious in-transit determinations
;6jul20    made final elongation-plot (xpflag) intereactive.
;24jul20  modified readlc for adaptation to latest ufit and renamed to rdnumtab. Moved all generic anciliary routines into utm_ancil.pro 
;1dec2020  added getrad function for extobj. First functional version for extobj.
;11feb21  first test version with transit depth calculated from exofast_occultquad_cel.pro in gettrluma
;15feb21 renamed gettranslum to gettrlump  (p for pixel, a for analytical)
;17feb21 added anlycalcflag and tests if conditions for analytical calc are met
;23feb21 perform occultation calc (zind array) now only for occulting objects. Removed
;         consideration of point-objects from gettrlump, ~luma
;24feb21 added occultcelflag (switch between exofast_occult.pro and exofast_occult_cel.pro)
;18mar21 fixed bug that caused crash in test if anlycalc is possible, whenever oobflag got raised 

 
;COPYRIGHT (C) 1999,2021 Hans J. Deeg;   
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    any later version.
;
;    As a special exception, this program may be compiled  with the
;    Interactive Data Language (IDL) and be linked to libaries
;    pertaining to IDL.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program in file 'licence.txt'; if not, write to
;    the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
;    Boston, MA 02110-1301 USA
;----------------------------------
;Known Bugs /Suggestions for improvement:
;    -modify findinfront to not use radii but distance from baryctrs to edges of tform/lform (at least for extobj)
;   - add oflag =3 for normalized luminosity (=dF/F (=oflag1) +1)
;    -remove limit of 16 objects and requirement of consecutive numbering
;    -correction for light-travel time: we could take the z -distance (form observer) and advance / delay a body by the light-travel time from the system-center.
; orbit plots with objects in real size (evtl. overlapping the lform -images in correct order)

;     -check on bug in inlcdigit. Does it remove the digits from inlcname?
;    -insert  3 param law from Sing 2010. See Kipping 2015. Also, his Eq. 73-75 to convert from alpha params -> c params and the criterion test of eq. 28. See 
;      https://github.com/davidkipping/LDC3 for corresponding code.
;    - permit sq law parametrizations q1=(u1+u2)^2 ; q2=0.5u1(u1+u2)^-1 from Kipping13
; - insert 'power-2' LD-law from Morello+2017 (easy)
; - evtl also Claret 4 coeffc law (see also Morello+17)
;   - interpolation in log g, Teff and Fe/H space to Claret's 4 LD parameters, with EXOFAST's QUADLD; see Eastman+ 2013
;   -use residlevel for display-offset of final output (for tflag=2,3); as in ufit

; - permit magnitudes as input for Nlum paramter (automatic if lunit = mag) (? is this a good idea?)
;
;  - distrib of UTM as .sav package?
; -modfy to use float-values for tformoc and lform arrays?
;    -interface for externally defined arbitrary occulters (e.g. tformoc array; then define in setup as external. External input either as passed as cmd-line variable or by reading from file. Define sizes as grid-size. Define center of object freely.) 
;
; CLEANUPS TO DO:
;incrempos methods should be based on actual time and not on tinc (makes self.etime obsolete)
;check CREATE_STRUCTURE to flexibilize on-fly struct. creation for extocc


@utm_ancil    ;compile anciliary routines
@setupfile ;compile setupfile

pro plotorb,posarr,tcount,maxdist,sunit, orbtrflag
;plot windows of orbital projections
  maxpl=maxdist*1.0             ;maxplot (size of plot-range)
;print,maxpl
  !x.range= [-maxpl,maxpl]
  !y.range= [-maxpl,maxpl]
  !x.style=2                    ;extend axis ranges
  !y.style=2
  if not(orbtrflag) then begin  ;plot current positions only
     
     wset,0                     ;plot x-y (sky plane)
     plot,posarr[*,0,tcount],posarr[*,1,tcount],psym=4,xtitle='X ('+sunit+')',ytitle='Y ('+sunit+')'


     wset,2                     ;plot z-y plane (side view)
     plot,posarr[*,2,tcount],posarr[*,1,tcount],psym=4,xtitle='Z ('+sunit+')',ytitle='Y ('+sunit+')'
     oplot,[!x.crange[1]],[0.0],psym=2,symsize=3 ;symbol for observer's direction
     
     wset,4                     ;plot x-z (top view) plane
     !y.range= [maxpl,-maxpl]   ;to put observer at bottom  
     plot,posarr[*,0,tcount],posarr[*,2,tcount],psym=4,xtitle='X ('+sunit+')',ytitle='Z ('+sunit+')'
     oplot,[0.0],[!y.crange[0]],psym=2,symsize=3 ;symbol for observer's direction
     

  endif else begin                                    ;plot traces
     tmp=size(posarr,/dimensions) & n_object = tmp[0] ;number of objects
                                ;plot x-y plane (sky plane)
     wset,0 
     for i=0,n_object-1 do begin
        if i eq 0 then plot,posarr[i,0,0:tcount],posarr[i,1,0:tcount],psym=3,xtitle='X ('+sunit+')',ytitle='Y ('+sunit+')'
        if i gt 0 then oplot,posarr[i,0,0:tcount],posarr[i,1,0:tcount],psym=3
     endfor

     wset,2                     ;plot z-y plane (side view)
     for i=0,n_object-1 do begin

        if i eq 0 then plot,posarr[i,2,0:tcount],posarr[i,1,0:tcount],psym=3,xtitle='Z ('+sunit+')',ytitle='Y ('+sunit+')'
        if i gt 0 then oplot,posarr[i,2,0:tcount],posarr[i,1,0:tcount],psym=3
        oplot,[!x.crange[1]],[0.0],psym=2,symsize=3 ;symbol for observer's direction
     endfor
     wset,4                     ;plot x-z plane (top view)
     !y.range= [maxpl,-maxpl]   ;to put observer (at +z) at bottom  
     for i=0,n_object-1 do begin

        if i eq 0 then plot,posarr[i,0,0:tcount],posarr[i,2,0:tcount],psym=3,xtitle='X ('+sunit+')',ytitle='Z ('+sunit+')'
        if i gt 0 then oplot,posarr[i,0,0:tcount],posarr[i,2,0:tcount],psym=3
        oplot,[0.0],[!y.crange[0]],psym=2,symsize=3 ;symbol for observer's direction
     endfor
     !y.range= [-maxpl,+maxpl]    
  endelse                       ;orbtrflag
end

pro plotlcutm,t,l,tunit,oflag,tflag,lunit,syslum0,lum0,fullc=fullc
;fullc is for interactive final plot of entire lc with different x-labelling 

;l is utm output (lumarr)
;m would be the model by utm
;if tflag <= 1 (lum0 = 0)  : l=m
;if tflag =2:  l=lum0+m   
;if tflag =3 : l=lum0-m    l is a residual. 
  case tflag of                 ;recovery of transit-model
     0:  m=l
     1:  m=l
     2:  begin
        if oflag ne 1 then m=l-lum0              ;recover model 
        if oflag eq 1 then m=(l+1)/(lum0+1) -1.  ;recover model in units dF/F
     end
     3:  begin
        if oflag ne 1 then m=lum0-l            ; ;recover model
        if oflag eq 1 then m=(lum0+1)/(l+1)-1. ;recover model in units dF/F
        if oflag eq 0 then l+=syslum0 ;add total syslum to move the residual to off-transit level of m (else, l around zero)
     end
  endcase  
  
  if tflag le 1 then begin
     miny=min(l)                ;minimum for y display
     maxy=max(l)                ;max
  endif
  if tflag ge 2 then begin
     miny=min(m)<min(lum0) <min(l) ;minimum for y display
     maxy=max(m)>max(lum0) >max(l) ;max
  endif

  if max(t) ge 100000 then t=(t mod 100000) ;to avoid display problems with long numbers

  if not(keyword_set(fullc)) then begin ; use non-interactive display for plots during UTM exec.
     case oflag of
        0: lumstr=lunit
        1: lumstr='dF/Fo'
        2: lumstr='delta-mag'
     endcase
     wset,5

     !x.range= [min(t),max(t)]
     if oflag ne 2 then !y.range=[miny,maxy] else !y.range=[maxy,miny] ;invert for mag-display
     !x.tickformat='(g10.6)'
     !x.ticks=3
     plot,t,m,xtitle='time ('+tunit+')',ytitle=lumstr,ystyle=2,xstyle=1 ;,xtickformat='(g10.6)',xticks=3;model as solid line
     if tflag ge 2 then oplot,t,lum0,psym=1                             ;input data as cross
     if tflag ge 2 then oplot,t,l,psym=4                                ;residuals (utm output) as diamond 
     !x.tickformat=''
     !x.ticks=0
     
  endif else begin              ; interactive final display 
     case oflag of
        0: lumstr=lunit
        1: lumstr='dF/F$_0$'
        2: lumstr='$\Delta$-mag'
     endcase
     if oflag eq 2 then yrange=[maxy+(maxy-miny)*0.1,miny-(maxy-miny)*0.1]                          ;invert for mag-display
     p=plot(t,m,'-g',thick=2,xtitle='time ('+tunit+')',ytitle=lumstr,ystyle=2,xstyle=1,yrange=yrange)       ;model as thick green line
     if tflag ge 2 then p=plot(t,l,'D',sym_color='CRIMSON',sym_thick=1,sym_size=0.7,/over) ;residuals (utm output) as red diamonds
     if tflag ge 2 then p=plot(t,lum0,'+',/over)                                                    ;input data as crosses
  endelse
end

;end

;###DEFINE INIT METHODS
function point::init,tinit,num,dflag ;inititalize values for all point-type objects; e.g all orbital paramters
  ival=1                             ;initialization value to be returned, 1 if init is sucessfull
  pi2 = 2*3.1415926535897932384626433832795028841971D 
  snum= str(num)
  self.num=num
  self.type=parasign(snum+"type")
  self.period=double(parasign(snum+"period",0.,/nowarn))
  self.etime=tinit              ; internal time counter
  eccinput=  float(parasign(snum+"ecc",0,/nowarn))   
  if eccinput lt 0.0 then begin 
     if dflag ge 2 then print,"WARN: negative eccentricity of ",eccinput," of object ",snum,". UTM ends with oobflag=1"
     ival =0                    ;optional OOB condition for negative ecc.
  endif
  self.ecc = abs(eccinput)      ; keep ecc always positve

  if self.ecc eq 0.0 then begin
     self.omega = float(parasign(snum+"omega",270,/nowarn))/!radeg ;arg of periapsis  (nowarn if circular)
  endif else begin
     self.omega = float(parasign(snum+"omega",270,warn='object '+snum+' is eccentric but no omega given; defaulted to 270deg' ))/!radeg ;arg of periapsi
  endelse
  
;if both parameters sqrtecsin = sqrt(ecc)*sin(omega) and sqrteccos= 
; sqrt(ecc)*cos(omega) are given, they override the eccentricity and omega
  eccosinyes = paratest(snum+"sqrtecsin") && paratest(snum+"sqrteccos")
  if eccosinyes then begin
     sqrtecsin = float(parasign(snum+"sqrtecsin")) ;sqrt(e) * sin(omega)
     sqrteccos = float(parasign(snum+"sqrteccos"))
     tmp=e_omega(sqrtecsin,sqrteccos)
     self.ecc=tmp[0]
     self.omega=tmp[1]
  endif                         ;if eccosinyes 

  self.omegadot = float(parasign(snum+"omegadot",0,/nowarn))/!radeg ;derivative of omega
  


; code that a given epoch overrides a given phase. if nothing given, set to 0
  trepochyes=paratest(snum+"trepoch") ;transit epoch: object at nue=3/2pi - omega 
  ecepochyes=paratest(snum+"ecepoch") ;eclipse epoch: object at nue=1/2pi - omeg
  phaseyes=paratest(snum+"phase")     ;phase  or Mean_anomaly / 2pi
  epochyes=paratest(snum+"epoch")     ;epoch of periapsis
  case 1 of                           ;assign epoch to object
     trepochyes : begin               ;transit epoch
        trepoch=double(parasign(snum+"trepoch"))
        te_tp=self.period*mean_anomal(pi2*3/4-self.omega,self.ecc)/pi2 ; timediff to periapsis
        if te_tp gt 0 then te_tp -= self.period                        ; keep to negative value (to get positive offset)
        self.epoch = trepoch-te_tp
        if ecepochyes or epochyes or phaseyes then if dflag ge 4 then print,"WARN: trepoch has overriden any of ecepoch, epoch, phase in object ",snum
        if dflag ge 3 then print,format='(a,i1,a,f15.5)',"Object ",snum,": Epoch for periapsis set to ",self.epoch
     end

     ecepochyes : begin         ;eclipse epoch
        ecepoch=double(parasign(snum+"ecepoch"))
        te_tp=self.period*mean_anomal(pi2/4-self.omega,self.ecc)/pi2 ; timediff to periapsis
        if te_tp gt 0 then te_tp -= self.period                      ;  keep to negative value (to get positive offset)
        self.epoch = ecepoch-te_tp
        if epochyes or phaseyes then if dflag ge 4 then print,"WARN: ecepoch has overriden any of epoch or phase in object ",snum
        if dflag ge 3 then print,format='(a,i1,a,f15.5)',"Object ",snum,": Epoch for periapsis set to ",self.epoch
     end
     phaseyes : begin           ;phase or Mean_anomaly / 2pi
        phase=double(parasign(snum+"phase"))
        self.epoch=tinit-(self.period*phase)
        if epochyes then if dflag ge 4 then print,"WARN: phase has overriden epoch in object ",snum
        if dflag ge 4 then print,format='(a,i1,a,f15.5)',"Object ",snum,": Epoch for periapsis set to ",self.epoch
     end
     epochyes :  self.epoch=double(parasign(snum+"epoch")) ;epoch of periapsis
     else : begin                                          ;none of them given: epoch=0
        self.epoch=0.                                      ;OK if period=0.0; warn only for moving objects
        if self.period ne 0.0 then if dflag ge 4 then print,"WARN: no epoch indicator given for object ",snum,". Epoch set to zero" 
     endelse 
  endcase
  self.inclin=double(parasign(snum+"inclin",90,/nowarn))
  self.rdist=double(parasign(snum+"rdist",0,/nowarn))
  self.pa=double(parasign(snum+"pa",270,/nowarn)) ;orbital position angle;dir: N towards E
  mainobj = parasign(snum+"mainobj","x",/nowarn)  ;main object. If absent or set to x or 99, no main object
  if STRMATCH(mainobj, 'x*', /FOLD_CASE) EQ 1 then self.mainobj=99 else self.mainobj=fix(mainobj) 
                                ; self -> geotcalc
  return, ival
end


function body::init,tinit,num,dflag ;inititalize body values
  snum= str(num)
;print,snum,self.mass,self.radi
  ival = self -> point::init(tinit,num,dflag) ;call init point
  if ival then begin
     self.radi=float(parasign(snum+"radi")) 
     self.oblate=float(parasign(snum+"oblate",0.,/nowarn))       ;oblateness = (R_eqat - R_pol) / R_eqat
     self.spinpa =  double(parasign(snum+"spinpa",0,/nowarn))    ;position angle of spin-axis dir: N towards E
     if self.radi gt 0.0 then begin
        self -> geotcalc
     endif else begin  
        if dflag ge 2 then print,"WARN: zero or negative radius of ",self.radi," of body-class object ",snum,". UTM ends with oobflag=1"
        ival=0
     endelse
  endif                         ;ival from point::init
  return, ival
end

function extobj::init,tinit,num,dflag  ;init external object
  snum= str(num)
  ival = self -> point::init(tinit,num,dflag) ;call init point
  self -> geotcalc,snum
  return, ival
end

function star::init,tinit,num,dflag
  snum= str(num)
  self.lum = float(parasign(snum+"lum",1))
  limblawstr=strlowcase(parasign("limblaw","lin")) ;limb-darken,default 'lin'
  case limblawstr of
     "lin": self.limblaw=1             ;linear limbd law
     "quad" : self.limblaw=2           ;quadratic laws
     "pmquad" : self.limblaw=2         ;quadratic laws
     "usquad" : self.limblaw=2         ;quadratic laws
     "root": self.limblaw=3            ;root-square law
     else: begin
        print,'unknown limbdarkening law is specified in parameter limblaw'
        stop
     end
  endcase

;assignement of input coefficients to u1, u2
  ldin1 = float(parasign(snum+"limbd",0.3))
  if self.limblaw gt 1 then ldin2 = float(parasign(snum+"limbd2",0.3)) else ldin2=0. ;second coef
  case limblawstr of
     "pmquad" : begin                     ;u+, u- is given
        self.limbd[0] = (ldin1+ldin2) /2. ;u1 
        self.limbd[1] = (ldin1-ldin2) /2. ; u2 
     end
     "usquad" : begin           ;convert q1, q2  to u1,u2. Kipping2013,  eq. 15,16
        self.limbd[0] =  (2* sqrt(ldin1) * ldin2)
        self.limbd[1]  = sqrt(ldin1)*(1-2*ldin2)
     end
     else : begin               ;other cases that use the input coeffcs directly 
        self.limbd[0] = ldin1   ;u1 
        self.limbd[1] = ldin2   ;u2
     end
  endcase
  
                                ;check if LD coeffc are physically plausible; for lin and quad law
  badldflag = 0 
  if self.limblaw eq 1 then if self.limbd[0] lt 0 or self.limbd[0] gt 1 then badldflag = 1
  if self.limblaw eq 2 then begin ;  criteria for u1, u2 from Kipping2013, eq. 8, inverted for out-of-boundary
     if (self.limbd[0] + self.limbd[1]) gt 1 then badldflag =1
     if self.limbd[0] lt 0 then badldflag =1
     if  (self.limbd[0] + 2* self.limbd[1]) lt 0 then badldflag =1
  endif
  if not(badldflag) then begin                  ;good ld's found
     ival = self -> body::init(tinit,num,dflag) ;call init 
     if ival then begin 
        self ->geolcalc                      ;calc luminosity distrib
        self.lformoc=float(self.lform)       ;copy init lum to occulted lum
     endif 
  endif else begin              ; bad LDs
     if dflag ge 2 then print,"WARN: Unphysical LD params of stellar object ",snum,". UTM ends with oobflag=1"
     ;print,'LDs are ', self.limbd[0],  self.limbd[1]
     ival=0                  
  endelse
  return, ival
end

function ring::init,tinit,num,dflag
  snum= str(num)
  self.transp = float(parasign(snum+"transp",0.5))
  ival = self -> body::init(tinit,num,dflag) ;call init body
return, ival
end


;####DEFINE PRINT METHODS 
pro body::print   ;method to print body obj	
 	print,'# general values of object ',str(self.num)
	print,'num',self.num
	print,'type ',self.type
	print,'radi',self.radi
;	print,'mass',self.mass
	print,'period',self.period
;	print,'phase',self.phase
	print,'inclination',self.inclin
        print,'pos angle',self.pa
	print,'rdist',self.rdist
	print,'positions',self.pos
	print,'mainobj',self.mainobj
end

pro star::print   ;method to print subtype star
	self -> body::print
 	print,'# star specific'
	print,'lum',self.lum
	print,'limbd',self.limbd
end

pro ring::print   ;method to print subtype ring
	self -> body::print
 	print,'# ring specific'
	print,'transparency',self.transp
end


pro body::helpobj  ;prints help about this object and stops
	help,/obj,self
	stop
end

pro point::printpos
	print,format='("nr: ",i3," ",3f14.4," ",a)',self.num,self.pos,self.type
end

pro body::displformocinit  ;initializes display of lform array 
;        gsize=size(self.lformoc,/dimension) & gs=gsize[0] 
	window,self.num+10,xsize=(200),ysize=(200),xpos=0,ypos=(self.num*240)+30,title=str(self.num)+' '+self.type+': luminosity'
end

pro body::disptforminit  ;initializes display of tform array 
;        gsize=size(self.lform,/dimension) & gs=gsize[0] 
	window,self.num+20,xsize=(200),ysize=(200),title=str(self.num)+' '+self.type+': transparency',xpos=220,ypos=(self.num*220)
end

pro star::displformoc  ;displays lformoc array 
	wset,self.num+10
	tvscl,0>(250<congrid(self.lformoc,200,200))
end

pro body::disptformoc  ;displays tform array 
	wset,self.num+20
	tvscl,0>(250<congrid(self.tform,200,200))
end


pro disp2Dform, form, title,plotsiz, winnr, xpos,ypos        ;displays a 2D array 'form' withg correct aspect ratio
  plotsiz=200                   ;max size of sidelength of plot
  formsize=size(form)
  aspect=float(formsize[1])/formsize[2]
  if aspect ge 1 then begin 
     xplotsiz=plotsiz
     yplotsiz=fix(plotsiz/aspect)
  endif else begin
     yplotsiz=plotsiz
     xplotsiz=fix(plotsiz*aspect)
  endelse
  window,winnr,xsize=xplotsiz,ysize=yplotsiz,title=title,xpos=xpos,ypos=ypos
  wset,winnr
  tvscl,0>(250<congrid(form,xplotsiz,yplotsiz))
end



;###DEFINE FUNCTIONS RELATED TO POSITON
;function subbody::getmainobjnum ;returns number of mainobj
;	return,self.mainobj
;end

function point::getpos   ;returns position of object
	return,self.pos
end


pro ring::incrempos,tinc,obj     ;position of ring is mainobj pos!
	mainobj=self.mainobj
	mainpos= obj(mainobj) -> getpos()   ;get pos of mainobj
	self.pos=mainpos  	    ;update pos
end


pro point::incrempos,tinc,obj   ;sets/increments phase and position
  pi2 = 2*3.1415926535897932384626433832795028841971D 
  self.etime=self.etime+tinc    ;internal time counter since we need t - t0 
;t0=self.epoch  ;epoch
;t=self.etime  ;current time
  Mt=(self.etime-self.epoch)*pi2/self.period ;mean anomaly
  if self.rdist ne 0. then begin             ;need to calculate something...
     ecc= self.ecc                           ;eccentricty
     if ecc ne 0. then begin                 ;
        Et=ecc_anomal(Mt, ecc)               ;eccentric anomal
        nue=true_anomal(Et,ecc)              ;true anomal
        radi=self.rdist*(1-ecc*cos(Et))      ;radius 
     endif else begin                        ;shortcut for circ. orbit
        nue=Mt
        radi=self.rdist
     endelse
;print,'etime,Mt,Et,nue,radi:',self.etime,Mt,Et,nue,radi
     xe=-radi*sin(nue)          ;ellipse's natural coord sys
     ye=radi*cos(nue)
     ze=0.0
;print,'ellipse coords: xe,ye',xe,ye
;print,'omega',self.omega

;projection into a coordsys rotated by arg of periapsis (omega)
     omega=self.omega + self.omegadot* (self.etime-self.epoch) ;current arg of periapsis
     xw=Xe*cos(omega)-Ye*sin(omega)
     yw=Xe*sin(omega)+Ye*cos(omega)
     zw=Ze
;xtest=radi*cos(omega+nue)  ;gives value of yw w/o rotation matrix
;print,'ellipse rotated by omega: xw,yw',xw,yw
;print,'xtest',xtest
;projection into coordsys rotated by inclination
     tet=self.inclin*!dpi/180.
     xi=xw*cos(tet)    	
     yi=yw           
     zi=xw*sin(tet)          

;rotation matrix for position angle resp. longitude of the ascending node
     del=self.pa*!dpi/180.
     self.pos[0]=cos(del)*xi-sin(del)*yi           
     self.pos[1]=sin(del)*xi+cos(del)*yi           
     self.pos[2]=zi                             
                                ;print,'final xyz',self.pos
  endif else begin              ;rdist eq 0
     self.pos[0]=0.        
     self.pos[1]=0.           
     self.pos[2]=0.
  endelse

  mainobj=self.mainobj 
  if mainobj ne 99 then begin
     mainpos= obj(mainobj) -> getpos()      ;get pos of mainobj
     self.pos=self.pos+mainpos  	    ;add pos
  endif
  
end



;###DEFINE FUNCTION RELATED TO GEOM FORM
pro body::geotcalc
                                ;calculates bodies 2D tranparency distribution
 
  gsize=size(self.tform,/dimension) ;get size of form-array
 
  radfac=(float(gsize[0]-2)/gsize[0])/2. ;factor that describes ratio between radius of object (in pix) and array dimension.
                                ;is designed so that an empty 1-pix wide frame remains around object.

 gs3=gsize[0] * 3                  ;array-size 3 times as big
  gbig=bytarr(gs3,gs3)              ;make array 3 times as big
  grad=gs3*radfac                   ;diameter will be ~90% of size of array
  gctr=gs3/2                        ;x,y position-index of array center
  grad2 = grad^2                    ;square of radius

  xdi=lindgen(gs3,gs3) mod gs3 - gctr ;array of xdist from ctr
  ydi=transpose(xdi)                  ;same for y-dists
  r2=(xdi^2+(ydi/(1-self.oblate))^2)  ;array with oblateness-corrected radial distances-squares
  gbig[where(r2 le grad2)]=255        ;assign opacity
  if self.spinpa ne 0.0 then begin    ;rotate
     gbig=rot(gbig,-self.spinpa,missing=0.,/interp) ;ccw rotation
  endif 
                                ;the big array gbig is resized to gsize, 
                                ;and transparency of border pixels will be interpolated 
  self.tform=rebin(gbig,gsize[0],gsize[1]) ;rebin to orig size
; derive size-scale of array
  a_tform_pix = total(self.tform,/double)/255D                   ;obj-cross-sect in pix-units
  a_tform_teo = !dpi * double(self.radi)^2 *(1-self.oblate)      ;expected obj-cross-sect in phys-units
  lscale= sqrt( a_tform_teo /  a_tform_pix )  
  self.scale=[lscale,lscale]    ;lenght-scales in sunit/pix; same for x and y
; print,'self.scale',self.scale

  self.bctr=gsize/2.            ;x and y value of barycenter within array. id of ctr-pixel is fix(self.bctr)

end

pro ring::geotcalc
                                ;calculates rings 2D tranparency distribution
  gsize=size(self.tform,/dimension) ;get size of form-array
  radfac=(float(gsize[0]-2)/gsize[0])/2. ;factor that describes ratio between radius of object (in pix) and array dimension.
                                ;is designed so that an empty 1-pix wide frame remains around object.

  gs3=gsize[0] * 3                  ;array-size 3 times as big
  gbig=bytarr(gs3,gs3)              ;make array 3 times as big
  grad=gs3*radfac                   ;outer diameter will be ~90% of size of array
  girad=grad*self.rdist/self.radi
  gctr=gs3/2 
  grad2 = grad^2                ;square of outer radius
  girad2 = girad ^ 2            ;sq of inner radius
  tet=(90.0-self.inclin)*!dpi/180.
  sintet=sin(tet)
;transp=exp(alog(self.transp)/sintet) ;transparency in depend of inclin
  transp=self.transp            ;transparency indepent of inclin
  transpc = fix(255*(1-transp)) ;here: 255: opaque, 0: transp
  xdi=lindgen(gs3,gs3) mod gs3 - gctr ;array of xdist from ctr
  ydi=transpose(xdi)                  ;same for y-dists
  r2=(xdi^2+(ydi/sintet)^2)           ;array with oblateness-corrected radial distances-squares
  gbig[where(r2 le grad2 and r2 ge girad2)]=transpc ;assign opacity
  if self.spinpa ne 0.0 then begin                  ;rotate
     gbig=rot(gbig,-self.spinpa,missing=0.,/interp) ;ccw rotation
  endif 
                                ;the big array gbig is resized to gsize, and border pixels
                                ;will be interpolated 
  self.tform=rebin(gbig,gsize[0],gsize[1]) ;rebin to orig size
  self.scale=self.radi/(radfac*gsize)      ;scale in x and y of array
  self.bctr=gsize/2.   ;x and y value of barycenter within array. id of ctr-pixel is fix(self.bctr)
end

pro extobj::geotcalc,snum
  gsext=[100,100]                    ;x-y size of 2D array of external objects
  extarrname=parasign(snum+"extname") ;name of .sav file with 2D transp.array named 'transarr'
  restore,extarrname                  ;restores transarr
  
  sizex=float(parasign(snum+"sizex"))      ;size of object along X in sunit
  sizey=float(parasign(snum+"sizey",sizex)) ;size along Y; use X if not given
  self.scale= [sizex / gsext[0], sizey/gsext[1] ] ;x and y scales, in length-per-pix

  bctrx=float(parasign(snum+"bctrx",0.5))       ;baryctr of object along X relative to array. default 0.5 is middle
  bctry=float(parasign(snum+"bctry",0.5))    ; along Y;
  self.bctr=[bctrx*gsext[0],bctry*gsext[1]]            ;position of object's barycenter within array, in pix-units

;  self.radi=pixsiz*(gsextx/2. > gsexty/2.)  ;size is the larger of the two dimensions  (Maybe use the diagonal, for rotations?)
  self.spinpa =  double(parasign(snum+"objpa",0,/nowarn)) ;position angle (rotation in sky plane) dir: N towards E
;  if self.spinpa ne 0.0 then begin                  ;rotate
;     transarr=rot(temporary(transarr),-self.spinpa,missing=0.,/interp) ;ccw rotation
;  endif 
  self.tform=transarr 
end


pro star::geolcalc;,badldflag=badldflag
                                ;calculates bodies 2D luminosity distribution


  gsize=size(self.tform,/dimension) ;get size of form-array
  radfac=(float(gsize[0]-2)/gsize[0])/2. ;factor that describes ratio between radius of object (in pix) and array dimension.
                                ;is designed so that an empty 1-pix wide frame remains around object.

  gs3=gsize[0] * 3                  ;array-size 3 times as big
;  gbig=bytarr(gs3,gs3)            ;make array 3 times as big
  grad=gs3*radfac               ;radius will be ~90% of size of array
  gctr=gs3/2 
  grad2 = grad^2                ;square of radius
  cl=self.limbd                 ;limbdarkening coefs, array of 2 values

  xdi=lindgen(gs3,gs3) mod gs3 - gctr ;array of xdist from ctr
  ydi=transpose(xdi)                  ;same for y-dists
  r2=(xdi^2+(ydi/(1-self.oblate))^2)  ;array with oblateness-corrected radial distances-squares
  
;  if badldflag then begin       ;bad ld; generate cross
;     gbig = bytarr(gs3,gs3)
;     for i=0,gs3-1 do gbig[i,i] = 255 
;     for i=0,gs3-1 do gbig[i,- i] = 255 
;  endif else begin
     case self.limblaw of       ;below expressions cause arithm. error from negative roots, which byte( converts to 0 values.
        1: gbig=byte(255*(1-cl[0]*(1-sqrt(1-(r2/grad2)))))                             ;lin law bright distr
        2: gbig=byte(255*(1-cl[0]*(1-sqrt(1-(r2/grad2)))-cl[1]*(1-sqrt(1-(r2/grad2)))^2))  ;brightness distr with quad limb darkening
        3: gbig=byte(255*(1-cl[0]*(1-sqrt(1-(r2/grad2)))-cl[1]*(1-sqrt(sqrt(1-(r2/grad2)))))) ;brightness distr with root-square limb darkening
     endcase
;  endelse
  if self.spinpa ne 0.0 then begin                  ;rotate
     gbig=rot(gbig,-self.spinpa,missing=0.,/interp) ;ccw rotation
  endif 
  
                                ;the big array gbig is resized to gsize, and border pixels
                                ;will be interpolated 
  self.lform=rebin(gbig,gsize[0],gsize[1]) ;rebin to orig size
                                ;below conversion factor from pixels surface-brightness units(0-255)
                                ;to luminosity unit (like Lsol)
  self.lconv=self.lum/total(self.lform)
end

;### INFORMAL ROUTINES

function point::getrad   ;returns radius of object
	return,0.
end

function body::getrad   ;returns radius of object
	return,self.radi
end

function extobj::getrad         ;returns radius of object
;calculates distance from the baryctr to most distant corner of rectangle in physical units
  tmp=size(self.tform)  
  xysize=tmp[1:2]*self.scale    ;physical x,y size
  xybary=self.bctr*self.scale   ;physical xy posn of baryctr within tform
  xydist=xybary > (xysize-xybary) ;get the distance from barycenter to the more distant limit in x,y
;rad=sqrt(total(xydist^2)) ; distance from the baryctr to most distant corner of rectangle
  rad=max(xybary)  ;this is a fudge to avoid errors from on-overlaping array and only ok for barycntr in center of array
  return,rad
end

function point::getrdist    ;returns distance of object
                               	return,sqrt(self.pos[0]^2+self.pos[1]^2+self.pos[2]^2)
end


function point::gettform   ;returns tform array of object
	return,[[0.]]
end

function body::gettform   ;returns tform array of object
	return,self.tform
end

function extocc::gettform   ;returns tform array of object
	return,self.tform
end


function star::getlum
	return,self.lum
end

;tests if object is suitable (=1) for analytical calculation
function point::anlyoktest  
        return,1    ;point always suitable
end

function body::anlyoktest   
  if self.oblate eq 0. then anlyok =1 else anlyok=0 ;suitable obj has to be circular
        return,anlyok
end

function star::anlyoktest  
      anlyok=self ->  body::anlyoktest()
      if self.limblaw ge 3 then anlyok = 0  ;only lin or square limblaw
        return,anlyok
end

function ring::anlyoktest  
        return,0    ;ring not suitable
end

function extobj::anlyoktest   
        return,0    ;external object not suitable
end


;### RELATED TO OCCULTATION OF LIGHTPATH
function body::findinfront,i,zind,pos,rad
     ;return indices of objects in front of it
n_object=n_elements(zind)
frind=-1   ;index of objects that are in front, -1: none
is=zind[i]                      ;sorted i index
for k=i+1,n_object-1 do begin   ;from current object to closest in z
  ks=zind[k]                     ;sorted k index 	
  if abs(pos[is,0]-pos[ks,0]) lt rad[is] + rad[ks] then begin  ;test in x
    if abs(pos[is,1]-pos[ks,1]) lt rad[is] + rad[ks] then begin  ;test in y
;      print,'object ',ks,' is in front of obj. ',is
      if frind[0] eq -1 then frind=[ks] else frind = [frind[*],ks]
    endif
  endif
endfor       
return,frind                  
end

function star::getlum  ;return intrinsic luminosity
       return, self.lum
end

pro star::resetlformoc   ;resets lformoc (occulted luminosity array)
       self.lformoc=float(self.lform)
end



function star::gettrluma,objfr,occultcelflag=occultcelflag
  
                                ;calc luminosity with 1 object from  objfr in front using M&A model
  debug=0                       ;set to 1 for display of lots of debug-info
  debug2=0                      ; if set to 1, plots tformoc etc
  debug3=0                      ;if set, prints more details on positions of the various arrays m ainly for finestepping
  if self.limblaw gt 2 then message,'UTM gettrluma: Only linear or square limblaws supported'
  lrad=self.radi                ;radius of lum-obj in xy units
  prad= objfr[0].radi           ;radius of occulter in xy units
  delpos=self.pos - objfr[0].pos            ;delta position in xy coord units
  bxy=sqrt(delpos[0]^2 + delpos[1]^2)       ;impact param  in xy units
  brs= bxy/lrad                             ; b in units of star size
  prs = prad/lrad                           ; p in units of star size
  ld=self.limbd
  if occultcelflag then begin
     exofast_occultquad_cel,brs,self.limbd[0],self.limbd[1],prs,normflux
 ;    print,'using exofast_occultquad_cel' & stop
  endif else begin
     exofast_occultquad,brs,self.limbd[0],self.limbd[1],prs,normflux
 ;   print,'using exofast_occultquad' & stop
  endelse
  return,self.lum*normflux      ;return total luminosity in units of lunit
end


function star::gettrlump,objfr,dispstarflag=dispstarflag,finestep=finestep
                                ;calc luminosity (pixel-array based) with objects in vector objfr in front 
  debug=0                       ;set to 1 for display of lots of debug-info
  debug2=0                      ; if set to 1, plots tformoc etc
  debug3=0                      ;if set, prints more details on positions of the various arrays m ainly for finestepping
  lsize=size(self.lform,/dimension)                                               ; size of lform in lform pix-units
  lscale=self.scale &  if debug or debug3 then print,'lscale(sunit/pix): ',lscale ; size of lform pixels, in (sunit)/pix
  lsizexy= lsize*lscale                                                           ; sizs of lform in xy coords
  lminxy=self.pos[0:1] - self.bctr * lscale                                       ;lower xy-coords of lform

  tmaxxy=dblarr(2) & tminxy=dblarr(2) ;min and max x-y extent of transiting body
  tindmin=intarr(2)                   ;minimum index of part of tform that is occulting lform  
  tindmax=intarr(2)                   ;maximum index  "   "
  lindmin=intarr(2)                   ;minimum index of subarray of lform corresponding to tform 
  lindmax=intarr(2)                   ;maximum index  "   "

  self.lformoc=float(self.lform) ;set lformoc array to intrinsic luminos.  ;use double for TINY occulters?
  n_object=n_elements(objfr)
  if n_object ge 2 or dispstarflag then fastflag=0 else fastflag =1 ;fastflag only for single-object occultation without display  
                                ;fastflag=0    ;outcomment for testing
  if debug then print,'fastflag: ',fastflag
  for i=0,n_object-1 do begin
     lumorg = total(self.lformoc,/double) ;get luminosity of lformoc before transit  
     if debug then print,'.......loop.for.occulter.......'
     if debug then print,'lumorg: emitter pix-flux before occ (unscaled): ',lumorg
     if debug then print,'lumorg*self.lconv: pixel-flux before occ (*self.lconv): ',lumorg*self.lconv
     tpos=objfr[i].pos
     tform=objfr[i].tform 
     tsize=size(tform,/dimension)                                                                   ;get size of tform 
     tscale=objfr[i].scale   & if debug or debug3 then print,'tscale(sunit/pix): ',tscale           ; xy-unit (sunit) /pixel
;       ;calc indices of tform that are occulting
     tsizexy = tsize*tscale                ;size of tform in xy (sunit) 
     lsize_in_t = round(lsizexy/tscale) ;size of lform in tform pix-units (this way, the diff tindmax-tindmin will be the same for all positions)
                                ; tind, lind:   tform-index, lform-index
     lmin_in_t_exact= (lminxy - tpos[0:1])/tscale + objfr[i].bctr         ;lower left corner of lform in tform-pix-units
     if debug or debug3 then print,'lmin_in_t_exact',  lmin_in_t_exact    ;lower left indize of lform in tform-pix-units
     lmin_in_t = round(lmin_in_t_exact)
                                ;lformoff_in_t  = lmin_in_t_exact - lmin_in_t ;residual offset of lform in tform-pix-units
                                ;if debug or debug3 then print,'lformoff_in_t',  lformoff_in_t
     tindmin=([0,0]> lmin_in_t) < (tsize-1)                                          ;lower index of tform-array that is occulting lform
     tindmax=( (tsize-1) < (lmin_in_t + lsize_in_t -1)) > [0,0]                      ;upper index 
     tindmax = tindmax >tindmin                                    ;fixes rare crashes where tindmax < tindmin. Origin is to be investigated
     tformoc=tform[tindmin[0]:tindmax[0],tindmin[1]:tindmax[1]]    ;subarray of tform that is occulting lform
;       a (trailing) y-dimens of size 1 might have been removed in prior subscript op. Need to revert tformoc to 2D array
;       (not needed for x-dimens, as lformbh remains 2D even if its x-size is 1 
     if (tindmax[1]-tindmin[1]) eq 0 then tformoc=reform(tformoc,tindmax[0]-tindmin[0]+1,1,/overwrite)



     if debug or debug3 then   print,'tsize,  objfr[i].bctr',tsize,  objfr[i].bctr
     if debug or debug3 then   print,'tindmin,tindmax',tindmin,tindmax
     if debug or debug3 then   print,'size(tformoc)',size(tformoc,/dim)
     


;       ;calc of indices of lform that are being occulted
     tlscale = tscale/lscale                     ;ratio tscale to lscale 
     tsize_inl = round(tsize*tlscale) > [1,1] ;size of tform in lform pix-units (this way, the diff lindmax-lindmin will be the same for all positions)
     if  total(tsize_inl mod 2) eq 0 then begin ;both dims of lformbh array would become even; increment to next odd (thus avoiding asymetric configs)
        tsize_inl += 1
        bctroff =-0.5           ;adjustment to barycenter position of occulter in lform array
     endif else bctroff = 0.
     tmin_in_l_exact= (tpos[0:1] - tscale* objfr[i].bctr- lminxy)/lscale +bctroff     ;lower left corner of tform in lform-pix-units
     tmin_in_l = round(tmin_in_l_exact)                                               ;lower left indice of tform in lform-pix-units
                                ; tformoff_in_l=tmin_in_l_exact - tmin_in_l                                     ;residual offset of tform in lform-pix-units

     lindmin=[0,0]> tmin_in_l                       ;lower index of lform-array that is behind occulter tform
     lindmax=(lsize-1) < (tmin_in_l+tsize_inl-1)    ;upper index
     if debug or debug3 then print,'lsize,  self.bctr',lsize,  self.bctr
     if debug or debug3  then print,'lindmin,lindmax,',lindmin,lindmax
     lformbh=self.lformoc[lindmin[0]:lindmax[0],lindmin[1]:lindmax[1]] ;subarray of lform that is behind occulting object 

     if debug or debug3 then   print,'size(lformbh)',size(lformbh,/dim)
;       a (trailing) y-dimens of size 1 might have been removed in prior subscript op. Need to revert lformbh to 2D array
     if (lindmax[1]-lindmin[1]) eq 0 then lformbh=reform(lformbh,lindmax[0]-lindmin[0]+1,1,/overwrite)

     lformbh_npix=product(size(lformbh,/dim)) ;number of pix in lformbh     
;
;  if debug then print,'total(lformbh) pre-subsh', total(lformbh,/double)          ;luminosity of lformbh before transit
     lumbhorg=total(lformbh,/double) ;luminosity of lformbh before subshift and transit
     if debug or debug3 then print,'lumbhorg: ', lumbhorg
     
     if finestep then begin     ;perform sub-pixel shift on tformoc if tlscale ge 1 (WORKS WELL) or on lformbh if tlscale lt 1 (SOME ISSUES)
                                ;diff between lower-left corners of lformbh and tformoc, in tform-pix-units
        lbhmin_in_t =  lmin_in_t_exact+lindmin/tlscale       ;lower-left corner of lformbh                                
        diff_ll_in_t =  tindmin - lbhmin_in_t                ;where tindmin is lower-left-corner of tformbh 

                                ;diff between upper right corners of lformbh and tformoc, in tform-pix-units
        lbhmax_in_t = lmin_in_t_exact+(lindmax+1)/tlscale       ;upper-right corner of lformbh
        diff_ur_in_t =  (tindmax+1) - lbhmax_in_t               ;where tindmax+1 is upper-right corner of tformbh
                                ;               tformoff_in_t = diff_ur_in_t

                                ; select the weighting between diff_ll and diff_ur
        ur_weight=tpos[0:1] ge self.pos[0:1]
        ll_weight = 1-ur_weight
        
        tformoff_in_t=ur_weight * diff_ur_in_t + ll_weight*diff_ll_in_t
;tformoff_in_t=diff_ll_in_t


        if debug3 then print,'tlscale',  tlscale
        if debug3 then print,'lbhmin_in_t',  lbhmin_in_t
        if debug3 then print,'tindmin',  tindmin
        
        if debug3 then print,'lbhmax_in_t',  lbhmax_in_t
        if debug3 then print,'tindmax+1',  tindmax+1
        
        if debug3 then print,'ll_weight',ll_weight
        if debug3 then print,'ur_weight',ur_weight
        if debug3 then print,'tformoff_in_t',  tformoff_in_t
        if tlscale[0] ge 1 then begin ;pixels in tform are larger than pixels in lform; shift tformoc
;           if tlscale[0] ge 1.5 then begin ;pixels in tform are larger than pixels in lform; shift tformoc  ;3jul19: raised temp. to 1.5 for EBfit test
           
           
           tformoc_dim=size(tformoc,/dim) 
           if min(tformoc_dim) gt 1 then begin ;check to avoid shifts on 1D arrays
              if max(tformoff_in_t)   gt 1. then message,'shift-value  tformoff_in_t larger than 1' 
              tformoc_tmp=subshift(tformoc,tformoff_in_t[0],tformoff_in_t[1]) ; output is float
              tformoc=byte(tformoc_tmp)
           endif
        endif else begin        ;pixels in lform are larger than pixels in tform; shift lformbh
           lformoff_in_l = tformoff_in_t * tlscale
           if debug3 then print,'lformoff_in_l',  lformoff_in_l

;alternative calc directly in lform-pix-units
;diff between lower-left corners of lformbh and tformoc, in lform-pix-units
           tocmin_in_l= tmin_in_l_exact+tindmin*tlscale
           diff_ll_in_l = tocmin_in_l - lindmin
;diff between upper-right corners of lformbh and tformoc, in lform-pix-units
           tocmax_in_l=tmin_in_l_exact+(tindmax+1.)*tlscale
           diff_ur_in_l = tocmax_in_l - (lindmax+1.)
           lformoff_in_l2=ur_weight * diff_ur_in_l + ll_weight*diff_ll_in_l
;CHECK WHY l2 and L do NOT CONICIDE
;effectof lumbhorg before or after subshift. the FUDGE factor. 
;CHECK ISSUES OF lformbhproj CONGRID (BUT CARE; MIGHT AFFECT TLCSALE le 1

           if debug3 then print,'lformoff_in_l2',  lformoff_in_l2

           lformbh_dim=size(lformbh,/n_dim) 
           if lformbh_dim eq 2 then begin                                                  ;check to avoid shifts on 1D arrays
              lformbh_tmp=subshift(lformbh,-lformoff_in_l[0],-lformoff_in_l[1],/padbil)    ; output is float
              lformbh=byte(lformbh_tmp)
           endif
           lumbh_shiftflux=(total(lformbh)-lumbhorg)*0.33 ;flux-variation of lumbh due to subshift. The 0.33 is a FUDGE FACTOR
           if debug3 then print,'lumbh_shiftflux',lumbh_shiftflux
           lumorg-= lumbh_shiftflux ;above flux-variation is subtracted from the unocculted flux  (only used for fastflag)
        endelse
     endif                                               ;finestep 
     if ~isa(lumbh_shiftflux) then lumbh_shiftflux=0.    ; define this offset here if not done in finesteping
     if debug2 then disp2Dform,lformbh,'lbh',220,27,950,50 

                                ; removed on 9nov15 the reverse gridding tformoc to lformbh, which was above here  
     
     
     if debug then print,'gridding lformbh to size of tformoc'                            
     lformbhproj=congrid(lformbh,tindmax[0]-tindmin[0]+1,tindmax[1]-tindmin[1]+1,/interp,/center) ;lformbh resample to size of tformoc     
     if debug or debug3 then   print,'size(lformbhproj)',size(lformbhproj,/dim)                          
     lumbhprojorg=total(lformbhproj,/double) ;original flux of occulted area
     if debug then print,'lumbhprojorg',lumbhprojorg
     if lumbhprojorg gt 0.0 then begin ;bypass a case (lumbhprojorg eq 0.0) that can happen at begin/end of transit 

        lformbhproj_npix=product(size(lformbhproj,/dim))           ;number of pix in lformbhproj
        proj_org_pixrat= float(lformbhproj_npix)/(lformbh_npix)    ;ratio Npixels resampled vs orig
        proj_org_fluxrat= lumbhprojorg/lumbhorg                    ;flux-ratio resampled vs orig
        proj_org_flux_vs_pix_rat = proj_org_fluxrat/proj_org_pixrat ;  flux/pix ratio resampled versus orig 
        if debug then print,'lfrombhproj/lfrombhp pix-rat, flux-rat: ',proj_org_pixrat,proj_org_fluxrat
        if debug then print,'flux-rat/pix-rat:',proj_org_flux_vs_pix_rat
        lformbhproj *= proj_org_flux_vs_pix_rat ;correct lformbhproj to flux expected from perfect regridding.
                                ;above multiplication always removes trailing 1-pixel length array dimensions. Get it back in next line
        if (tindmax[1]-tindmin[1]) eq 0 then lformbhproj=reform(lformbhproj,tindmax[0]-tindmin[0]+1,1,/overwrite)

        lumbhprojorgcor= lumbhprojorg *  proj_org_flux_vs_pix_rat ;correct also the pre-occultation total flux
        if debug then print,'lumbhprojorgcor',lumbhprojorgcor

        if debug2 then disp2Dform, tformoc, 'tformoc',200, 30, 220, 50
        if debug2 then disp2Dform,lformbhproj,'lbproj_pre',220,31,450,50 
        lformbhproj= temporary(lformbhproj)*(255-tformoc)/255. ;APPLY TRANSPARENCY OF OCCULTER TO LFORMBHPROJ 
        lumbhprojafter=total(lformbhproj,/double)             
        if debug or debug3 then print,'lhprojafter: ',lumbhprojafter        
        if debug2 then disp2Dform,lformbhproj,'lbproj_aft',220,28,700,50    
                                ;          rellosproj=(lumbhprojorgcor-lumbhprojafter)/lumbhprojorg ;relative flux-loss that occured in the projected array lformbhproj (which is correct value)
        rellosproj=(lumbhprojorgcor-lumbhprojafter)/lumbhprojorgcor ;relative flux-loss that occured in the projected array lformbhproj (which is correct value)   ;fixed version 9 may19
        
        if debug then print,'rellosproj',rellosproj
;           relbrite = lumbhprojafter/lumbhprojorg                   ;relative brightness after occultation in lformbhproj
;           ;calc for correction of DF/F from area-variation of occulter due to regridding between integer-sized arrays
        lformbhxy = lscale*(lindmax-lindmin+1)       ;xy-size of star area that is used in occultation-calc
        tformocxy = tscale*(tindmax-tindmin+1)       ;xy-size of occulter area (correct value)
        if debug then print,'lformbhxy',lformbhxy
        if debug then print,'tformocxy',tformocxy 
        lscalecorrfac = product(tformocxy)/product(lformbhxy) ;correction factor 
        if debug then print,'lscalecorrfac',lscalecorrfac     
        rellosdesired=rellosproj*lscalecorrfac ;that's the loss we want in lformbh
        if debug then print,'rellosdesired: ', rellosdesired
        if rellosdesired gt 1 then totecflag=1 else totecflag=0    ;flag for total (or near total eclipse)
        if rellosdesired gt 1 then rellosdesired=1D                ;set to compute flux loss


        desiredlumbh=lumbhorg*(1-rellosdesired)
        if debug then print,'desiredlumbh: ', desiredlumbh
        if not(fastflag) then begin ;full calc that projects the lformbhproj array back into self.lformoc

;          ;scale lformbhproj back to orig size of lformbh  and revert the flux-correction
           if total(size(lformbh,/dim)) eq 2 then lformbhoc = mean(lformbhproj)/proj_org_flux_vs_pix_rat else $ ;if lformbh is 1-pix size, use mean          
              lformbhoc=congrid(lformbhproj,lindmax[0]-lindmin[0]+1,lindmax[1]-lindmin[1]+1,/interp,/center) /proj_org_flux_vs_pix_rat ;re-grid and re-correct flux 
;     lformbhoc=congrid(lformbhproj,lindmax[0]-lindmin[0]+1,lindmax[1]-lindmin[1]+1,/inter,/minus) /proj_org_flux_vs_pix_rat ;scale lformbhproj back to orig size of lformbh;  
           
           lumbhafter=total(lformbhoc,/double)
           if  lumbhafter eq 0 then begin ;fudge to recognize bad congrid from ignored non-zero pixels at borders of lformbhproj
              lformbhoc=congrid(lformbhproj,lindmax[0]-lindmin[0]+1,lindmax[1]-lindmin[1]+1,/interp,/minus_one) /proj_org_flux_vs_pix_rat ;interpolating now
              lumbhafter=total(lformbhoc,/double)
           endif
           if lumbhafter eq 0 then totecflag=1 ;so it was really a total eclipse 

;          ;correct lformbhoc to have correct flux, which is desiredlumbh
;          ;apply correction algos pending on area ratio occulter/emitter              
           tvsl_area=product(tformocxy) / product(lsizexy)             & if debug then print,'area ratio tform/lform',tvsl_area
           if tvsl_area ge 0.002 then bigflag=1 else bigflag=0 ;use correction over full self.lformoc only for large occulters     
           if not(bigflag) then begin ;corr. over occulted zone; less sensitive to float-pt errors; visible squares in display
              if debug or debug3 then print,'SMALL OCCULTER' 
              if lumbhafter gt 0 then lformbhoc=lformbhoc*((desiredlumbh-lumbh_shiftflux)/lumbhafter)
           endif 
           self.lformoc[lindmin[0]:lindmax[0],lindmin[1]:lindmax[1]]=lformbhoc ;glue lformbhoc back into lformoc; line needed with/without bigflag code
           if debug2 then disp2Dform,lformbhoc,'lformbhoc',220,29,1150,50    
           if bigflag then begin ;corrects over full array; does not show squares but would have precision problems for small occulters.
              if debug  or debug3 then print,'BIG OCCULTER' 
              lumcorroff=desiredlumbh-lumbhafter                                ;lumin offset between desired and real one
              lumafteroccult=total(self.lformoc,/double)                        ;actual luminosity of lformoc after occultation
              lumcorrfac=double(lumafteroccult+lumcorroff-lumbh_shiftflux)/lumafteroccult ;correction factor to multiply into lformoc (close to 1 for small occulters,)
              if totecflag  then lumcorrfac =0.0
              if debug then print,'lumcorrfac: ',lumcorrfac

              self.lformoc=self.lformoc*lumcorrfac ;errors occur here since self.lformoc is float that is multiplied by value very close to 1           
              if debug then print,'emitter flux after occ (unscaled): ',total(self.lformoc)
           endif                ;bigflag
           if debug2 then self-> displformoc
        endif else lumoc_fast=lumorg-lumbhorg+desiredlumbh ;else: fastflag eq 1
      ;if totecflag  then stop
     endif else if fastflag then lumoc_fast = lumorg ;if lumbhprojorg ge 0.0   
  endfor                        ;loop over objfr[i]
  if not(fastflag) then   lumoc=total(self.lformoc,/double) $    ; luminosity of 'self' with occultations
  else  lumoc=lumoc_fast                                         ;fastflag eq 1
  if debug then print,'final emitter flux after occ (unscaled): ',lumoc
  if FINITE(lumoc, /NAN) then print,'An N/A otuput flux happened. Needs debugging:-('
  if FINITE(lumoc, /NAN) then stop
  if debug then print,'final emitter flux after occ (*self.lconv): ',lumoc*self.lconv
;  if debug then print,'final dF: ',(lumoc- lumorg) / lumorg  

  if debug then print,'-----------end of loop occulted obj----------' 
  return,lumoc*self.lconv       ;return total luminosity in units of lunit
end



;================================================
;============ UTM MAIN PROGRAM  =================
;================================================

pro utm,setfname,dflag,tarr,lumarr,intrarr=intrarr,modeloffvalue=modeloffvalue ;,oobflag=oobflag
  common parrs,pvarr,pnarr
;parameters
;  setfname  name of setupfile
;  dflag     flag for level of display details
;(arrays tarr and lumarr will be taken from command-line if defined there. Else, corresponding files need definition in setupfile.)
;  tarr      input array with time-points  (only used if tflag ge 1)
;  lumarr    ouput array (if tflag=1); input/output array (if tflag ge 2) of luminosity values of same kind as set in oflag.
;keywords:
;  intrarr   output array that is 1 if any eclipsing occurs and 0 else.
;  modeloffvalue output scalar with off-eclispe model value
; oobflag    output flag that is 1 if an out-of bounds condition was encountered
;save status of graphics sys variables that are modified
  xsysorg=!x  &ysysorg=!y & psysorg=!p
;read setupfile
  readsetupfile,setfname
  oobflag =0                              ; out of bounds flag
  if total(size(dflag)) eq 0 then dflag=1 ;default display-flag

  
;set display flags either from dflag or from overrides if -1 or 1

  tmp= fix(parasign("plcflag",'0',/nowarn)) ; plot entire lc at program end in interactive window
  if tmp eq 0 then  if dflag ge 1 then plcflag=1 else plcflag =0  
  if tmp ne 0 then plcflag = (1+tmp)/2 ;map flag values -1,1 to 0,1

  tmp= fix(parasign("plcsimflag",'0',/nowarn)) ;plot last 100 pts of lightcurve during simulation.
  if tmp eq 0 then  if dflag ge 2 then plcsimflag=1 else plcsimflag =0  
  if tmp ne 0 then plcsimflag = (1+tmp)/2 ;map flag values -1,1 to 0,1

  tmp= fix(parasign("plotorbflag",'0',/nowarn)) ;plot orbital positions
  if tmp eq 0 then  if dflag ge 3 then plotorbflag=1 else plotorbflag =0  
  if tmp ne 0 then plotorbflag = (1+tmp)/2 ;map flag values -1,1 to 0,1

  tmp= fix(parasign("dispstarflag",'0',/nowarn)) ;show figures of stars
  if tmp eq 0 then  if dflag ge 3 then dispstarflag=1 else dispstarflag =0  
  if tmp ne 0 then dispstarflag = (1+tmp)/2 ;map flag values -1,1 to 0,1

; further display flags without predefinition in dflag
  xpflag = fix(parasign("xpflag",0,/nowarn))       ;plot elongations (Xpositons)
  orbtrflag = fix(parasign("orbtrflag",0,/nowarn)) ;plot traces of orbital positions

  paradd,"oobflag","0"          ;reset oobflag if raised by previous invoke of utm 

;check if setupfile-preprocessing requested and execute it.
  preproc=parasign("preproc","none",/nowarn) ;name of function that will modify setup-params
  if preproc ne "none" then prepflag=1 else prepflag=0
  if prepflag then begin        ;do preprocessing of setupfile
     if dflag ge 4 then print,'Setupfile preprocessing with: ',preproc
     call_procedure, preproc                    ; CALL PREPROCESSING PROGRAM
     save_tmpset=fix(parasign("saveprepset",0)) ;flag to save preprocessed setup parameters
     if save_tmpset then begin
;create unique name for processed setup that is based on systime to 1/100s sec
        stime=frac(systime(1)/10000)*10000
        seed=stime               ;this, since randomu will modify seed
        ranu=randomu(seed)*10000 ;random number
        setmp=strcompress(string(format='(f15.2,f7.1)',stime,ranu),/remove_all)+'preprocd.set'
        savesetupfile,setmp
     endif
     prep_oobflag=parasign("oobflag",0,/nowarn) ;check if an oobflag was raised during preprocessing
     if prep_oobflag eq 1 then begin
        if dflag ge 2 then print,'WARN: Out of bound params while executing preprocessor ',preproc,'. UTM ends with oobflag=1' 
        oobflag=1
;        paraput,"oobflag","0"
     endif
  endif

;get units
  tunit=parasign("tunit","day")        ;get time-unit
  lunit=parasign("lunit","Lsun")       ;get luminosity-unit
  sunit=parasign("sunit","Rsun")       ;get size-unit
;  munit=parasign("munit","Msun")       ;get mass-unit


;initialize time-array and luminosity array
  tflag = fix(parasign("tflag",0)) 
  pflag=n_params()-2            ;For ufit: Nr params in cmd-line call to utm. Is 1 if tarr, 2 if tarr, lumarr.
  case 1 of
     tflag eq 0 :begin          ;make up times from tinit,tinc,tfin
;clock
        tinc=double(parasign('tinc',0))
        tinit=double(parasign('tinit',0.1))
        tfin=double(parasign('tfin',10))
        nlcpts=long((tfin-tinit)/tinc) +1 ;number points to simulate
        tarr=dblarr(nlcpts)
        tcount = 0L & t=tinit
        for tcount=0L,nlcpts-1 do begin
           tarr(tcount)=t
           t=t+tinc
        endfor 
     end
     
     (tflag ge 1)  : begin       ; input lc from file or cmd-line params
        if pflag eq 0 then begin ;input lc from file
; get input-lc, using one of several keywords
           if paratest('indatafile') then indatafile=parasign('indatafile')  else $
              if paratest('inlcname') then indatafile=parasign('inlcname') else $
                 if paratest('infilename') then indatafile=parasign('infilename') else $
                    print,'WARN: need to give input lightcurve' ;compatibility with old param-name'

           tmarr=rdnumtab(indatafile,tflag<2)
           tarr=tmarr[*,0]
           if paratest('inlcdigi') then begin
              ndigi=fix(parasign('inlcdigi')) ;number of precomma digits to keep
              tarr=digirem(tarr,ndigi)        ;keep only ndigi precomma/postcomma digits in time-column
           endif
           if tflag eq 2 or tflag eq 3 then lum0arr=tmarr[*,1] ;initial model vals
           nlcpts=n_elements(tarr)
        endif else begin            ;input lc in params tarr,lumarr
           nlcpts=n_elements(tarr)  ;tarr is already defined as param
           if ( tflag eq 2 or tflag eq 3 ) then lum0arr=lumarr ;initial  model vals
        endelse
     end
  endcase

  if total(size(lum0arr)) eq 0. then lum0arr=dblarr(nlcpts) ;if not defined yet, do so; needed in plotlcutm

;resampling stuff
  resample=fix(parasign('resample',1,/nowarn))
  if resample gt 1 then begin   ;resample time-points, lum0arr
     sample_span = float(parasign('sample_span',0.02043,warn='WARN: Resampling span set to 0.02043 (Kepler long cadence in days)'))
     nlc_resa= nlcpts*resample  ; number of resample pts
     tarr_resa=dblarr(nlc_resa)
     lum0_resa=dblarr(nlc_resa)
     for i=0,nlcpts-1 do begin  ;resampling loop
        for j=0,resample-1 do begin
           tarr_resa[j+i*resample]=tarr[i]+ sample_span*((j*2)-(resample-1))/(resample*2) ; equidistant pts within sample_span
           lum0_resa[j+i*resample]=lum0arr[i]                                             ; assign same data-val to all pts in same sample_span
        endfor
     endfor
;rename stuff 
     tarr_org=tarr
     lum0arr_org=lum0arr
     tarr=tarr_resa
     lum0arr=lum0_resa
     nlcpts = nlc_resa          ;use this larger value in remainder
  endif                         ;resample



  lumarr=fltarr(nlcpts)
  intrarr=intarr(nlcpts)        ;array that is 1 in-transit and 0 offtransit

  if oobflag then goto,endwithoob ;jump to end if oobflag raised


;define structures for the object-types

  gs=fix(parasign("gsize",36))  ;base-size of geometric arrays
  gs=gs-(gs mod 2)              ;decrement to next even nr
                                ;sizes of geometric arrays (tform,lform) to be used for various bodies 
                                ;they are forced to be odd numbers, to have a central pixel:
  gsizemultstar=float(parasign("gsizemultstar",4,/nowarn)) ;size-multiplier for stars
  gsizemultplan=float(parasign("gsizemultplan",1,/nowarn)) ;size-multiplier for planets
  gsizemultring=float(parasign("gsizemultring",2,/nowarn)) ;size-multiplier for rings

  gsstar=2*fix((gs*gsizemultstar)/2)+1                 ;star
  gsplan=2*fix((gs*gsizemultplan)/2)+1                 ;planet
  gsring=2*fix((gs*gsizemultring)/2)+1                 ;ring

  anlycalcflag= fix(parasign("anlycalcflag",'0',/nowarn)) ;if -1 always pixel-based calculation; 0 does analytical when possible, +1 only analytical and stops w message if not possible

  if anlycalcflag ge 1 then begin ;initalize only dummy 2x2 tform and lform arrays
     gsplan = 2
     gsstar = 2
  endif
  occultcelflag = fix(parasign("occultcelflag",'0',/nowarn))  ;switch between exofast_occult.pro and exofast_occult_cel.pro

  gsext=[100,100]               ;x-y size of 2D array of external objects
  p={point, num:0,type:'point',period:0D,epoch:0D,ecc:0D, omega:0D,omegadot:0D,etime:0D,inclin:0D,pa:0D,rdist:0D,pos:dblarr(3), mainobj:0} ; contains all orbiting objects
;  ep={eccpoint,ecc:0D, omega:0D,etime:0D, inherits point}
  a={body, radi:0.,oblate:0.,spinpa:0.,scale:dblarr(2),bctr:fltarr(2),inherits point} ; covers all round bodies with a defined radius. A body on its own does not cause eclipses (no tform)
  b={star,lum:0D,limblaw:0,limbd:fltarr(2),lform:bytarr(gsstar,gsstar),lformoc:fltarr(gsstar,gsstar),lconv:0D,tform:bytarr(gsstar,gsstar),inherits body}
  c={planet,tform:bytarr(gsplan,gsplan),inherits body}
;  d={subbody,mainobj:0,inherits body}
  e={moon,tform:bytarr(gsplan,gsplan),inherits body} ;this is an obsolete obj-class
  f={ring,transp:0.0,tform:bytarr(gsring,gsring),inherits body}
  g={extobj,tform:bytarr(gsext),spinpa:0.,scale:dblarr(2), bctr:fltarr(2),inherits point} ;external object
  h={extocc, inherits extobj}                                                             ;external occulter
;check CREATE_STRUCTURE to flexibilize on-fly struct. creation

;ii={extlum,lum:0.0,lform:bytarr(gsext,gsext),lformoc:fltarr(gsext),lconv:0.0,tform:bytarr(gsext),inherits extobj} ;ext luminous obj
;create the objects and assign initial values to them

;first,find number of objects by looking for largest param Ntype
;that is defined in setupfile
  n_object=0&while paratest(str(n_object)+"type") do n_object=n_object+1
  if dflag ge 3 then print,str(n_object)+' objects will be defined'
  posarr=fltarr(n_object,3, nlcpts) ;array that keeps x,y,z positions of objects

;objects of type(i) are created and initialized
  type=strarr(n_object)
  obj=objarr(n_object)
  for i=0,n_object-1 do begin
     type[i]=parasign(str(i)+"type")
     obj[i]=obj_new(type(i),tarr[0],i,dflag) ;call initialization method
;print,'obj_valid',obj_valid(obj[i])
     if obj_valid(obj[i]) then begin         ;valid object was created
        obj[i] -> incrempos,0.,obj           ;set initial position
        if obj_isa(obj[i],'star') then if dispstarflag then obj[i] -> displformocinit
     endif else begin 
        oobflag=1               ;raise flag if creation  failed due to OOB
;        print,'i=',i,': oobflag raised'
     endelse
  endfor                        ;i=0,n_object
  if oobflag then goto,endwithoob ;jump to end if oobflag raised

  occultobid=where(type ne "point") ;array with indizes of the occulters (all except point)
  n_occult=n_elements(occultobid)  ;number of occulters

;flag-vector for pixel-based / anlytical calculation
  anlyokflag= intarr(n_object)  ;bin-array, if 1 then object is suitable for analytical calc. 

  if anlycalcflag ge 0 and dispstarflag eq 0 then begin 
     for i=0,n_object-1 do anlyokflag[i] = obj[i] -> anlyoktest() ;test if objs compatible with analyt calc
  endif                                                           ; (else: anylokflag =[0,0,...] 

  if total(anlyokflag) lt n_object and anlycalcflag eq 1 then begin 
     message,'UTM: anlycalcflag = 1 but purely analytical calculation impossible (2D object display ON, or unsuitable objects)'  
  endif

  if dflag ge 1 then begin
     if anlycalcflag eq 1 then print,'UTM: using purely analytical occult. calc.'
     if anlycalcflag le 0 then begin  ;total(anlyokflag) has to indicate at least 2 objects suitably for anlycalc
        if total(anlyokflag) ge 2 then print,'UTM: using hybrid mode (analytical occult. calc. whenever possible)'
        if total(anlyokflag) le 2 then print,'UTM: using pixel-based occult. calculation'
     endif
  endif
;     print,'anlyokflag: ',anlyokflag

;  generate windows 
  if plotorbflag then window,0,xsize=400,ysize=320,title='sky (x-y) plane',xpos=0
  if plotorbflag then window,2,xsize=400,ysize=320,title='side-view (z-y) plane',xpos=400
  if plotorbflag then window,4,xsize=400,ysize=320,title='top-view (x-z) plane',xpos=800
  if plcsimflag  then window,5,xsize=800,ysize=320,title='light curve',xpos=400,ypos=200

;get some values of the objects that are needed in time-loop
  zind=intarr(n_object)         ;sort-index in z-direction
  pos=dblarr(n_object,3)        ;positions of all objects
  rad=dblarr(n_object)          ;radii of all objects
  for i=0,n_object-1 do rad[i] = obj[i] -> getrad()  
  maxdist=0D
;for i=0,n_object-1 do print, obj[i] -> getrdist()
;  for i=0,n_object-1 do maxdist = (maxdist > obj[i] -> getrdist())  

;initialize display
  for i=0,n_object-1  do if obj_isa(obj[i],'star') then if dispstarflag then obj[i]-> displformoc

;prepare header of output model-lc
  oflag = fix(parasign("oflag",0))
  case oflag of
     0: lumstr=lunit
     1: lumstr='dF/F'
     2: lumstr='delta-mag'
  endcase
  onoise=float(parasign("onoise",0,/nowarn))    ;opt noise to be added to output
  ooff=float(parasign("ooff",0,/nowarn))        ;opt offset added to output
  oslope=double(parasign("oslope",0,/nowarn))   ;opt slope added to output
  oquad=double(parasign("oquad",0,/nowarn))     ;opt quadrat. term added to output 
  if dflag ge 2 then nowarnflag=0 else nowarnflag=1 ; print warn-texts only for dflag >=2
  if paratest("oslope") or paratest("oquad ") then tozero=double(parasign("tozero",tarr[0],nowarn=nowarnflag, warn="WARN: 'oslope' or 'oquad' is set, but 'tozero' isn't. Offsets are calculated relative to first time-point")) else tozero=0D

  wflag = fix(parasign("wflag",0,/nowarn))
  tmpstr=(parasign("lchead",0,/nowarn))    ;check if setupfile (header) should be preprended to lc.
  case 1 of                                ;to deal with old setup files were lchead was specified with y / n 
     strupcase(tmpstr) eq 'Y': lchflag =1 
     strupcase(tmpstr) eq 'N': lchflag =0
     else: lchflag = fix(tmpstr) ;a number is given
  endcase
  if wflag eq 1 then begin      ;write header of output-lc
     lcoutname=parasign("outlcname","model.lc",/nowarn)
     if paratest("lcfilename") then lcoutname=parasign("lcfilename") ;compatibility with old param-name
     if lchflag then begin
        savesetupfile,lcoutname,/comm
        openw,unout,lcoutname,/get_lun,/append
        printf,unout,'#  '
        printf,unout,'# simulated lightcurve follows'
     endif else openw,unout,lcoutname,/get_lun
     printf,unout,'#   time ('+tunit+') '+ lumstr
  endif

;get total system luminosity
  contlum=float(parasign("contlum",0,/nowarn)) ;optional additional light-source to be added to syslum0
  syslum0 = contlum
  for i=0,n_object-1  do if obj_isa(obj[i],'star') then syslum0+=obj[i] -> getlum()
  
  lasttrflag=bytarr(n_object)   ;flags if there was a transit in last loop

  finestep= fix(parasign("finestep",'1',/nowarn)) ;do finestepping in gettrlump

;START MAIN LOOP
  onoiseed=frac(systime(1)/10000)*10000 ;seed for adding of noise (onoise)  
  tcount=0L
  for tcount=0L,nlcpts-1 do begin
     t=tarr[tcount]
     syslum=contlum                                         ;reset syslum(t)
     for i=0,n_object-1 do pos[i,*] = obj[i] -> getpos()    ;returns position of objects
     posarr[*,*,tcount]=pos[*,*]                            ;keeps xyz-positions of all objects

     if dflag ge 4 then for i=0,n_object-1 do obj[i] -> printpos
     zind=occultobid(sort(pos[occultobid,2]))  ;sequence with IDs of occulters, with farthest obj (-z) first
     for i=0,n_occult-1 do begin
         izs=zind[i]                               ;current index within zind array
        if obj_isa(obj[izs],'star') then begin    ;look for other objects in front of it
                                ; get array 'frind' with indexes of the objects that are in front
           frind= obj[izs] -> findinfront(i,zind,pos,rad)
                                ;frind contains either -1 (none in front) or the indexes
           if frind[0] eq -1 then begin  
              syslum+=obj[izs] -> getlum() ;no transits 
              if lasttrflag[i] eq 1 then if dispstarflag then obj[izs] -> displformoc
              lasttrflag[i]=0
           endif else begin     ;transit case
              if dflag ge 4 then begin
                 print,'in front of obj ',str(izs),' are obj(s) ',str(frind)
                 if anlycalcflag eq 0 then print,'anlyokflags: star: ',str(anlyokflag[izs]),' occulters: ',anlyokflag[frind]
              endif
              if n_elements(frind) eq 1 and (anlyokflag[izs] eq 1 and anlyokflag[frind[0]] eq 1 ) then begin
                 if dflag ge 4 then print,'analytical occult calc'
                 syslum+=obj[izs] ->gettrluma(obj[frind],occultcelflag=occultcelflag)    ;analytical calc
              endif else begin                               ;do pixel-calc
                 if anlycalcflag gt 0 then begin 
                    print,'UTM: anlycalcflag = 1 but encountered condition incompatible with analytical calculation.'
                    print,'     (Most likely multiple occulters on a star)'
                    stop
                 endif
                 if dflag ge 4 then print,'pixel-based occult calc'
                 syslum+=obj[izs] ->gettrlump(obj[frind],dispstarflag=dispstarflag,finestep=finestep) ;orig UTM command
              endelse       
              lasttrflag[i] = 1
              if dispstarflag then obj[izs] -> displformoc
              obj[izs] -> resetlformoc
           endelse
        endif
     endfor     
     if syslum +1e-10 lt syslum0 then intrarr[tcount]=1 ;set to 1 if any eclipsing happens. 1e-10 offset to avoid supurious in-transit determination
;    rloss=(syslum0-syslum)/syslum0
;   definitions of lumfin (output luminosities)    
     case oflag of
        0:lumfin=syslum                         ;system luminosity
        1: lumfin=(syslum-syslum0)/syslum0      ;relative flux-variation
        2: lumfin=(-2.5)*alog10(syslum/syslum0) ;magnitude changes
     endcase
     lumfin=lumfin+ooff +oslope*(t-tozero)+oquad*(t-tozero)^2    ;add offsets, slope, parabole
     if onoise ne 0.0 then begin                                 ;add noise to output
        lumfin=lumfin+randomu(onoiseed,/normal)*onoise
     endif
     if (tflag eq 2 or tflag eq 3) then begin ;combine input data and model
        case oflag of
           0: begin                                                 ;output in system luminosity or absolute flux
              if (tflag eq 2) then lumfin=lum0arr[tcount]+lumfin    ; add lumfin to input lc
              if (tflag eq 3) then lumfin=lum0arr[tcount]-lumfin    ; subtract lumfin from input lc
           end
           1: begin                                                             ;output in relative flux units
              if (tflag eq 2) then lumfin=(1+lum0arr[tcount])*(1+lumfin) -1.    ; multiply by relative flux loss
              if (tflag eq 3) then lumfin=(1+lum0arr[tcount])/(1+lumfin) -1.    ; divide by relative flux loss
;              if (tflag eq 2) then lumfin=lum0arr[tcount]+lumfin ; add lumfin to input lc
;              if (tflag eq 3) then lumfin=lum0arr[tcount]-lumfin ;   subtract lumfin from input lc           
           end
           2: begin                                                 ;output in mag variations
              if (tflag eq 2) then lumfin=lum0arr[tcount]+lumfin    ; add lumfin to input lc
              if (tflag eq 3) then lumfin=lum0arr[tcount]-lumfin    ; subtract lumfin from input lc
           end
        endcase
     endif

     lumarr[tcount]=lumfin      ;define the output array
     
     if dflag eq 1 or dflag ge 4 then begin
        if oflag eq 1 then print,format='("t= ",f16.8," ",a,"= ",e12.5)',t,lumstr,lumfin $
        else print,format='("t= ",f16.8," ",a,"= ",f15.5)',t,'brightn',lumfin
     endif

     if plotorbflag eq 1 then begin                                  ;plot in rescaled display 
        for i=0,n_object-1 do maxdist = (maxdist > max(abs(pos)))    ;keeps max distance of objects from coord-ctr  
        plotorb,posarr,tcount,maxdist,sunit,orbtrflag                ;orbital plane plots
     endif
     if plcsimflag and tcount ge 1 then plotlcutm,tarr[0>(tcount-100):tcount],lumarr[0>(tcount-100):tcount],tunit,oflag,tflag,lunit,syslum0,lum0arr[0>(tcount-100):tcount] 
     if dflag ge 4 then begin 
        tmp='' & read,tmp,prompt='--------hit return to continue--------------'
     endif 
     if tcount lt nlcpts-2 then begin ;position increment
        for i=0,n_object-1 do obj[i] -> incrempos,(tarr[tcount+1]-t),obj
     endif    
  endfor                        ;t MAINLOOP


;rebin output if resampling was done
  keep_resamp_flag = fix(parasign('keep_resample',0,/nowarn))
  if resample gt 1 then if not(keep_resamp_flag) then begin
     nlcpts_bin=nlcpts/resample
                                ;  tarr_bin=rebin(tarr,nlcpts_bin)
     lumarr_bin=rebin(lumarr,nlcpts_bin)
     lum0arr_bin=rebin(lum0arr,nlcpts_bin)
     intrarr_bin=byte(rebin(float(intrarr),nlcpts_bin)+0.99999) ; to get 1 in the binned array if just one element in resampled array is 1.
     nlcpts=nlcpts_bin                                          ;rename
                                ; tarr=tarr_bin
     tarr=tarr_org              ;UTM should never modify the tarr array, as it would mess up UFIT
     lumarr=lumarr_bin
     lum0arr=lum0arr_bin
     intrarr=intrarr_bin

  endif

  if wflag then begin           ;write outlc
     for i=0,nlcpts-1 do printf,unout,format='(f14.6," ",g13.7)',tarr[i],lumarr[i]
     free_lun,unout
  endif

  if plcflag then begin         ;final plot of entire lightcurve
     plotlcutm,tarr,lumarr,tunit,oflag,tflag,lunit,syslum0,lum0arr,/fullc
  endif                         ;plcflag

  if xpflag then begin          ;final elongation plot of x-positions
     if resample gt 1 then if not(keep_resamp_flag) then print,'WARN: final plot of x-positions not supported for resampling with final rebinning.' ; it will not work since the posarr has wrong dimension
     radarr=fltarr(n_object)
     for i=0,n_object-1 do radarr[i] = obj[i] -> getrad() ;get radii of all objects
     pxp= plot(tarr,posarr[0,0,*]-radarr[0],yrange=[min(posarr[*,0,*]),max(posarr[*,0,*])],xrange=[min(tarr),max(tarr)],xstyle=1, xtitle='time ('+tunit+')',ytitle='x-pos ('+sunit+')',title='Elongation (X-position)',dimensions=[800,320],location=[400,450])
     pxp=plot(tarr,posarr[0,0,*]+radarr[0],/over) ;first obj + radius

     for i=1,n_object-1 do begin
        pxp=    plot(tarr,posarr[i,0,*]-radarr[i],/over) ;other objects
        pxp=    plot(tarr,posarr[i,0,*]+radarr[i],/over)
     endfor
  endif                         ;xpflag

  case oflag of                   ; off-eclipse value of output model for utm output param
     0: modeloffvalue= syslum0    ;+ ooff ;sys luminosity
     1: modeloffvalue = 0.        ;ooff
     2:  modeloffvalue = 0.       ;ooff
  end 
;  endif                         ;oobflag = 1
  endwithoob: if oobflag then begin
     if resample gt 1 then begin ;need to do this as else larger arrays would be returned
        tarr=tarr_org
        lumarr=lum0arr_org
     endif
     if dflag ge 1 then print,'UTM early end with oobflag raised'
     paradd,'oobflag', '1'
  endif                                  ;oobflag
  if isa(obj) then obj_destroy,obj       ;kill object array
  !x= xsysorg &!y= ysysorg & !p= psysorg ;reset graphics sys variables 
end

