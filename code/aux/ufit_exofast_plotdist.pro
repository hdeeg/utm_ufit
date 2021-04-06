;+
; NAME:
;   ufit_exofast_plotdist
; PURPOSE:
;   Plots the distributions from MCMC chains (Probability Distribution
;   Functions (PDFs); returns median values and 68%
;   confidence intervals.
;
; DESCRIPTION:
;
; CALLING SEQUENCE:
;   PLOT = ufit_exofast_plotdist( pars [,MEDIANPARS=, ANGULAR=, PARNAMES=, UNITS=,
;                 BESTPARS=, PROBS=, /DEGREES], /FOURSIGMA, NBINS= )
;
; INPUTS:
;   PARS       - Array of NPARS x NSTEPS for each step in a Markov
;                Chain (with all chains concatenated together).
;
; KEYWORDS:
;   ANGULAR    - An array that indexes the parameters that are
;                angular (radians unless /DEGREES is set). This will enable
;                special ways to calculate the median and plots, which
;                may otherwise fail. Default is none.
;   PARNAMES   - An NPARS string array containing the parameter names
;                to be plotted on the axes.
;                HINT: TeXtoIDL is awesome.
;   UNITS      - An NPARS string array containing the units to be
;                displayed on the X-axis of the PDFs. 
;   BESTPARS   - An NPARS input array that specifies the "best" parameters
;                of each value. If specified, a red line will be drawn on
;                each plot at these values (this doesn't have to be
;                the best values).
;   NOPLOTMEDIAN -If set, plotting of the median values and errors is suppressed. Else, 
;                 median, +err, -err (68% confidence interval) of the distribution will be 
;                 shown as vertical solid (median) and dashed (errors) lines.   
;   PROBS      - Array of probability contours to draw on the covariance
;                plots. Default is erf([1d,2d]/sqrt(2d0)) ~ [0.68,0.95].   
;   DEGREES    - If set, ANGULAR parameters are in degrees, not radians.
;   FOURSIGMA  - If set, the x-range will extend to 4times 
;                the +- 1 sigma values (This is the default behaviour of 
;                the original exofast_plotdist code)
;   NBINS      - If set, the default of 100 bins per histogram is overrriden
;   NCOL       - IF set, the default of 2 plots side-by-side is overriden 
;   EXTRA      - Any further keywords accepted by PLOT, to set font-size, colours, tickmarks etc.

; OUTPUTS:
;   PLOT       - plot-structure of distribution plot
;   MEDIANPARS - A 3 x NPARS array of median values, +err, -err (68%
;                confidence interval).
;
;; Example:
;   Try the main-level example program at the end of this file::
;     IDL> .run ufit_exofast_plotdist.pro
;
; REVISION HISTORY:
;   2009/11/24 - Written: Jason Eastman - The Ohio State University
;   2015/09/16 - H.Deeg: version for UFIT. 
;   2019/03/04   H.Deeg: Major modifications: Changed into a function using newer IDL graphics.
;                Removed covariance plots (scattermatrix.pro should be used instead) 
;                Removed options for output file names; use .save method instead, see example program.
;                Added _EXTRA, NBINS and NOPLOTMEDIAN keywords.
;                Added main-level example program at end.
;   2019/03/08  Added NCOL keyword
;   2019/03/19   Added FOURSIGMA keyword and changed default to plot 
;               full range of values
;
; TODO: use more sophisticated panel positonining, like in scattermatrix.pro
; code for the angular keyword has not been revised over orignal exofast code 

function ufit_exofast_plotdist, pars, medianpars=medianpars, angular=angular, $
                                parnames=parnames, units=units, $
                                bestpars=bestpars, probs=probs, degrees=degrees,nbins=nbins, ncol=ncol, _extra=e

  nparnames = n_elements(parnames)
  nunits = n_elements(units)

;; 68% and 95% probability contours
  if n_elements(probs) eq 0 then probs = erf([1d,2d]/sqrt(2d0))

  sz = size(pars)
  if sz[0] ne 2 then message, 'ERROR: PARS must be an NPARS x NSTEPS array'
  npars = sz[1]
  nsteps = double(sz[2])

  if nparnames ne npars then begin
     parnames = strarr(npars)
     if nparnames ne 0 then message, /continue, $
                                     'WARNING: PARNAMES does not match parameter array; ignoring labels'
  endif

  if nunits ne npars then begin
     units = strarr(npars)
     if nunits ne 0 then message, /continue, $
                                  'WARNING: UNITS does not match parameter array; ignoring labels'
  endif

  medndx = nsteps/2
  halfsigma = erf(1d/sqrt(2d0))/2d
  lowsigndx = round(nsteps/2.d0 - nsteps*halfsigma)
  hisigndx = round(nsteps/2.d0 + nsteps*halfsigma)
  medianpars = dblarr(3,npars)


;; for angular parameters, center about the mode first
  nangular = n_elements(angular)
  for i=0, nangular-1 do begin

     hist = histogram(pars[angular[i],*],nbins=100,locations=x,/nan)
     max = max(hist,modendx)
     mode = x[modendx]

     ;; choose degrees vs radians
     if keyword_set(degrees) then halfrange = 180d0 $
     else halfrange = !dpi

     toohigh = where(pars[angular[i],*] gt (mode + halfrange))
     if toohigh[0] ne -1 then pars[angular[i],toohigh] -= 2.d0*halfrange

     toolow = where(pars[angular[i],*] lt (mode - halfrange))
     if toolow[0] ne -1 then pars[angular[i],toolow] += 2.d0*halfrange

  endfor

;number of panels to put side-by side

 ;number of panels to put side-by side
if not(keyword_set(ncol)) then ncol=2  ;number of columns
  nrow=(npars-1)/ncol +1        ;number of rows

  if isa(nbins) eq 0 then nbins =100      ;set default of 100 bins
  if isa(font_size) eq 0 then font_size=8 ;set default font smaller than the usual 12
  current=0                               ;plot first panel in new window
  for i=0, npars-1 do begin

     ytitle='Probability'

     ;; create the x axis titles
     xtitle = parnames[i]
     
     ;; add the units
     if units[i] ne '' then xtitle += ' (' + units[i] + ')'
     
     bad = where(~finite(pars[i,*]),complement=good)
     if bad[0] ne -1 then message, $
        "ERROR: NaNs in " +  parnames[i] + " distribution"    

     ;; get the median and +/- 1 sigma values
     sorted = sort(pars[i,*])
     medianpars[0,i] = pars[i,sorted[medndx]]
     medianpars[1,i] = pars[i,sorted[hisigndx]] - medianpars[0,i]  ; + error
     medianpars[2,i] = medianpars[0,i] - pars[i,sorted[lowsigndx]] ; - error

if keyword_set(foursigma) then begin; original exofast code: plot to +- 4 times +-1 sigma distance
     xmax = (medianpars[0,i] + 4*medianpars[1,i]) > min(pars[i,*],/nan)
     xmin = (medianpars[0,i] - 4*medianpars[2,i]) < max(pars[i,*],/nan)
endif else begin  ; plot full xrange
 xmax=max(pars[i,*]) 
xmin=min(pars[i,*])   
  endelse

     if xmin eq xmax then begin
        print, 'WARNING: Parameter ' + strtrim(i,2) + ' (' + parnames[i] + $
               ') is singularly valued; this is often a symptom of a bigger problem. Skipping covariances, which would fail.'
                                ;nocovar = 1
     endif

     ;; plot the probability distribution
     hist = histogram(pars[i,sorted],nbins=nbins,locations=x,min=xmin,max=xmax,/nan)
     distplot=plot( x, hist/double(nsteps), xtitle=xtitle, ytitle=ytitle,xstyle=1, xrange=[xmin,xmax],$
                    font_size=font_size,/histogram,/stairstep, layout=[ncol,nrow,i+1], current=current,_extra=e)
     
     ;; if the best parameters are given, overplot them on the distributions as red line
     ;; give the confidence interval around them

     if not(keyword_set(NOPLOTMEDIAN)) then begin ;,plot median and +- 68% confidences
        distplot=plot([medianpars[0,i],medianpars[0,i]],distplot.yrange,/over)
        distplot=plot([medianpars[0,i]+medianpars[1,i],medianpars[0,i]+medianpars[1,i]],distplot.yrange,linesty='dash',/over)
        distplot=plot([medianpars[0,i]-medianpars[2,i],medianpars[0,i]-medianpars[2,i]],distplot.yrange,linesty='dash',/over)
     endif
     if isa(bestpars) then begin ;plot 'best' value
        distplot=plot([bestpars[i],bestpars[i]],distplot.yrange,'r',/over)
     endif
     current=1
  endfor
  return, distplot
end




; main-level example program
;generates npar distributions, each one consisting of 2 gaussion distribs with ndat points and calls the plotting routine
npar=7
npts=100000  ;number of points 

ndat=npts/2  ;number of pts in each gaussian distrib
pname='par'+strcompress(string(indgen(npar)+0),/remove_all)
seed=!null
mult1=10^(randomu(seed,npar))  ;npar multipliers between 1..10 used for sdev
mult2=10^(randomu(seed,npar))  ;
offs1=(randomu(seed,npar)-0.5)*100  ;npar offsets from -50,50,
offs2=offs1+(randomu(seed,npar)-0.5)*10  ;npar offsets added to previous offs
multall=10^(randomu(seed,npar)*10-5)  ;npar multipliers between 10e-5, 10e5, used for entire distrib, to test for issues of very small/large numbers

datarr=fltarr(npar,2*ndat)
bestpars=fltarr(npar)  ; for 'best parameters'
for i=0,npar-1 do begin
k=multall[i]*[randomn(seed,ndat)*mult1[i] +offs1[i],randomn(seed,ndat)*mult2[i] +offs2[i]]
bestpars[i]=multall[i]*(offs1[i]+offs2[i])/2.  ;we use the center between the two distibs.
datarr[i,*] = k
endfor

medianpars=0.
dplot=ufit_exofast_plotdist(datarr, parnames=pname,bestpars=bestpars, nbins=50, medianpars=medianpars) ;plots of distributions
savename='example_distributions.png'
dplot.save,savename
print,'saved plot to file ',savename
  for i=0,npar-1 do begin
     print,pname[i],': best:',bestpars[i],' distr.mean:',medianpars[0,i],' +',medianpars[1,i],' -',medianpars[2,i]
  endfor
end


