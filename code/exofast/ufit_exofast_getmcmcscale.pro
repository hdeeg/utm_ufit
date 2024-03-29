;+
; NAME:
;   exofast_getmcmcscale
; PURPOSE:
;   Returns the optimal scale for MCMC chains by varying it until the
;   delta chi^2 = 1, or until delta chi^2 = DCMULT
;
; PROCEDURE:
;   Calculates the chi^2 of the best fit model. Then, for each
;   parameter, takes a small positive step (defined by seedscale) and
;   re-calculates the chi^2. It brackets the value that yields delta
;   chi^2 = 1 by doubling the seedscale until delta chi^2 > 1, then
;   does a binary search to find the value that is as close as
;   possible to delta chi^2 = 1. It repeats with a small negative
;   step. The average of the positive and negative excursion is
;   returned as the optimal scale.
;
;   If the best fit is accurate and the errors are
;   gaussian, uncorrelated, and accurately describe the data
;   (chi^2/dof = 1), this will result in the optimal acceptance rate
;   for MCMC fitting, 44% for 1 parameter, approaching 23% as the
;   number of parameters is large.
;
; INPUTS:
;   BESTPARS   - An NPARS element array of the best fit parameters
;   CHI2FUNC   - A string of the named function that calculates the
;                chi^2 of your model given the parameters specified by
;                BESTPARS.
;
; OPTIONAL INPUTS:
;   SEEDSCALE  - An NPARS arrray that contains the small step that is
;                initially taken. If the chi^2 surface is a smooth,
;                monotonically increasing surface, then it only
;                affects execution time of the program. Since the
;                chi^2 surface is rarely so forgiving, it's better to
;                err on the side of too small. The default is 1d-3 for
;                each parameter.
;   BESTCHI2   - The chi^2 of the best fit parameters. If not given, it
;                is calculated from BESTPARS using CHI2FUNC.
;   ANGULAR    - An array of indices (with up to NPARS elements)
;                indicating which elements are angular values. If
;                specified, and the step sizes determined are larger
;                than 2*pi, the step size will be set to pi.
;                NOTE: angular values must be in radians.
;   DCMULT  - changes the stopping criterium from delta chi2=1 to 
;                 delta chi^2= DCMULT. DCMULT should be larger than 1
;
; OPTIONAL KEYWORDS:
;  DEBUG       - If set, it will print progress to the screen that
;                will help identify problems.
;
; REVISION HISTORY:
;   2009/11/23 - Written: Jason Eastman - The Ohio State University
;   2013/02/26 - Change print statement to message statement
;   2015/09/16 -  H.Deeg: version for ufit: repeat statement at end: lt 1d-4 
;   2019/02/22 -  H. Deeg  modified end conditions in last repeat statement
;                 to 1d-3 and in line 159 to 1d-7. 
;   2019/03/09 -  H. Deeg: added keyword DCMULT to avoid the problem of 
;                stopping criterium deltachi2=1 often being too fine for UTM
;   2019/04/28 -  H. Deeg: added hook to deal with NAN chis-squares from out-of-bounds condition
;   2020/02/19 -  H. Deeg: added more useful feedback if unconstrained scales encountered
;   2020/07/27 -  H. Deeg: added autoseedscale kw for better seedscales of bestpar*1d-3 (instead of constant 0.001)
;                reversed conditon (l 171) for smallest scale-values from 1d-7 back to original 1d-14 
;                mod repeat statement at end (l 218): lt 1d-1  (from previous 1d-3; orig. 1d-7) for precision in reaching deltachi2= dcmult         

 
function ufit_exofast_getmcmcscale, bestpars, chi2func, tofit=tofit, $
                                    seedscale=seedscale, bestchi2=bestchi2,$
                                    angular=angular0, dcmult=dcmult, debug=debug, autoseedscale=autoseedscale
;common disp_vars, dflag 
  if not(keyword_set(dcmult)) then dcmult=1d
  npars = n_elements(bestpars)
  if n_elements(tofit) eq 0 then tofit = indgen(npars)
  nfit = n_elements(tofit)
  if n_elements(seedscale) ne nfit then begin
     if ~keyword_set(autoseedscale) then seedscale = dblarr(nfit) + 1d-3 else seedscale=dblarr(nfit) + bestpars[tofit]*1d-3
  endif
  maxiter = 1d4
  if n_elements(bestchi2) eq 0 then $
     bestchi2 = call_function(chi2func, bestpars)
;; which values are angular values (error if step > 2*pi)
  if n_elements(angular0) eq 0 then angular = [-1] $
  else angular = angular0

  if keyword_set(debug) then begin
     print, "   Par           Minimum Step                " + $
            "Maximum Step                Value                       " + $
            "Chi^2                    Best Chi^2"
  endif

  mcmcscale = [[seedscale],[seedscale]]

  for i=0, nfit-1 do begin

     for j=0,1 do begin

        testpars = bestpars
        
        minstep = 0
        maxstep = 0
        niter = 0
        bestdeltachi2 = !values.d_infinity

        repeat begin

           chi2changed = 0
           
           ;; an infinite step size means it's not constrained
           if ~finite(bestpars[tofit[i]] + mcmcscale[i,j]) or $
              ~finite(bestpars[tofit[i]] - mcmcscale[i,j]) then begin
              message, "EXOFAST_GETMCMCSCALE: Parameter "+strtrim(tofit[i],2) + $
                       " is unconstrained. Check your starting conditions"
           endif            

           ;; add the offset to a parameter
           if j eq 0 then begin
              testpars[tofit[i]] = bestpars[tofit[i]] + mcmcscale[i,j]
           endif else begin
              testpars[tofit[i]] = bestpars[tofit[i]] - mcmcscale[i,j]
           endelse

           ;; determine the new chi^2
           chi2 = call_function(chi2func, testpars)
           ;;added by HD: if stepped into boundary, increment chi2 by moderate
           ;;size that permits determin. of useful midpoint in next steps
           if ~finite(chi2) then chi2=bestchi2+2.d0*dcmult    
           ;; determine the next step size based on deltachi2
           if (chi2 - bestchi2) ge 1.d0*dcmult then begin
              ;; if deltachi2 is too large, set max to current value
              maxstep = mcmcscale[i,j]
              ;; the next step will be the midpoint of min/max
              mcmcscale[i,j] = (maxstep + minstep)/2.d0
           endif else if (chi2-bestchi2 ge 0) then begin
              ;; if deltachi2 is too small, set min to current value
              minstep = mcmcscale[i,j]
              ;; if a bound on the max hasn't been determined double the step
              ;; otherwise, take the midpoint of min/max
              if maxstep eq 0 then begin
                 mcmcscale[i,j] *= 2
              endif else mcmcscale[i,j] = (maxstep + minstep)/2.d0
           endif else begin
              if keyword_set(debug) then begin
                 print,'WARNING: better chi2 found by varying parameter '+$
                       strtrim(tofit[i],2) + ' from ' + $ 
                       strtrim(string(bestpars[tofit[i]],format='(f40.10)'),2)+$
                       ' to ' + $
                       strtrim(string(testpars[tofit[i]],format='(f40.10)'),2)+$
                       ' (' + strtrim(chi2,2) + ')'
              endif
              
              ;; chi2 is actually lower! (Didn't find the best fit)
              ;; attempt to fix
              bestpars = testpars

              ;; could be way off, double the step for faster convergence
              mcmcscale[i,j] *= 2d0 
              bestchi2 = chi2
              chi2changed = 1
              niter = 0
           endelse

           deltachi2 = chi2-bestchi2          
           ;; in case we chance upon a better match than we bracket
           ;; (implies chi^2 surface is rough)
           if abs(deltachi2 - 1d0*dcmult) lt abs(bestdeltachi2 - 1d0*dcmult) then begin
              bestdeltachi2 = deltachi2
              bestscale = mcmcscale[i,j]
           endif

           ;; can't always sample fine enough to get exactly
           ;; deltachi2 = 1 because chi^2 surface not perfectly smooth
           if abs(minstep - maxstep) lt 1d-14 or niter gt maxiter then begin
;          if abs(minstep - maxstep) lt 1d-7 or niter gt maxiter then begin  ;this might be better for UTM
              if not chi2changed then begin
                 if abs(bestdeltachi2 - 1.d0*dcmult) lt 0.75*dcmult then begin
                    mcmcscale[i,j] = bestscale
                    chi2 = bestchi2 + 1*dcmult
                 endif else begin
                    message, 'Convergence Error: cannot find the value' + $
                             ' for which deltachi^2 = 1*dcmult',/continue
                    mcmcscale[i,j] = !values.d_nan
                    chi2 = bestchi2 + 1*dcmult
                 endelse
              endif
           endif
           
           ;; if the parameter has no influence on the chi2
           if abs(chi2 - bestchi2 - 1.d0) eq 1.d0 and niter gt maxiter then $
              message, 'ERROR: changing parameter ' + strtrim(tofit[i],2) + $
                       ' does not change the chi^2. Exclude from fit.' 

           ;; if angle is so poorly constrained 
           ;; no value has delta chi^2 = 1
           if (where(angular eq tofit[i]))(0) ne -1 and $
              minstep gt (2.d0*!dpi) then begin 
              message, 'WARNING: no constraint on angle',/continue 
              mcmcscale[i,j] =  !dpi
              chi2 = bestchi2 + 1.d0*dcmult
           endif

           ;; near a boundary
           ;; exclude this direction from the scale calculation
           if ~finite(chi2) then begin
              chi2 = bestchi2 + 1.d0*dcmult
              mcmcscale[i,j] = !values.d_nan
           endif
           niter++

           ;; print the steps and chi^2 it's getting
           if keyword_set(debug) then begin
              if j eq 0 then str = '(hi)' $
              else str = '(lo)'                
              print, tofit[i], str, minstep, maxstep, $
                     testpars[tofit[i]],chi2,bestchi2, $
                     format='(i3,x,a4,x,f27.18,f27.18,f27.18,f27.18,f27.18)'
           endif
                                ;  endrep until abs(chi2 - bestchi2 - 1.d0) lt 1d-8
;HD: MODIFIED THE LIMIT BELOW to to be less demanding
        endrep until abs(chi2 - bestchi2 - 1.d0*dcmult) lt 1d-2 *dcmult

     endfor    
  endfor
;; replace undefined errors with the other
  bad = where(~finite(mcmcscale[*,0]))
  if bad[0] ne -1 then mcmcscale[bad,0] = mcmcscale[bad,1]
  bad = where(~finite(mcmcscale[*,1]))
  if bad[0] ne -1 then mcmcscale[bad,1] = mcmcscale[bad,0]
  bad = where(~finite(mcmcscale),nbad)
  if nbad gt 0 then begin       ; Failed: parameter(s) unconstrained
     print,'ufit_exofast_getmcmcscale: Error: unconstrained scales encountered'
     print,' printing mcmcscale variable'
     print,format='(f19.8)',mcmcscale
     stop
                                ;return, -1 
  endif
;; return the average of high and low
  finalscale=total(mcmcscale,2)/2.d0
  return, finalscale
end
