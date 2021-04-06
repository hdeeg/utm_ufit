;pro ufit,setfname,disflag,sdev,a
;fitting program to call model-generating programs that read
;setup-files.
;This version is build to call utm (universal transit modeler) or freely definable modellers
;written by Hans J Deeg
; 
;version history
; 3/5/99 first working version
; 5/5/99 includes wfit handling
; 3/8/00 dbrat and chitol are now optionals in parameter file
; 30/9/00 clarified screen-output of chi and stdev
; 9/10/00 inserted output params a,var; dflag 0 doesn't plot anymore
; 18/10/00 changed repmax from 10 to 3 and itmax from 20 to 10
; 8/11/00 repmax and itmax are now optional paramfile values without 
; warning on absence. Temporary fitsets have now a unique
; filename. This allows to run more then 1 ufit process at the same
; time in the same directory.
; 9/11/00 fixed bug introduced yesterday
; 7/12/00 fixed minor bug preventing saving of model-lc if wflag=1
; 26/8/2005 added copyright and licence
; 19/9/2006 calcs reduced ChiSq and Q value for entire dataset and in-transit
; 9/5/2015 minor fix due to error from undefined onoise param 
; 12/5/15 name of output setupfile now related to input setup.
; 13-15/May/15  profound revision with adding of AMOEBA, changes of graphics and text output, scale-param, etc.
; 18may15 added nofit keyword to give only the starting chisq and plot.
; 21may15 added fitsetnamext  to give output setups a distinct name
; 22may15  changed max number of fit-params to 16. Untested on LM fit, but warns if >8.
; 23may15 changed dflag to go from 0..4 and reassigned levels of info displayed. Changed output param 'var' to sdev.
;24may15 undid fitsetnamex as it was confusing. Added savenextfinflag
;25may15  changed calc of stddev to use formula without the mean
;30may15 added graceful exit with saving of last fitsetmp as output if nmax is exceeded in AMOEBA. Also adds output params on number of amoeba calls and if an exceed happened.
;10jun15 added randominit to randomize intial values. Fit's standard deviation is written to output setup as 'sdev'
;11jun15  added version_nr for overriding of  version numbering.
;17jun15 fixed bug that fitted paramters weren't displayed correctly
;1jul15 changed ftol and nmax params to amoeba_ftol and amoeba_nmax (ftol remains recognized)
;4jul15 added ramdom_uniform keyword for uniformly distributed offsets to initial values.
;4sep15 replaced call to  ameoba by exofast_amoeba
;8sep15 cleanup: removed several duplicate arrays and variables. Inserted ndat, ymod into common block. Moved variables specific to fit-algos into common blocks am_vars, lm_vars
;9sep15a added meas-error to chisq-calc in amoeba and ufit in/output
;9sep15b moved chisq calc and stats-display stuff into function chisq_utm
;10sep15 added chisq to output-setupfile. Changed default fit-algorithm to Amoeba
;12sep15 several changes to display, added residlevel paramters
;16sep15 added mcmc fit
;17sep15  fix crash on  /nofit and dflag=1. Added mc_stopscale setup-param
;19sep15 MCMC keeps now also the best-fit parameters in output setupfile
;9oct15 fixed: bug on undefined savenextfinflag for non-standard setupfile-names; missing axes titles in final plot; reversing of y-axes for oflag=2 (delta-mags)
;15oct added inlcdigi to keep fewer digits in input-lc time-stamps, like in utm.
;23oct15 added AICc and BIC to stats output
;21feb19 incremented disflags 1 to 5 by 1 while disflag=1 now shows only initial and final statistics. Renamed disflag from debflag
;22feb19 implemented saving of scales from GETMCMCSCALE to setup. Moved dflag 
;      to separate common block disp_vars
;9mar19  added mc_dcmult setup-param to avoid the problem of 
;       default stopping criterium of delta chi^2 =1 in get_mcmcscale 
;        often being too fine. 
;17mar19 changed default mc_dcmult to 5
;11apr19  changed behaviour of plcflag to adapt to utm version of today
;26apr19  added framework for out of boundary (oobflag) conditions. Can be used
;         in UTM preprocessor to define uniform priors.
;28apr19  changed the out-of-boundary chisq to NAN for integration with demc.pro
;4may19   added initmult keyword for initial scatter of fit-parameters in MCMC
;9may19   fix, so that model lightcurves are written if /nofit is used and wflag =1 
;19feb20  amoeba handles now correctly out-of-bounds  (UTM raises oobflag) conditions
;15mar20  modified oobflag to be now handled only through pnarr, pvarr
;24mar20  moved execution of utm to be always within chisq_utm, and chisq_utm is now called indirectly through call_function
;10apr20  added framework for priors being defined from setupfile. 
;23jul20  added ability to execute other modeller routines than UTM (for now: ufit_polynome and similar ones). 
;28jul20  removed all LM-fit code. Moved ancilliary routines (some of them shared with UTM) into utm-ancil.pro.
;20aug20  ufit3: test version dealing with more complex indatafiles
;25aug20  added treatment of auxiliary columns in indatafile and options for reading from indatafile (delchar etc.) 
;28aug20  improved error-message if a fit parameter is not in the setup or numbers of scale and fit params don't agree
;8sep20   added _extra to chisq_utm and fixed plot-ranges for mag-display and with and without erroplot
;6feb20 marginal fix in output (removed duplicated printing of 'best fit:...')
;20mar21 added absolute difference prior and fixed error on one-sided priors

;TODO  
;   utm: permit individual yerrs   
;

;COPYRIGHT (C) 1999-2019 Hans J. Deeg;   
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    any later version.
;
;    As a special exception, this program may be compiled  with the
;    Interactive Data Language (IDL) and be linked to libaries
;    pertaining to IDL.
;
;    This program is distributed in the hope that it will be ueful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program in file 'licence.txt'; if not, write to
;    the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
;    Boston, MA 02110-1301 USA

@utm_ancil   ;compile anciliary routines
@setupfile ;compile setupfile



function chisq_model, print=print, sigall=rms_res, _extra=ex    
;returns chi-sq for non-UTM models. NANs in ymod are ignored. Use _extra to pass further keywords to modeller
  common sharxy,ndat,xdat,ymod,ydat,yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue,testlist
  common fitparnamecom,fitparname,fitsetmp,adim
  common parrs,pvarr,pnarr
  common disp_vars, dflag

  model=parasign("model")
  call_procedure,model,fitsetmp,dflag,xdat,ymod,ydat=ydat,yerr=yerr,auxdat1=auxdat1,auxdat2=auxdat2,auxdat3=auxdat3, _extra=ex   ;use of ydat, _extra (e.g. for display-stuff) depends on modeller
  yresid=ydat-ymod 

  oobflag=int(parasign("oobflag",0,/nowarn)) ;check if an oobflag was raised
  if oobflag then begin 
     rss=!VALUES.F_INFINITY     ;out of bound flag is raised
     chisq =!VALUES.F_INFINITY  
  endif else begin
gi=where(finite(ymod))  ;good IDs
ngood=n_elements(gi)
     rss=total(yresid[gi]^2)                  ;RSS, residuals sum square 
     chisq=total((yresid[gi]/yerr[gi])^2)
     rms_res=sqrt(rss/ngood)     ;root mean square deviation (close to =stddev if mean =0) of residuals
  endelse

  if keyword_set(print) then begin
     chisq_red=chisq/(ngood- adim-1) ;reduced chi-sq
     BIC=chisq+adim*alog(ngood)
     AIC=chisq+adim*2
     AICc=AIC+2*adim*(adim+1)/(ngood-adim-1) ; AIC corected for small sample size
     if oobflag then begin
        print,'WARN: oobflag: a parameter is out-of-bounds, chisq = INF'
     endif else begin
        print,'chisq (all pts)    :',chisq, ' avg.meas_err: ',avg(yerr)
        print,'chisq_red (all pts):',chisq_red 
        print,'AICc               :',AICc,  '   BIC       : ',BIC
        print,'rms_residuals      :',rms_res,   '   npt_all   : ',ngood
     endelse
  endif
  

;below block to write the fitted model to disk
  wflag=fix(parasign('wflag',0,/nowarn) )
  if wflag eq 1 then begin      
     lcoutname=parasign("outlcname","model.txt",/nowarn)
     openw,unout,lcoutname,/get_lun
     for i=0,n_elements(xdat)-1 do printf,unout,format='(f14.6," ",g13.7)',xdat[i],ymod[i]
     free_lun,unout
  endif

  return, chisq
end


function chisq_utm, sigall=sigall,sigontr=sigontr,sigofftr=sigofftr, print=print,chisqadd=chisqadd, _extra=ex   
; executes utm and
;returns chisq of utm (excluding priors) 
;optional calc and printing of on/off-transit stats
;input keyw:
;   print   generates stats display
;output keyw
;    sigall, sigontr, sigofftr  standard deviations of all pts, on-transit, off-transit
;    chisqadd   addition to chisq calculated in preprocessor (usually for priors)


  common sharxy,ndat,xdat,ymod,ydat,yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue
  common fitparnamecom,fitparname,fitsetmp,adim
  common parrs,pvarr,pnarr
  intrarr =0                    ;array needed for calc of  off-transit stddev
  utm,fitsetmp,0,xdat,ymod,intrarr=intrarr,modeloffvalue=modeloffvalue 

  yresid=ydat-ymod 
  oobflag=int(parasign("oobflag",0,/nowarn)) ;check if an oobflag was raised
  if oobflag then begin 
     rss=!VALUES.F_INFINITY     ;out of bound flag is raised
     chisq =!VALUES.F_INFINITY  
  endif else begin
     rss=total(yresid^2)                  ;RSS, residuals sum square 
     chisq=rss/yerr^2               ;simple chi-sq of all points
  endelse
  sigall=sqrt(rss/(ndat-1))     ;stddev

  if keyword_set(print) then begin
     chisq_red=chisq/(ndat- adim-1) ;reduced chi-sq
     BIC=chisq+adim*alog(ndat)
     AIC=chisq+adim*2
     AICc=AIC+2*adim*(adim+1)/(ndat-adim-1) ; AIC corected for small sample size
     if oobflag then begin
        print,'WARN: oobflag: a parameter is out-of-bounds, chisq = INF'
     endif else begin
        print,'chisq (all pts)    :',chisq, ' meas_error : ',yerr
        print,'chisq_red (all pts):',chisq_red 
        print,'AICc               :',AICc,  '  BIC       : ',BIC
        print,'stdev (all pts)    :',sigall,'  npt_all   : ',ndat
     endelse
  endif

  if keyword_set(sigofftr) then begin
     offtr_id=where(intrarr eq 0,offct)
     if offct ne 0 then begin   ; do these calcs only for curves with off-transit part
        sigofftr=sqrt(total(yresid[offtr_id]^2)/(offct-1)) 
        if keyword_set(print) then print,'stdev (off-eclipse):',sigofftr,'  npt_off   : ',offct
     endif else sigofftr=!values.f_nan
  endif                         ;keywd sigofftr

  if keyword_set(sigontr) then begin
     intr_id=where(intrarr eq 1,onct) ; indexes on and off-eclipse
     if onct ne 0 then begin ; do these calcs only for curves with  on and off-transit part
        sigontr=sqrt(total(yresid[intr_id]^2)/(onct-1))
;        chisq_ontr= total(yresid[intr_id]^2/sigofftr^2)
        chisq_ontr= total(yresid[intr_id]^2/yerr^2)
        chisq_red_ontr = chisq_ontr/(onct-adim-1) 
        if keyword_set(print) then begin
           print,'stdev (on-eclipse) :',sigontr,'  npt_on    : ',onct
           print,'chisq (on-eclipse) :',chisq_ontr 
           print,'chisq_red (on-ecl) :',chisq_red_ontr
        endif                   ;print
     endif  else sigontr=!values.f_nan
  endif                         ;keywd sigontr

  chisqadd=double(parasign("chisqadd",0,/nowarn,foundflag=foundflag)) ;addit value from chisqadd (priors calculated in preproc)
if foundflag eq 0 then chisqadd=[]   ;delete value  if chisqadd absent in setup
  return, chisq
end



function amoeba_utm,atest
; interface from call by amoeba to chisq_utm
;returns chisq of utm (lumarr=model) as fit-fct for amoeba. atest are current vals of fitted params
  common parrs,pvarr,pnarr
  common sharxy,ndat,xdat,ymod,ydat,yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue
  common fitparnamecom,fitparname,fitsetmp,adim
  common am_vars, itct, chisqmin
  common disp_vars, dflag
  common fitfunc, fitfunct, priorfunct

                                ; lumarr=ydat                   ;data from inlc
;update 'atest' into parameterfile
  for i=0,adim-1 do begin
     paraput,fitparname[i],atest[i]
  endfor
  savesetupfile,fitsetmp
  if dflag ge 4 then print,'AMOEBA_UTM  test params:',atest
  if dflag ge 1 then begin sigall=1. & sigontr=1. & endif     ;set to 1 to execute keywords
   chisqmod=call_function(fitfunct,sigall=sigall,sigontr=sigontr,chisqadd=chisqadd) ; returns chisq from utm-model
   chisq=call_function(priorfunct,chisqadd=chisqadd, chisqmod=chisqmod)
   if finite(chisq,/nan) then chisq = chisqmin+1e6 ;raise chisq if UTM raised oobflag

;DISPLAY STUFF
   if itct eq 0 then chisqmin = chisq                         ;initalize current best chisq
   chisqmin = chisqmin < chisq                                ;update it
   
   if dflag eq 2 or dflag eq 3 then begin
      print,format='($,i0," ")',itct  ;count of interations
      if chisq eq chisqmin then begin ;print stats only if new best value
         print,' '
         print,'AMOEBA_UTM new min params:',atest
         print,'AMOEBA_UTM new min: chisq,stddev,stddev(on-tr)', chisq,sigall,sigontr
         if dflag eq 3 then plotlcymod,/update        ;update disp with new best solu
      endif                                           ; chisq eq chisqmin
   endif                                              ;  dflag eq 2 or dflag eq 3
   if dflag ge 4 then begin                           ;display every iteration
      print,'AMOEBA_UTM itct, chisq,chi/chimin,stddev,stddev(on-tr):',  itct, chisq, chisq/chisqmin,sigall,sigontr
      if dflag ge 5 then plotlcymod,/update ;update disp with temp. solu
   endif
   itct+=1
   return,chisq
end


function mcmc_utm,atest,determinant=determinant
; interface from call by demc to chisq_utm
;returns chisq of utm (lumarr=model) as fct called by exofast_demc. atest are current vals of fitted params
  common parrs,pvarr,pnarr
  common sharxy,ndat,xdat,ymod,ydat,yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue
  common fitparnamecom,fitparname,fitsetmp,adim
  common am_vars, itct, chisqmin
  common disp_vars, dflag
  common fitfunc, fitfunct, priorfunct
; lumarr=ydat                   ;data from inlc
;update 'atest' into parameterfile
  for i=0,adim-1 do begin
     paraput,fitparname[i],atest[i]
  endfor
  savesetupfile,fitsetmp
if dflag ge 1 then begin sigall=1. & sigontr=1. & endif     ;set to 1 to execute keywords

   chisqmod=call_function(fitfunct,sigall=sigall,sigontr=sigontr,chisqadd=chisqadd) ; returns chisq from utm-model
   chisq=call_function(priorfunct,chisqadd=chisqadd, chisqmod=chisqmod)


;DISPLAY STUFF
   if itct eq 0 then chisqmin = chisq                         ;initalize current best chisq
   chisqmin = chisqmin < chisq                                ;update it
   if dflag ge 3 then begin                                   ;display every iteration
      print,'MCMC_UTM  test params:',atest
      print,'MCMC_UTM itct, chisq,chi/chimin,stddev,stddev(on-tr):',  itct, chisq, chisq/chisqmin,sigall,sigontr
;     if dflag ge 4 then plotlcymod ;update disp with temp. solu
   endif
   itct+=1
   determinant=1.               ;dummy value expected in exofast_demc
   return,chisq
end

pro plotlcymod,first=first,update=update, last=last,xtitle=xtitle,ytitle=ytitle,errplot=errplot ;plot data and fit
  common sharxy,ndat,xdat,ymod,ydat, yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue
  common  plotstruc, plc, residlevelflag, residlevel  
;from common block 
;xdat  (xdata)
;ydat  (ydata)
;ymod  (model from utm)
  yerrorg=yerr ;keep orig as this routine might modify yerr w undesirable effects
  resid=ydat - ymod             ;this is a residual centered on zero
  
  if keyword_set(residlevelflag) then begin ; a level for the residuals was given in setup
     residmean = mean(resid)
     resid -= (residmean - residlevel) ;average the residuals around residlevel
  endif else begin
     resid+=modeloffvalue       ;add off-eclipse model-value to residual 
  endelse


  if keyword_set(errplot) then begin
     if n_elements(yerr) eq 1 then yerr=replicate(yerr,ndat) ;if meas-err is used, set all yerr to the same error
  endif
  if ~isa(yerr) then yerr=0.0   ;use dummy for the yrange-calc

  if keyword_set(errplot) then begin ;errobars plotted
     if oflag ne 2 then begin        ;common case
        ymin=min([ydat,ymod-yerr,resid])
        ymax=max([ydat,ymod+yerr,resid])
     endif else begin           ;delta-mags; revert y-scale
        ymax=min([ydat,ymod-yerr,resid])
        ymin=max([ydat,ymod+yerr,resid])
     endelse
  endif else begin              ;no erropars plotted
     if oflag ne 2 then begin   ;common case
        ymin=min([ydat,ymod,resid])
        ymax=max([ydat,ymod,resid])
     endif else begin           ;delta-mags; revert y-scale
        ymax=min([ydat,ymod,resid])
        ymin=max([ydat,ymod,resid])
     endelse
  endelse

;extend yrange by 3%=1/33 above, below
        ymin=ymin-(ymax-ymin)/33.
        ymax=ymax+(ymax-ymin)/33.

                                ;finding the offset to residual by searching for most common value of ymod (= off-eclipse)
;does not work realiably of there is no or very few  off-transit points in lc
;     hist=histogram(ymod,nbin=n_elements(ymod),locat=locat)
;     tmp=max(hist,max_id)       ;get indize of most frequent bin in hist                            
;     offset=locat(max_id)
;     resid+=offset              ;add to residual
;  endelse

  if keyword_set(first) then begin    
                                ;inital plot with input data
     if ~keyword_set(errplot) then begin ; no error-bars
        plc=  plot(xdat,ydat,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax]) ;input data as crosses
     endif else begin ; with error-bars                                                                      ;do erroplot
        plc=errorplot(xdat,ydat,yerr,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax])
     endelse
  endif
  if keyword_set(update) then begin ;add intermediate models
     plc.yrange=[ymin,ymax]
     plc=plot(xdat,resid,':',over=plc)       ;resids as dot-line
     plc=plot(xdat,ymod,over=plc)            ;model solid line
  endif
  if  keyword_set(last) then begin ;do final plot in color
     plc.yrange=[ymin,ymax]                                                 
     plc=plot(xdat,resid,'D',sym_color='CRIMSON',sym_thick=1,sym_size=0.7,over=plc) ;resids as red squares
     plc=plot(xdat,ymod,'-g',thick=2,over=plc) ;final model as thick line
  endif
  yerr=yerrorg                  ;revert to orig
end

function genprior::init,prinum,prpnid,prparmax ;generic prior initialization 
;assignment of variable name and position with validity checks) 
  common parrs,pvarr,pnarr                     ;vars for setupfile 
  self.prinum=prinum
  foundflg=0  
     if pnarr[prpnid+1] eq 'prvar' then begin ;assign the var-name and check if existing
        self.prvar=pvarr[prpnid+1]
        foundflg = paratest(self.prvar,parid=parid) ;check if  prvar is in setup
        if foundflg eq 0 then message,'UFIT: prior nr. '+str(prinum+1)+' is not corresponding to any parameter: prvar = '+self.prvar 
        self.prvid=parid        ;index of prior-variable in setup
     endif else message,"UFIT: setup of prior nr. "+str(prinum+1)+"'prtype' needs to be followed by 'prvar'"
  return,foundflg
end

function gprior::init,prinum,prpnid,prparmax  ;init gauss prior and also abs-difference prior
  common parrs,pvarr,pnarr      ;vars for setupfile 
  ival=self -> genprior::init(prinum,prpnid,prparmax)
  foundflg=ival
  self.prsign=!VALUES.D_INFINITY   ;initiualize as equivalent to uniform prior
  self.prsigp=!VALUES.D_INFINITY 
  onesidedflg=0                 ;flag for onesided sigma
  for i=prpnid+2, prpnid+prparmax+1 do begin 
     case pnarr[i] of
        'prctr': begin
           self.prctr=double(pvarr[i])
           foundflg+=1    
        end
        'prsigp': begin
           self.prsigp=double(pvarr[i])
           onesidedflg+=1 
           foundflg+=1    
        end
        'prsign': begin
           self.prsign=double(pvarr[i])
           onesidedflg+=1 
           foundflg+=1    
        end
        'prsig': begin          ;set positive and neg. sigmas to the same
           self.prsign=double(pvarr[i])
           self.prsigp=double(pvarr[i])
           foundflg+=2    
        end
        else: dummy=0             ;ignore statement
     endcase
  endfor
  if onesidedflg eq 1 then foundflg+=1 ;if either prsigp or prsign given, don't need the other
  if foundflg ne 4 then message,'UFIT: prior nr. '+str(prinum+1)+' ('+self.prvar+'): incorrect parameter-set is given, revise setup' else ival =1
  return,ival
end

function uprior::init,prinum,prpnid,prparmax ; init uninformed prior
  common parrs,pvarr,pnarr                   ;vars for setupfile 
  ival=self -> genprior::init(prinum,prpnid,prparmax)
  foundflg=ival
  onesidedflg=0                 ;flag for onesided sigma
  limits=0
  for i=prpnid+2, prpnid+prparmax+1 do begin 
     case pnarr[i] of 
        'prloli': begin
           self.prloli=double(pvarr[i])
           onesidedflg+=1 
           limits=1             ;lower limit
           foundflg+=1    
        end
        'prhili': begin
           self.prhili=double(pvarr[i])
           hilimit=1
           limits=2             ;upper limit
           onesidedflg+=1 
           foundflg+=1    
        end
        else: dummy=0           ;ignore statement
     endcase
  endfor
  if onesidedflg eq 1 then begin ;define the other side as Nan
     case limits of
        1: begin
           self.prhili=!VALUES.D_NAN
           foundflg+=1         
        end
        2:  begin
           self.prloli=!VALUES.D_NAN
           foundflg+=1         
        end
        else: dummy=0
     endcase
  endif
  if foundflg ne 3 then message,'UFIT: prior nr. '+str(prinum+1)+' ('+self.prvar+'): incorrect parameter-set is given, revise setup' else ival =1
  return,ival
end

function genprior::chisqprior
;returns chisq of a prior
return, 0D
end

function gprior::chisqprior
;returns chisq of a symetric or asymetric gaussian prior
  testval=double(parasign(self.prvar)) ;get value of variable 
  if self.prsigp eq self.prsign then chisqpr= ((testval-self.prctr)/self.prsigp)^2 else begin
     if testval lt self.prctr then chisqpr= ((testval-self.prctr)/self.prsign)^2 else $
        chisqpr= ((testval-self.prctr)/self.prsigp)^2
  endelse
  return, chisqpr
end


function aprior::chisqprior
;returns chisq of a symetric or asymetric absolute-difference  prior
  testval=double(parasign(self.prvar)) ;get value of variable 
  if self.prsigp eq self.prsign then chisqpr= abs((testval-self.prctr)/self.prsigp) else begin
     if testval lt self.prctr then chisqpr= (self.prctr-testval)/self.prsign else $
        chisqpr= (testval-self.prctr)/self.prsigp
  endelse
  return, chisqpr
end


function uprior::chisqprior
;returns chisq of a prior
  testval=double(parasign(self.prvar)) ;get value of variable 
  case 1 of ;; also if self.prloli or hili are NaN, cases are not selected
     testval ge self.prhili: chisqpr=!VALUES.F_INFINITY
     testval le self.prloli: chisqpr=!VALUES.F_INFINITY
     else: chisqpr=0D  
  endcase
  return, chisqpr
end

function chisqallprior,chisqadd=chisqadd, chisqmod=chisqmod, print=print
  common priors,prior
;calculates chisq of all defined priors and adds it to any pre-existing chisq (chisqmod)
; keyw:
; chisqadd  : additional  chisq from modeller program (typically a chisq from a prior in preprocessor) that is added to output chisq 
; chisqmod : the chisq of the model-fit, which is added to output chisq 
; print  : prints chisquares ONLY IF some prior is defined

  if isa(prior) or isa(chisqadd) then priorflag =1 else priorflag=0 ;ANY priors defined
  if priorflag then begin       ;any priors are there
     if isa(chisqadd) then chisqpri=chisqadd else chisqpri=0D ; add priors from model preproc 

     if isa(prior) then begin   ;a prior was defined in setup

        for i=0,n_elements(prior)-1 do begin
           chisqpri += prior[i].chisqprior() ;chisq of priors
        endfor
     endif


     if keyword_set(print) then  print,'chisq (priors)     :',chisqpri
     if isa(chisqmod) then begin ;add chisq of priors to chisq of model
        chisqtot=chisqpri+chisqmod
        if keyword_set(print) then print,'chisq (total)      :',chisqtot
     endif else chisqtot=chisqpri
  endif else begin   ;priorflag =0; use input model chisq 
     if isa(chisqmod) then chisqtot=chisqmod else chisqtot = 0D
  endelse
  return, chisqtot
end


;=================================================================================================
;=================== UFIT MAIN PROGRAM ===========================================================
;=================================================================================================
pro ufit,setfname,disflag,sdev,a,nofit=nofit

;parameters
;  setfname   name of setupfile
;optional paramters:
; ;disflag   level of graphic display, goes from 0 to 6, default 3
;  sdev     standard deviation (output param)
;  a        array of best fit parameters (output param)
;keywords 
;      nofit: no fit is made, only statistics of the input model is displayed

  common parrs,pvarr,pnarr                                                                            ;vars for setupfile
  common sharxy,ndat,xdat,ymod,ydat,yerr,auxdat1,auxdat2,auxdat3,oflag,intrarr,modeloffvalue,testlist ;principal data variables 
  common am_vars, itct, chisqmin                                                                      ;vars for Amoeba fit
  common lm_vars, yarr,xind, pdrat, alast                                                             ;vars for LM fit
  common fitparnamecom,fitparname,fitsetmp,adim                                                       ;related to fit params
  common  plotstruc, plc, residlevelflag, residlevel                                                  ;plot-structure
  common disp_vars, dflag
  common fitfunc, fitfunct, priorfunct ;fit and prior functions
  common priors,prior
  !except=0                     ;don't report math errors

  array1 = [2,3,4,5]
  array2 = [2.0,3.0,4.0,5.0]
  testlist = LIST(array1, array2)

;interpreting setup filename, extract eventual version numbers and prep for output setup-name
  readsetupfile,setfname
  setfnameflag= strmatch(setfname,'*.utm') ; is 1 if ending in .utm
  if setfnameflag then begin 
     setfnamecore0 = strmid(setfname,0,strlen(setfname)-4) ;extract initial core name
     setfnamecorele = strpos(setfnamecore0,'.f_in')        ;check for presence of .fin = fit-input
     if strpos(setfnamecore0,'.fout') ne -1 and not(keyword_set(nofit)) then begin ;a .fout file given
        print,'WARN: output setup ',setfname,' given as input.'
        tmps='' & read,'It will be overwritten by next output. Proceed (Y/N)? ',tmps
        if strupcase(tmps) ne 'Y' then stop    
     endif 
     if setfnamecorele ne -1 then begin                       ;version number found
        setfnamecore = strmid(setfnamecore0,0,setfnamecorele) ;extract final core name
        setfnameversion = int(strmid(setfnamecore0,setfnamecorele+5))
     endif else begin           ;no version nr = initial setup-file
                                ;      fitsetnamext = parasign('fitsetnamext','')
;fitsetnameext=''   ;turn this off again  
        setfnameversion = 0     
        setfnamecore = setfnamecore0 ;add eventual extension to name
     endelse
  endif                                                 ;setfnameflag
  override_version_nr= parasign('version_nr',0,/nowarn) ;override version nr from setup
  if override_version_nr gt 0 then setfnameversion=override_version_nr
  if isa(disflag) then dflag=disflag else dflag=3 ;default display-flag. (dflag cannot be defined in UFIT cmd-line due to being prior to common block)

;select fit-function
  model=strlowcase(parasign('model','utm',/nowarn)) ;function or program that generates model
  case model of
     'utm' : fitfunct='chisq_utm' 
;     'ufit_rvcurve' : fitfunct='chisq_rv' 
     else:  fitfunct='chisq_model' ;use this for generic modellers
  endcase
;  fitfunct='chisq_utm'          ;could later be replaced by read-in from setupfile
  priorfunct='chisqallprior'    ;could evtl. be replaced by read-in from setupfile

; get input-lc, using one of several keywords and assign data-vectors
  if paratest('indatafile') then indatafile=parasign('indatafile')  else $
     if paratest('inlcname') then indatafile=parasign('inlcname') else $
        if paratest('infilename') then indatafile=parasign('infilename') else $
           print,'WARN: need to give input lightcurve' ;compatibility with old param-name'

  yerrflag = fix(parasign('yerrflag',0,/nowarn)) ;indicates presence of column with y-errors in indatafile
;define up to 3 columns in indatafile with auxiliary data
  if yerrflag then ncols=3 else ncols=2 ;number of columns to read in indatafile

;if paratest('auxcol') then begin
  auxcol1=fix(parasign('auxcol',0,/nowarn)) ;note that auxcol counts from 1; if 0 then not used
  if auxcol1 ne 0 then ncols=ncols>auxcol1
  auxcol2=fix(parasign('auxcol2',0,/nowarn)) ;note that auxcol counts from 1; if 0 then not used
  if auxcol2 ne 0 then ncols=ncols>auxcol2
  auxcol3=fix(parasign('auxcol3',0,/nowarn)) ;note that auxcol counts from 1; if 0 then not used
  if auxcol3 ne 0 then ncols=ncols>auxcol3

  delchar=parasign('delchar'," ",/nowarn)    ;delimiter character in indatafile, default is 'space'
  comchar=parasign('comchar',"#",/nowarn)    ;comment character in indatafile, default is #
  firstrow=fix(parasign('firstrow',1,/nowarn)) ;first row to read from indatafile
  maxdatrow=fix(parasign('maxdatrow',0,/nowarn)) ;max number of rows to read from indatafile. Reads all rows if 0

;  lcarr=rdnumtab(indatafile,ncols)
  lcarr='dummy'
  rdtab,indatafile,datarr=lcarr,maxcol=ncols,comc=comchar,delstr=delchar,firstrow=firstrow,maxdatrow=maxdatrow

  xdat=double(lcarr[*,0])       ;x,y arrays with data to be fitted
  ydat=double(lcarr[*,1])
  if yerrflag then yerr = abs(double(lcarr[*,2])) else yerr = float(parasign('meas_error',1)) ;mesurement error for estimate of chi-square
;read up to 3 columns from indatafile with auxiliary data
  if auxcol1 ge 1 then auxdat1=lcarr[*,auxcol1-1]
  if auxcol2 ge 1 then auxdat2=lcarr[*,auxcol2-1]
  if auxcol3 ge 1 then auxdat3=lcarr[*,auxcol3-1]


  if paratest('inlcdigi') then begin ;opt removal of leading digits in timestamps
     ndigi=fix(parasign('inlcdigi')) ;number of precomma digits to keep
     xdat=digirem(xdat,ndigi) 
  endif

  fitalgo=parasign("fitalgo","am",/nowarn) ; type of fit algorithm
  
;---- 
;make up vector of initial fit parameters from setup file
;first, get NAMES of fit parameters from setup
  fitparname=strarr(16)         ;max of 16 fit parmaeters
  adim=0
  for i=0,n_elements(pnarr)-1 do begin
     if pnarr[i] eq 'fit' then begin ;search 'fit' keywords
                                ;   if adim eq 16 then print,"Ufit: Too many fit parameters specified!"
        fitparname[adim] = pvarr[i]  
        adim = adim+1           ;number of params
     endif
  endfor
  if adim eq 0 then print,"Ufit: At least 1 fit parameter needs to be specified!"
  fitparname = fitparname[0:adim-1]
;now get VALUES of fit params from setup
  a=dblarr(adim)
  for i=0,adim-1 do begin
     a[i]=double(parasign(fitparname[i],warn="UFIT: ERROR: value for fit parameter '"+fitparname[i]+"' not found in setupfile"))
  endfor
  aorg=a                        ;keep copy of original values
;  if dflag ge 1 then begin print,'ufit:  input start parameters--------------'
;     for i=0,adim-1 do print,fitparname[i]," :",a[i]
;  endif



;extract scale vals from setup
  scaleyes=paratest('scale')
  scale=fltarr(adim)            ;read in scale values
  if scaleyes then begin
     j = 0                      ;index in scale array
     for i=0,n_elements(pnarr)-1 do begin
        if pnarr[i] eq 'scale' then begin ;search for 'scale' keywords
       if j eq adim then message,"More scale parameters than fit parameters are given in setup. Needs to be the same."
           scale[j] = pvarr[i]  
           j += 1   
        endif
     endfor
     if j ne adim then message,"Fewer scale parameters than fit parameters are given in setup. Needs to be the same"
  endif else begin
     scale+=1
     if fitalgo ne 'mc' then print,'WARN: no scale parameters given, using scale=1 for all fit paramters'
  endelse

;define objects for priors
  pri={genprior,prinum:0,prvar:'name of variable',prvid:0}     ;generic prior (= no prior)
  gpri={gprior,prctr:0D,prsign:0D,prsigp:0D,inherits genprior} ;gauss prior
  apri={aprior,inherits gprior}                                ;absulute difference prior
  upri={uprior,prloli:0D,prhili:0D,inherits genprior}          ;uninformative or uniform prior. Initialize with NaN for 1-sided priors

;get number of priors from setup file and collect stuff to initialize them
  prtype=[]                     ;array with types of priors in pnarr, pvarr   
  prpnid=[]                     ;id of prtype statement in pnarr
  for i=0,n_elements(pnarr)-1 do begin
     if pnarr[i] eq 'prtype' then begin ;prior variable name
        prtype = [prtype,pvarr[i]]      ;add type of prior to array
        if pvarr[i] eq 'nprior' then prtype[-1] = 'genprior' ;rename 'nprior' (no prior) 
        prpnid = [prpnid,i] 
     endif
  endfor
  nprior=n_elements(prtype)                     ;number of priors
  if nprior ge 1 then prflag = 1 else prflag =0 ;flag on existence of priors
  if prflag then begin                          ;priors are in setup and will be defined
     prparmax= intarr(nprior)+3                 ; 3 is max number of prior-params (beyond prtype, prvar)
     prparmax[-1] = 3 < (n_elements(pnarr) - prpnid[-1] -2) ;for last prtype statement, prparmax is at most until end of pnarr
     prior=objarr(nprior)
     for i = 0, nprior-1 do begin ;prior objects of type prtype are created and initialized
        prior[i]=obj_new(prtype[i],i,prpnid[i],prparmax[i])
     endfor
  endif                         ; prflag eq 1
;assign values to fit paramenters
  if adim eq 0 then print,"Ufit: At least 1 fit parameter needs to be specified!"
  fitparname = fitparname[0:adim-1]
;now get VALUES of fit params from setup
  a=dblarr(adim)
  for i=0,adim-1 do begin
     a[i]=double(parasign(fitparname[i]))
  endfor
  aorg=a                        ;keep copy of original values
;  if dflag ge 1 then begin print,'ufit:  input start parameters--------------'
;     for i=0,adim-1 do print,fitparname[i]," :",a[i]
;  endif

                                ;randomize initial values
  randominit= float(parasign('randominit',0,/nowarn)) ;multiplier for randomizing input values
  if randominit ne 0. then begin
     if fix(parasign('random_uniform',0,/nowarn)) eq 1 then uniform=1 else uniform =0
;uniform dist or gauss
;create unique name for fitsetmp that is based on systime to 1/100s sec
     stime=frac(systime(1)/10000)*10000
     seed=stime                 ;this, since randomu will modify seed
     if uniform then begin
        offset=(2*randomn(seed,adim,/uniform)-1.) *randominit*scale ;unif. dist from -scale to +scale
     endif else offset=randomn(seed,adim)*randominit*scale          ;gauss dist with sigma=scale
     a+=offset                                                      ;add to a
     for i=0,adim-1 do begin
        paraput,fitparname[i],a[i]
                                ;print,aorg[i],scale[i],offset[i],a[i]
     endfor
     paradd,'randominit',0      ;set this to zero for output setup and the intermedite setup in next lines
     if override_version_nr ge 1 then begin
        setfname_actual=setfnamecore+'.f_in'+str(setfnameversion)+'.utm' ; setup-file with actaully used intial values
        savesetupfile,setfname_actual
     endif                      ;override_version_nr ge 1
  endif                         ;randomize


;initialize remaining common-block vars
;adim=n_elements(atrue)
  ndat=n_elements(xdat)         ;number of data-points to be fitted
  ymod=dblarr(ndat)             ;will keep fitted y values
  oobflag = 0
;  oobchisq=float(parasign('oobchisq',1e6,/nowarn))  ;default chisq for out-of bounds condition
;
;-------
;override some param values in the set-up arrays
  tflag_orig=fix(parasign('tflag',1,/nowarn))       ;keep original flag setting
  paraput,'tflag','1' ,/nowarn                      ;use only time-points (x-vals in utm)

  plcflag_orig= fix(parasign('plcflag',99,/nowarn)) ; ;keep original flag setting (if absent, give 99)
  paradd,'plcflag','-1'                             ;don't generate final lc display

  if paratest('onoise') then paraput,'onoise','0'     ;set onoise to zero if there
  if paratest('chisqadd') then paraput,'chisqadd','0' ;reset chisqadd to zero if there

  wflag_orig=fix(parasign('wflag',99,/nowarn)) ;Determines if last fitted model-lc is being saved. Keep original, absent=99
  paradd,'wflag',0                             ;don't save intermediate lc's 
  

;create unique name for fitsetmp that is based on systime to 1/100s sec
  stime=frac(systime(1)/10000)*10000
  seed=stime                    ;this, since randomu will modify seed
  ranu=randomu(seed)*10000      ;random number

  fitsetmp=strcompress(string(format='(f15.2,f7.1)',stime,ranu),/rem)+'ufit.tmp'
;print,fitsetmp
  savesetupfile,fitsetmp


;prepare plots

  if keyword_set(nofit) and dflag lt 3 then dflag =3 ;to force plotting

  residlevelflag = paratest('residlevel')                                              ;check if residlevel is defined
  if residlevelflag then residlevel = float(parasign('residlevel')) else residlevel=0. ;level at which to display residuals

  tmp= fix(parasign("showfinalfitflag",'0',/nowarn))                                ;show separate plot of final fit 
  if tmp eq 0 then  if dflag ge 2 then showfinalfitflag=1 else showfinalfitflag =0  ;uses defaults from dflag
  if tmp ne 0 then showfinalfitflag = (1+tmp)/2                                     ;map flag values -1,1 to 0,1

  tmp= fix(parasign("errplotflag",'0',/nowarn))                                ;plot error-bars for data 
  if tmp eq 0 then  if yerrflag then errplotflag=1 else errplotflag =0         ;uses default (from yerrflag)
  if tmp ne 0 then errplotflag = (1+tmp)/2                                     ;map flag values -1,1 to 0,1

  if dflag ge 2 or showfinalfitflag eq 1 then begin ;a plot is made at some moment; do preps 
     if paratest('xunit') then xtitle=parasign('xunit')  else $
        if paratest('tunit') then xtitle=parasign('tunit') else $ ;compatibility w UTM convention of tunit
           xtitle=''


     oflag=fix(parasign('oflag',0,/nowarn)) ;type of y-axis (reversed for 2=delta-mags)
     case oflag of
        0: if paratest('yunit') then ytitle=parasign('yunit')  else $
           if paratest('lunit') then xtitle=parasign('lunit') else $ ;compatibility UTM convention of lunit
              ytitle=''
        1: ytitle='$\Delta$ F/F$_0$'
        2: ytitle='$\Delta$ mag'
     endcase
  endif

  if dflag ge 3  then begin     ;make initial plot     
     modeloffvalue=0.           ;need this for first plotlcymod invoke
     plotlcymod,/first,xtitle=xtitle,ytitle=ytitle,errplot=errplotflag ;Plot the initial data
  endif
  if dflag ge 1  then begin
     print,'running utm for stats on initial setup'
     paraput,'tflag',1,/nowarn                             ;only generate model
     if keyword_set(nofit) then paraput,'wflag',wflag_orig ; generate output-lc if wanted
     savesetupfile,fitsetmp
                                ;display initial params  
     print,'------------ufit: initial parameters--------------'
     for i=0,adim-1 do print,fitparname[i]," :",a[i]
     if dflag ge 2 then doplot=1 else doplot=0
     chisqmod=call_function(fitfunct,/sigall,/sigontr,/sigofftr,/print,chisqadd=chisqadd, doplot=doplot) ;print initial stats of model-fit
     chisq=call_function(priorfunct,chisqadd=chisqadd, chisqmod=chisqmod, /print)                        ;adds chisq of priors
  endif
;  if dflag ge 3 then plotlcymod,/update         ;display initial model
  if keyword_set(nofit) then  plotlcymod, /last ;do 'final' plot already now
  if keyword_set(nofit) then goto,cleanup
  case fitalgo of
     
;------------------------Amoeba fit--------------------------------
     'am' : begin    
        if dflag ge 2 then print,'Amoeba fitting'
        
        ftol=float(parasign('amoeba_ftol',1e-5,/nowarn)) ;ftol param see help amoeba. New param name; checks for old one only if isn't given
        if ftol eq 1e-5 then ftol=float(parasign('ftol',1e-5,/nowarn)) ;old ftol param
        
        func_vals=dblarr(adim+1) ; array for output params
        ncalls=0                 ;number calls to amoeba (output param)
        nmax=uint(parasign('amoeba_nmax', 65535,/nowarn)) ; max number calls to amoeba
        itct=0L
        aout = exofast_amoeba(ftol,FUNCTION_NAME = 'amoeba_utm',function_val=func_val,p0=a,scale=scale,nmax=nmax,ncalls=ncalls)
                                ;update 'aout' into parameterfile
        tmp=size(aout)  
        if tmp[0] eq 0 then begin ; Max number of AMOEBA function calls exceeded
                                ;in this case the last fitsetmp internal to amoeba is taken as output
           print,"Max number of AMOEBA function calls exceeded (nmax paramter)"
           paradd,'amoeba_nmax_exceed',1
        endif else begin
           
           for i=0,adim-1 do begin
              paraput,fitparname[i],aout[i]
           endfor
           paradd,'amoeba_nmax_exceed',0 ;put this in case an input setup with a previous exceed=1 was used
        endelse
     if dflag ge 1 then begin print & print,'AMOEBA finished. Number of calls to fit-fct: ',ncalls & end
        paradd,'amoeba_ncalls_done',ncalls
     end

;------------------------------MCMC run---------------------------------
     "mc" : begin
        if dflag ge 1 then print,'Executing EXOFAST_DEMC'
        outmcname=setfnamecore+'.mc'+str(setfnameversion)+'.sav' ;name MCMC .sav file         
        overwrite=0
        while file_test(outmcname) and overwrite eq 0 do begin ;mc output already exists
           warnstr= 'WARNING: MCMC output file '+outmcname+' already exists. Overwrite Y/N: '
           tmps=' ' & read,prompt=warnstr,tmps
           if strupcase(tmps) eq 'Y' then overwrite =1 else begin ;ask for new version nr
              newvers=0
              read,prompt='Give other version number (current: '+str(setfnameversion)+' ) : ',newvers
              setfnameversion=newvers
              outmcname=setfnamecore+'.mc'+str(setfnameversion)+'.sav' ;new MCMC .sav file   
           endelse
        endwhile
        itct=0L
        mc_out_pars=1.
        mc_out_chisq =1.
;check for presence of parameters that may override defaults in exofast_demc
        mc_maxsteps = long(parasign('mc_maxsteps',50000)) ;overrides exofast_demc defaults
        mc_nchains = long(parasign('mc_nchains',5))
        mc_getscale=int(parasign('mc_getscale',1,/nowarn))
        mc_stopscale=int(parasign('mc_stopscale',0,/nowarn)) ;if 1, UFIT ends after determining scales
        mc_initmult=float(parasign('mc_initmult',2,/nowarn)) ; initial parameter scatter
        if mc_getscale then begin
           print,'MCMC stepping scales to be determined automatically' 
           mc_dcmult=float(parasign('mc_dcmult',5,/nowarn)) ;multiplier for get_mcmcscale
        endif else begin 
           print,'MCMC stepping scales taken from setup file'
           mc_scale=scale       ;use those of setupfile
           mc_stopscale=0
        endelse
        if not(mc_stopscale) then print,'MCMC chain will be written to ',outmcname     
        ufit_exofast_demc,a,'mcmc_utm',mc_out_pars,chi2=mc_out_chisq,nchains=mc_nchains,maxsteps=mc_maxsteps,burnndx=burnndx,scale=mc_scale,dcmult=mc_dcmult,initmult=mc_initmult,stopscale=mc_stopscale,autoseedscale=1

        if mc_getscale then begin  ;new scales were determined; update in setup array           
           if scaleyes  then begin ;scale params are in setup
              j=0
              for i=0,n_elements(pnarr)-1 do begin
                 if pnarr[i] eq 'scale' then begin ;search for 'scale' keywords
                    pvarr[i] = mc_scale[j]    
                    j += 1   
                 endif
              endfor
           endif else begin     ;no scale params in setup; add them to setup
              for j=0, adim -1 do begin
                 pnarr=[pnarr,'scale']
                 pvarr=[pvarr,string(mc_scale[j])]
              endfor
           endelse
        endif                   ;if mc_getscale

        if (mc_stopscale) then begin ;no MCMC chain was run but update fit-params by better ones that might have been found during scale determination
           aout=a       
        endif else begin                ;a full MCMC chain was run
           paradd,'mc_burnndx',burnndx  ;add burn-index to output setups


           save,mc_out_pars,mc_out_chisq,burnndx,fitparname,filename=outmcname 
           print,"MCMC output: analyze it with: ufit_mcshow,'"+outmcname+"'"
           print,"MCMC 1-sigma model lims, use: ufit_limitshow,'"+setfname+"','"+outmcname+"'"

; extract useful part of MC chains
           lastid=n_elements(where(MC_OUT_CHISQ[*,0] gt 0))-1 ;id of last filled element
           MC_OUT_PARS = MC_OUT_PARS[*,burnndx:lastid,*]      ;remove burn-in part, unused one
           MC_OUT_CHISQ = MC_OUT_CHISQ[burnndx:lastid,*]      ;remove burn-in part
           nstep=lastid-burnndx+1                             ;number of saved steps

;do some processing to find best chi-sq 
           j_pars = reform(MC_OUT_PARS, adim, long(nstep*mc_nchains)) ;join all chains
           j_chisq = reform(MC_OUT_CHISQ , long(nstep*mc_nchains))
           minchi=min(j_chisq,idmin)
;           print,'Best fit: chisq= ',minchi
           aout=j_pars[*,idmin]   ;best-fit parameters          
        endelse                   ; if(mc_stopscale)

        for i=0,adim-1 do begin
           paraput,fitparname[i],aout[i]
        endfor
     end                        ;MCMC fit
     
     else: print,"invalid fit-algorithm specified in setup parameter fitalgo"
  endcase

;rerun UTM a last time with best fit params due to one or more reasons:
;- last invoke of UTM by Amoeba often not with final parameters, hence y-model data not correct
;-  to save final model-lc
;-  to get the intrarr
  if dflag ge 3 then  print,'-----------------------------------------------------------'
  if dflag ge 3 then print,'running utm once more with best-fit parameters'
  if wflag_orig eq 1 then begin ;generate output model-lc name
     if setfnameflag then  begin
        outlcname=setfnamecore+'.fout'+str(setfnameversion)+'.lc' 
     endif else outlcname=setfname+'.lc'
     paradd,'outlcname',outlcname      ;add or overwrite outlcname to setup
     paraput,'wflag',wflag_orig        ; generate output-lc if wanted
     if dflag ge 1 then print,'saving model light-curve to file: ',outlcname
  endif
  paradd,'plcflag',-1                   ;don't generate final lc display
  paraput,'tflag',1,/nowarn             ;only generate model
  savesetupfile,fitsetmp
  intrarr =0                       ;array needed for calc of  off-transit stddev
  rms_res=1.                       ;to initalize keyword
  chisqmod=call_function(fitfunct,sigall=rms_res,chisqadd=chisqadd)
  chisq=call_function(priorfunct,chisqadd=chisqadd, chisqmod=chisqmod) ;chisq of priors
  paradd,'chisq',chisq                                                 ;add chi-sq to output setups
  paradd,'rms_res',rms_res                                             ;add rms of resids to output setups


;save setup with fitted params
  paraput,'tflag',tflag_orig,/nowarn ;return flag to initial val 
  if wflag_orig eq 99 then pardel,'wflag'
  if plcflag_orig ne 99 then  paraput,'plcflag',plcflag_orig ;return flag to initial val if it was there
  if setfnameflag then  begin
     setfnameout=setfnamecore+'.fout'+str(setfnameversion)+'.utm' 
     savesetupfile,setfnameout
     savenextfinflag=int(parasign('savenextf_in',1,/nowarn))
     if savenextfinflag ge 1 then if override_version_nr ne 0 then begin
        savenextfinflag = 0
        print,'WARN: an output version_nr was specified, flag savenextf_in set to 0'
     endif
     if savenextfinflag then begin
        setfname_nextin=setfnamecore+'.f_in'+str(setfnameversion+1)+'.utm' ;next input file; is copy of .fout but with next version nr. Disabled if override version nr has been set.
        savesetupfile,setfname_nextin
     endif
  endif else begin
     setfnameout=setfname+'.fit'
     savenextfinflag = 0
  endelse
  savesetupfile,setfnameout
  if dflag ge 1 then begin
     print,'saved parameters to setup: ',setfnameout
     if savenextfinflag then  print,'duplicate setup for next fit: ',setfname_nextin
  endif

;display fitted params and final plot 

  if dflag ge 1 or showfinalfitflag then begin
;  if dflag ge 1 and not(mc_getscale) then begin
     print,'------------ufit: best-fit parameters--------------'
     if aout[0] ne -1 then for i=0,adim-1 do print,fitparname[i]," :",aout[i]
                                ; print,'number iterations: ',iterct
     chisqmod=call_function(fitfunct,/sigall,/sigontr,/sigofftr,/print,chisqadd=chisqadd,doplot=1) ;print final stats
     chisq=call_function(priorfunct,chisqadd=chisqadd, chisqmod=chisqmod, /print)                  ;chisq of priors
  endif
  if showfinalfitflag then begin
     if ~isa(modeloffvalue) then modeloffvalue=0.
     plotlcymod,xtitle=xtitle,ytitle=ytitle,/first,/last,errplot=errplotflag ;do a clean final plot
  endif


;remove temp fitsetfile
  cleanup:  keepfitempset_flag = fix(parasign('keepfitempset',0,/nowarn))
  if not(keepfitempset_flag) then begin ;delete temporary setup files
     exestr='\rm '+fitsetmp       
     spawn,exestr
  endif
  
  if isa(prior) then obj_destroy,prior ;remove objects of prior class 
  prior=[]                             ;to also destroy the empty heapvar
  if dflag ge 1 then print,'------------------------ufit finished-------------------------'
  if dflag ge 6 then stop

end





