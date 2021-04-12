@setupfile

pro ufit_limitshow, setupfile,mcsavfile,nshow=nshow,nxmod=nxmod,plotxdat=plotxdat,showreject=showreject,chain=chain,sigma=sigma ;,contbin=contbin,burn=burn,contour=contour,nosave=nosave,chimax=chimax,nitermax=nitermax,_extra=ex 

;read mc chain and  corresponding setupfile and
;generates plot with models from mcmc iteration that are within 1 sigma of the best-fit one. Also plots the best-fit case.

;keyword: 
;  nshow:  gives or limits the number of overplotted models. Is done by displaying only every n'th one so that not more than nshow are displayed. Defaults are 100 plots. nshow=-1 forces display of all.
; nxmod:  gives number of x-points for which model is calculated. Default is 50 equidistant pts between min and max of x-data values
; plotxdat : instead of calculating & plotting the model at the nxmod data points, plots it at the input x-values
; showreject :  rejected (non-accepted) iteration steps are also included in extracted list of parameters 
; chain : vector of chain-numbers to be analyzed. Permits removal of stale chains.
; sigma  : maximum delta-chisq against the minimum chisq for which models are plotted. Default is 1 (= 1 sigma)
  
;HJD 28jul2020, mostly glued together from ufit_mcshow and from ufit itself.
;     6feb2021 fixed bug in labelling of plot axes
;    12apr2021 forced all utm plotting flags to be turned off

;example use: IDL> ufit_limitshow,'polyfit.fout1.utm','polyfit.mc1.sav',nshow=50,sigma=2


  readsetupfile,setupfile

  restore, mcsavfile
  rootname=strmid(mcsavfile,0,strlen(mcsavfile)-4) ;extract core name

  lastid=n_elements(where(MC_OUT_CHISQ[*,0] gt 0))-1-1 ;id of last filled element. Discard last step, since it might be incomplete (if generated from stop during demc.pro execution)

  if keyword_set(nitermax) then lastid = lastid < nitermax

  print,'number of executed steps',lastid
  if isa(burnndx) then   print,'burnndx: ',burnndx else begin
     print,"burnndx not yet defined"
     burnndx =0
  endelse

  if isa(chain) then begin      ;analyze only the chains given in kw
     MC_OUT_PARS = MC_OUT_PARS[*,*,chain]
     MC_OUT_CHISQ = MC_OUT_CHISQ[*,chain]

;if n_elements(chain) eq 1 then chain=reform(chain,1)  ;convert number to array  

  endif
  if keyword_set(burn) then begin             ;remove unused last part of chain
     MC_OUT_PARS = MC_OUT_PARS[*,0:lastid,*]  ;array of [npar,niter,nchain]
     MC_OUT_CHISQ = MC_OUT_CHISQ[0:lastid,*]
  endif else begin              ;remove burn-in and unused last part
     if lastid le burnndx then print,'nitermax needs to be bigger than burnndx'
     MC_OUT_PARS = MC_OUT_PARS[*,burnndx:lastid,*] 
     MC_OUT_CHISQ = MC_OUT_CHISQ[burnndx:lastid,*] 
  endelse
  tmp= size(MC_OUT_PARS,/dimensions)
  npar=tmp[0]
  nstep=long(tmp[1])
  if n_elements(tmp) eq 3 then nchain=tmp[2] else nchain=1

  niter=nchain*nstep            ; total number of MCMC iterations used
  print,'number of chains: ',nchain
  print,'length of chains: ',nstep
  print,'total number of iterations: ',niter

  if total(size(fitparname)) eq 0 then fitparname=sindgen(npar)+'par' ;create param-names if not in .sav

;  print,'i: evoluton of parameters in MCMC chains'
;  parv=parevol_plot(fitparname,mc_out_pars,statvals=mc_out_chisq,chimax=chimax)
  
  if ~isa(chain) then chain=indgen(nchain)  

; keep only chain-elements that correspond to accepted changes of parameters
; replace all unaccepted chain elements by NaNs
  nans = make_array(npar,/float,value=!values.f_nan)
  if ~keyword_set(showreject) then begin
     for i = 0,nchain-1 do begin
        nacc=nstep                  ;counter of accepted links
        for j=nstep-1,1,-1 do begin ;go from back of chain to keep links that were accepted
           if  array_equal(MC_OUT_PARS[*,j-1,i], MC_OUT_PARS[*,j,i]) then begin
              MC_OUT_PARS[*,j,i] = nans
              MC_OUT_CHISQ[j,i]=!values.f_nan
              nacc -= 1
           endif
        endfor
        print,'chain ',chain[i],' number of accepted iterations: ',nacc,' =',float(nacc)/nstep*100,'%'
     endfor
  endif

  if nchain eq 1 then begin  
     MC_OUT_PARS = reform(MC_OUT_PARS,npar,nstep,1) ;change to a 3D array
     MC_OUT_CHISQ = reform(MC_OUT_CHISQ,nstep,1)    ;change to 2D array
  endif

  j_pars = reform(transpose(MC_OUT_PARS,[0,2,1]), npar, niter) ;join all chains and put dims in order par, step, chain
  j_chisq = reform(transpose(MC_OUT_CHISQ) , niter)            ; order step, chain
  chainid= lindgen(niter) mod nchain                           ; vector with ids of chains


  minchi=min(j_chisq,idmin)
  medianpars=0.
  bestpars=j_pars[*,idmin]      ; set of best parameters
  print,'best chisq: ',minchi, ' at step (incl burnndx):',long(idmin/nchain)+burnndx

  if keyword_set(sigma) then chimax=minchi+sigma else chimax=minchi+1. ;keep only pars below chimin+1 
;put upper limit to chi-sq
;  if keyword_set(chimax) then begin
  goodchi=where(j_chisq le chimax)
  j_chisq=j_chisq[goodchi]
  j_pars=j_pars[*,goodchi]
  chainid=chainid[goodchi]
  print,'number of steps below chimax: ',n_elements(goodchi)
                                ; endif

                                ;remove NaNs
  if ~keyword_set(showreject) then begin
     finiteid=where(finite(j_chisq))
     naccept=n_elements(finiteid)
     j_chisq=j_chisq[finiteid]
     j_pars=j_pars[*,finiteid]
     chainid=chainid[finiteid]
                                ;    print,'total number of accepted iterations: ',naccept,' =',float(naccept)/niter*100,'%'
  endif


  if not(keyword_set(nshow)) then nshow=100L                         ; max number of pts that are plotted
  if nshow eq -1 or naccept le nshow then stride =1 else begin       ;plot all pts
     stride = ((naccept-1)/nshow) +1                                 ;stride in plotting the datarr
     print,"Sample of ",strcompress(string(naccept))," points. Every ",strcompress(string(stride)),"'th plotted"
  endelse




  
                                ; print,'ii: scatter plots of params against stdev'
                                ; psdev=plot_par1_vs_params(fitparname,j_pars,statname='chisq',statvals=j_chisq,yminor=0,xminor=0,symbol='dot', nshow=nshow, _extra=ex)

;prep execution of modeller and plotting of data
  model=strlowcase(parasign('model','utm')) ;function or program that generates model

; get input-lc, using one of several keywords and assign data-vectors
  if paratest('indatafile') then indatafile=parasign('indatafile')  else $
     if paratest('inlcname') then indatafile=parasign('inlcname') else $
        if paratest('infilename') then indatafile=parasign('infilename') else $
           print,'WARN: need to give input lightcurve' ;compatibility with old param-name'

  yerrflag = fix(parasign('yerrflag',0,/nowarn)) ;indicates presence of column with y-errors in indatafile
  if yerrflag then   lcarr=rdnumtab(indatafile,3) else   lcarr=rdnumtab(indatafile,2)
  xdat=lcarr[*,0]               ;x,y arrays with data to be fitted
  ydat=lcarr[*,1]
  if yerrflag then yerr = abs(lcarr[*,2]) else yerr = float(parasign('meas_error',1)) ;mesurement error for estimate of chi-square
  if paratest('inlcdigi') then begin                                                  ;opt removal of leading digits in timestamps
     ndigi=fix(parasign('inlcdigi'))                                                  ;number of precomma digits to keep
     xdat=digirem(xdat,ndigi) 
  endif

  tmp= fix(parasign("errplotflag",'0',/nowarn))                                ;plot error-bars for data 
  if tmp eq 0 then  if yerrflag then errplotflag=1 else errplotflag =0         ;uses default (from yerrflag)
  if tmp ne 0 then errplotflag = (1+tmp)/2                                     ;map flag values -1,1 to 0,1

  if paratest('tflag') then paraput,'tflag','1'  ;reset to 1 for UTM

;create unique name for fitsetmp that is based on systime to 1/100s sec
  stime=frac(systime(1)/10000)*10000
  seed=stime                    ;this, since randomu will modify seed
  ranu=randomu(seed)*10000      ;random number
  fitsetmp=strcompress(string(format='(f15.2,f7.1)',stime,ranu),/rem)+'ufit.tmp'
;print,fitsetmp



;do first plot of data
  if paratest('xunit') then xtitle=parasign('xunit')  else $
     if paratest('tunit') then xtitle=parasign('tunit') else $ ;compatibility w UTM convention of tunit
        xtitle=''
  if paratest('yunit') then ytitle=parasign('yunit')  else $
     if paratest('lunit') then ytitle=parasign('lunit') else $ ;compatibility UTM convention of lunit
        ytitle=''

  if ~errplotflag then begin
     ymin=min([ydat,ydat])
     ymax=max([ydat,ydat])

;     plc=  plot(xdat,ydat,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax])    ;input data as crosses
  endif else begin                                                                                 ;do erroplot
     if n_elements(yerr) eq 1 then yerr=replicate(yerr,ndat)                                       ;set to same error
     ymin=min([ydat,ydat-yerr])
     ymax=max([ydat,ydat+yerr])
;     plc=errorplot(xdat,ydat,yerr,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax])
  endelse

;extend yrange by 3% above, below
  ymin=ymin-(ymax-ymin)/33.
  ymax=ymax+(ymax-ymin)/33.


  if ~keyword_set(plotxdat) then begin ;plot model at the nxmod points in x
     if ~isa(nxmod) then nxmod=50      ;number of model pts
     xmod=dindgen(nxmod)/(nxmod-1)*(max(xdat)-min(xdat))+min(xdat)
  endif else xmod=xdat          ;plot model at input x-vales

  print,'Plotting data-set number: '
  for i=0, naccept-1,stride do begin
     for j=0,npar-1 do begin    ; assign paramters
        paraput,fitparname[j],j_pars[j,i]
     endfor

if model eq 'utm' then begin ;supress any plotting
paradd,'plcflag',-1 
paradd,'plcsimflag',-1 
paradd,'plotorbflag',-1 
paradd,'dispstarflag',-1
paradd,'xpflag',0 
endif
     savesetupfile,fitsetmp
     call_procedure,model,fitsetmp,0,xmod,ymod
     if i eq 0 then over =0 else over=1
     plc=plot(xmod,ymod,'-',color='silver',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax],thick=2,over=over) ;overplotting of model-fcts
;     print,'i=',i
     print,format='($,i0," ")',i
  endfor
  print,' '
                                ;stop

;overplot best fit
  for j=0,npar-1 do begin       ; assign paramters
     paraput,fitparname[j],bestpars[j]
  endfor
  savesetupfile,fitsetmp
  call_procedure,model,fitsetmp,0,xmod,ymod
  plc=plot(xmod,ymod,'g-',thick=2,/over)

; overplot data
  if ~errplotflag then begin
;     ymin=min([ydat,ydat])
;     ymax=max([ydat,ydat])
     plc=  plot(xdat,ydat,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,/over)                 ;input data as crosses
  endif else begin                                                                                 ;do erroplot
;     if n_elements(yerr) eq 1 then yerr=replicate(yerr,ndat)                                       ;set to same error
;     ymin=min([ydat,ydat-yerr])
;     ymax=max([ydat,ydat+yerr])
     plc=errorplot(xdat,ydat,yerr,'+',xtitle=xtitle,ytitle=ytitle,ystyle=2,xstyle=3,yrange=[ymin,ymax],/over)
  endelse
  


;saving of plots
  if not(keyword_set(nosave)) then begin
     outrootn=strmid(mcsavfile,0,strlen(mcsavfile)-4)+'.png'
     outfile='sigma_model_'+outrootn
     plc.save,outfile
     print,'saved limit-plot to file ', outfile 

  endif                         ;nosave
end


