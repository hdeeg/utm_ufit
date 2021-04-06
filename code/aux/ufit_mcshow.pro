pro ufit_mcshow, restorefile,contbin=contbin,burn=burn,contour=contour,nosave=nosave,chimax=chimax,nshow=nshow,nitermax=nitermax,chain=chain,nochain=nochain,showchaincol=showchaincol,showreject=showreject,_extra=ex 

;generates various plots from mcmc iteration made by ufit_exofast_demc
; restorefile  : .sav file from exofast-demc or ufit_exofast-demc 
;
; Keywords:
; nbatch  : size of step-increments for 'running plot' showing advances 
; burn  : if set, the burnndx value from the .sav file is overriden and the 
;         chain starting at burn is analyzed. Use -1 to show full chain.
; contbin : number of bins in X and Y for contours in correlation plots
; contour  : if set, correlation plot with 1,2-sigma contours instead of scatterplot
; chimax  : consider only points with a chisq below or equal to chimax
;  nosave  :if set, plots are not saved to png files
;  nshow: parevol_plot, chisq-scatterplot and correlation scatterplot limit the number of points, by displaying only every n'th point so that not more than nshow are displayed. Defaults are 1000, 5000 and 10000 pts. nshow=-1 forces display of all points
; nitermax : the data shown will be truncated to nitermax iterations.
; chain : vector of chain-numbers to be analyzed. Permits removal of stale chains.
; nochain : vector of chain-numbers to be supressed. Permits removal of stale chains. Is inverse of chain keyword and overrides it.

; showchaincol : if set, shows additional window indicating colors of chains
; showreject :  rejected iteration steps are also includeed in histogram of parameters and in calc of parameter mean.

;
;written by Hans J Deeg
;
;version history
; 28feb2019 first dated version
; 10mar2019 added plot of params against chisq; added nosave
; 13mar 2019 added _extra keyword
; 25apr2019 removed obsolete parid cmd-line param
; 27apr2019 added burn keyword, moved runp routine to end
; 29apr2019  added chimax
; 3may19 added parevol
; 9may19 routine to remove 'unaccepted' chain-elements that are equal to previous one. Before, frozen chains affected stronlgy the stats.  Also added showreject keyw 
; 11may19 propagated chimax to parevol_plot
; 12mar20 made scatterplot default, instead of contour. Removed the scatter kw 
;    and added contour kw
; 22mar20 added chain kw
;  7aug20 modified burn kw to permit input of starting index.
;         added showchaincol kw to show window with chain colors
;         added nochain kw
; 22mar21 added outcommented example about converting parameters from one into an other for display and stats, see around line 140


  restore, restorefile
  rootname=strmid(restorefile,0,strlen(restorefile)-4) ;extract core name

  lastid=n_elements(where(MC_OUT_CHISQ[*,0] gt 0))-1-1 ;id of last filled element. Discard last step, since it might be incomplete (if generated from stop during demc.pro execution)

                                ;get original Nr of chains of(MC_OUT_PAR
  tmp= size(MC_OUT_PARS,/dimensions)
  if n_elements(tmp) eq 3 then nchainorg=tmp[2] else nchainorg=1

  if keyword_set(nitermax) then lastid = lastid < nitermax

  print,'number of executed steps',lastid
  if keyword_set(burn) then burnndx=(burn>0) ;set burnndx to 0 if -1 given 
  if isa(burnndx) then   print,'burnndx: ',burnndx else begin
     print,"burnndx not defined in .sav file and is set to 0. Override w burn keyword"
     burnndx =0
  endelse

  if keyword_set(nochain) then begin ;generate inverse of chain keyword
     chainid=bytarr(nchainorg)+1
     chainid[nochain]=0
     chain=where(chainid)       ; verctor with indizes of valid chains
  endif

  if isa(chain) then begin      ;analyze only the chains given in kw
     MC_OUT_PARS = MC_OUT_PARS[*,*,chain]
     MC_OUT_CHISQ = MC_OUT_CHISQ[*,chain]
  endif else chain=indgen(nchainorg)
  
                                ;remove burn-in and unused last part
  if lastid le burnndx then print,'nitermax needs to be bigger than burnndx'
  MC_OUT_PARS = MC_OUT_PARS[*,burnndx:lastid,*] 
  MC_OUT_CHISQ = MC_OUT_CHISQ[burnndx:lastid,*] 

  tmp= size(MC_OUT_PARS,/dimensions)
  npar=tmp[0]
  nstep=long(tmp[1])
  if n_elements(tmp) eq 3 then nchain=tmp[2] else nchain=1

  niter=nchain*nstep            ; total number of MCMC iterations used
  print,'number of chains: ',nchain
  print,'length of chains: ',nstep
  print,'total number of iterations: ',niter

  if total(size(fitparname)) eq 0 then fitparname=sindgen(npar)+'par' ;create param-names if not in .sav


;plot additional window with colors of chain-nrs in parevolplot
  if keyword_set(showchaincol) then begin
     w=window(dimensions=[200,400],/no_toolbar)
     loadct,39,rgb_table=rgbtab,ncolor=nchain+1 ;use  nchain +1 to avoid use of the last color index (is white in table 39) 
     for i=0, nchain-1 do begin
        j=chain[i]
        rgbcol=transpose(rgbtab[i,*])
        p=plot([0,1],[j,j],'-',thick=2,sym_color=rgbcol,color=rgbcol,over=(i<1),ytickinterval=1,yminor=0,xshowtext=0,xmajor=0,xminor=0, title='chain-nr versus color',yrange=[-0.9,nchainorg-0.1],/current)
     endfor
  endif

  print,'i: evoluton of parameters in MCMC chains'
  parv=parevol_plot(fitparname,mc_out_pars,statvals=mc_out_chisq,chimax=chimax)

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
  j_chisq = reform(transpose(MC_OUT_CHISQ) , niter) ; order step, chain
  chainid= lindgen(niter) mod nchain                ; vector with ids of chains



  
;below example to convert MCMC paramters to others
;here going from MCMC output using q1, q2 (Kipping2013 uniform samnpling LD coeffc) to the usual u1, u2:
;print,fitparname[5:6]
;q1=j_pars[5,*]
;q2=j_pars[6,*]
;u1=2*sqrt(q1)*q2      
;u2=sqrt(q1)*(1-2*q2)                        
;j_pars[5,*]=u1
;j_pars[6,*]=u2
;fitparname[5]='u1'
;fitparname[6]='u2'


;put upper limit to chi-sq
  if keyword_set(chimax) then begin
     goodchi=where(j_chisq le chimax)
     j_chisq=j_chisq[goodchi]
     j_pars=j_pars[*,goodchi]
     chainid=chainid[goodchi]
     print,'number of steps below chimax: ',n_elements(goodchi)
  endif

                                ;remove NaNs
  if ~keyword_set(showreject) then begin
     finiteid=where(finite(j_chisq))
     naccept=n_elements(finiteid)
     j_chisq=j_chisq[finiteid]
     j_pars=j_pars[*,finiteid]
     chainid=chainid[finiteid]
     print,'total number of accepted iterations: ',naccept,' =',float(naccept)/niter*100,'%'
  endif

  minchi=min(j_chisq,idmin)
  medianpars=0.
  bestpars=j_pars[*,idmin]      ; set of best parameters
  print,'best chisq: ',minchi, ' at step (incl burnndx):',long(idmin/nchain)+burnndx

  
  print,'ii: scatter plots of params against stdev'
  psdev=plot_par1_vs_params(fitparname,j_pars,statname='chisq',statvals=j_chisq,yminor=0,xminor=0,symbol='dot', nshow=nshow, _extra=ex)

  print,'iii: parameter histograms'
; chisq distrib. plots and calc. of median values and 68% confidence intervals

  parhi= ufit_exofast_plotdist(j_pars,medianpars=medianpars,parnames=fitparname,bestpars=bestpars)
  for i=0,npar-1 do begin
     print,fitparname[i],': best:',bestpars[i],' distr.mean:',medianpars[0,i],' +',medianpars[1,i],' -',medianpars[2,i]
  endfor

  print,'iv: parameter correlation plots'
  if not(keyword_set(contbin)) then contbin=10
  if keyword_set(contour) then contour=1 else contour=0
  sm= scattermatrix(j_pars,fitparname, dimension=[600,600], bestpars=bestpars,contour=contour, contbin=contbin,yminor=0,xminor=0,axis_title_size=7,font_size=7, nshow=nshow, _extra=ex)


;saving of plots
  if not(keyword_set(nosave)) then begin
     outrootn=strmid(restorefile,0,strlen(restorefile)-4)+'.png'
     outfile='sdev_vs_par_'+outrootn
     psdev.save,outfile
     print,'saved sdev vs params plot to file ', outfile 


     outfile='parhist_'+outrootn
     parhi.save,outfile
     print,'saved histogram plot to file ', outfile 


     outfile='scatterplot_'+outrootn
     sm.save,outfile
     print,'saved parameter correlation plot to file ', outfile 
  endif                         ;nosave
end


