function scattermatrix,datarr,pname, bestpars=bestpars, nshow=nshow, contour=contour, contlevel=contlevel, contbin=contbin, text_orientation=text_orientation,axis_title_size=axis_title_size, noticklabel=noticklabel, nocorr=nocorr, _extra=e
;Create a matrix of scatter plots.
;based on correlplot vers  24oct2015
; input paramters:
;    datarr: array with [npar,ndat] elements, where npar is the number of parameters and ndat the number of data points
;   pname  : text-array with npar parameter-names
;output: plot-structure of scatterplot
; Keywords:
;    bestpar: vector with npar values that are inicated as red symbols. Intended to show the best-fit parameters but can be anything else
;    nshow: scatterplot will show only every n'th point so that not more than nshow are displayed. Default 
;            is 5000. nshow=-1 forces display of all points
;    contour: if set, contours are plotted instead of a scatterplot
;    contlevel: array of contour levels [0...1]. Default are 1 and 2 sigma probability contours, erf([1d,2d]/sqrt(2d0)) ~ [0.68,0.95]  
;    contbin: number of bins in X and Y direction that are used for 2D histograms in generation of conours; default = 100. 
;    text_orientation  : inclination of labels; default 45deg
;    axis_title_size  : font size of axis titles and of the correlation coefficients. Default is 12
;    noticklabel : if set, tick labels are supressed
;    nocorr : by default, correlation cefficients are printed in every panel. nocorr supresses this. 
;   _extra : other keywords accepted by plot function. Note that axis titles are drawn by text function, not by xtitle or ytitle keywords to plot
;
;  DEPENDENCIES:  (only for contour plot)
;  ufit_exofast_erell.pro
;
;
; Example:
;   Try the main-level example program at the end of this file::
;
;     IDL> .run scattermatrix
;written by Hans J Deeg
;
;version history
; 28feb2019 first dated version. 
; 8mar2019  further fixes; modified defaults for plots
; 11mar2019 added contbin keyword
; 13mar2019 added printing of correlation coefficients and nocorr keyword
; 26apr2019  fix in the logic of nshow
;  TODO: try contour plot from within errell, that might permit to color areas within contours

;
  siz=size(datarr)
  ndat=siz[2]

  datmin=min(datarr,dimens=2)   ; get min and max of each paramter
  datmax=max(datarr,dimens=2) 

  if not(keyword_set(contour)) then begin ; limit number of points in scatterplot
     if not(keyword_set(nshow)) then nshow=5000L ; max number of pts that are plotted
     if nshow eq -1 or ndat le nshow then stride =1 else begin    ;plot all pts
        stride = ((ndat-1)/nshow) +1 ;stride in plotting the datarr
        print,"Sample of ",strcompress(string(ndat))," points. Every ",strcompress(string(stride)),"'th plotted"
     endelse

  endif else begin              ;contour plot preps
     ;; 68% and 95% probability contours
     if n_elements(contlevel) eq 0 then contlevel = erf([1d,2d]/sqrt(2d0))
  endelse
  npar=siz[1]
  ngx=npar-1                    ;number of panels in x direction
  ngy=npar-1                    ;number of panels in y direction
                                ;    w=window(dimensions=[850,850])

  if keyword_set(noticklabel) then leftm=0.05 else leftm=0.1 ;margins, pending on presence of tick-labels
  if keyword_set(noticklabel) then botm=0.05 else botm=0.1
  rightm=0.03
  topm=0.03
  panelwi=(1.-leftm-rightm)/ngx
  panelhi=(1.-botm-topm)/ngy

;set non-standard defaults for plot (can be overriden through identical extra keywords in calling routine)
  font_size=6                   ;set default font size
  dimension=[600,600]           ;default window dimension 

  if not(keyword_set(text_orientation)) then text_orientation = 45 ;default text orientation
  if not(keyword_set(axis_title_size)) then   axis_title_size = 8 ;default axis title size


  corr = correlate(datarr)      ;correl. coefc

  current=0                     ;for first plot
  for j= npar-1, 1, -1  do begin ;j=4,3,2,1  for npar=5
     for i=0, j-1 do begin ;plot distributon of params against each other  ;i=0,1,2  for npar=4
        x0=leftm+i*panelwi
        x1=x0+panelwi
        y0=botm+(npar-j-1)*panelhi
        y1=y0+panelhi
        position=[x0,y0,x1,y1]
;       print, 'i=',i,', j=',j,', position=',position
        if keyword_set(contour) then begin ; contour plot
           ufit_exofast_errell, datarr[i,*],datarr[j,*],xpath=xpath,ypath=ypath,$
                                prob=contlevel,nxbin=contbin, nybin=contbin
           pd=plot(xpath,ypath,xrange=[datmin[i],datmax[i]],yrange=[datmin[j],datmax[j]],$
                   position=position,dimension=dimension,font_size=font_size,$
                   xshowtext=0,yshowtext=0,current=current,_extra=e) 
           if isa(bestpars) then pd=plot([bestpars[i]],[bestpars[j]],'+r', /over) ;overplot best paramter
        endif else begin        ;scatterplot
           pd=plot(datarr[i,0:*:stride],datarr[j,0:*:stride],'.',xrange=[datmin[i],datmax[i]],$
                   yrange=[datmin[j],datmax[j]],position=position,dimension=dimension,font_size=font_size,xshowtext=0,yshowtext=0,current=current,_extra=e)
           if isa(bestpars) then pd=plot([bestpars[i]],[bestpars[j]],'+r', /over) ;overplot best paramter
        endelse

; put axis-labels and titles depending on the panels' positions
; uses text function to write titles as it permits better positioning than in plot fct.
        if j eq npar -1 then begin ; lowest line, write labels
           if keyword_set(noticklabel)  then begin
              textypos=botm*0.5 ;adapt title position
           endif else begin      
              textypos=botm*0.2
              pd['axis0'].text_orientation = text_orientation ;hash-style referncing, see example to axis
              if i lt j-1 then begin ;supress last tick-label to avoid overlaps except for rightmost panel
                 tmp= pd['axis0'].tickvalues 
                 pd['axis0'].tickvalues=tmp[0:-2]
              endif
              pd['axis0'].showtext=1 ;show labels
           endelse
           t1=text(x0+panelwi*0.5,textypos,pname[i],align=0.5,vertical_align=0.5,font_size=axis_title_size) 
        endif
        if i eq 0 then begin    ; left column, write labels 
           if keyword_set(noticklabel) then begin
              textxpos=leftm*0.5 ;adapt title position
           endif else begin     
              textxpos=leftm*0.2 
              pd['axis1'].text_orientation = text_orientation
              if j gt 1  then begin ;supress top tick-label to avoid overlaps except for uppermost panel
                 tmp= pd['axis1'].tickvalues 
                 pd['axis1'].tickvalues=tmp[0:-2] 
              endif
              pd['axis1'].showtext =1 ;show labels
           endelse
           t1=text(textxpos,y0+panelhi*0.5,pname[j],align=0.5,vertical_align=0.5,orient=90,font_size=axis_title_size) 
        endif
        if not(keyword_set(nocorr)) then begin
           tc=text(x0+panelwi*0.5,y0+panelhi*0.9,string(format='(F0.3)',corr[i,j]),align=0.5,vertical_align=0.5,font_size=axis_title_size)
        endif
        current=1 
     endfor                     ;i
  endfor                        ; j
;  stop
  return,pd
;pd.save,'test.png'
end                             ;pro scattermatrix



; main-level example program
;generates npar distributions, each one consisting of 2 gaussion distribs with ndat points
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

sm=scattermatrix(datarr,pname, bestpars=bestpars, nshow=2000) ;scatterplot
sm2=scattermatrix(datarr,pname, xminor=0, yminor=0,  noticklabel=1, bestpars=bestpars, dimension=[700,700],axis_title_size=10, font_size=8,/contour ) ;with contours, showing some more options
savename='example_scatter_matrix.png'
stop
sm2.save,savename
print,'saved second plot to file ',savename
stop
end
