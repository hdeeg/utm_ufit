function plot_par1_vs_params, pname_in,parr_in,statvals=statvals,statname=statname,ncol=ncol, nshow=nshow,_extra=ex
;generates scatterplots of  the first parameter in parr[] array (intended to be chisq or stddev) versus all other parameters' values. 
;
; input:
;pname_in: a string vector with  N parameter names. pname[0] gives the type of statistical parameter ( 'chisq' or 'sdev', but can be anything) 
;parr_in:  an N x M array, of M values for each of the N paramenters
;
;keywords:
; statvals : input vector with M statisticval values; all values of 
;    parr are plotted against it  
; statname : input string with name of the statistical values ( 'chisq' or 'sdev', but can be anything). statvals and statname must be set together 
; 
;  ncol : the number of columns in which panels are plotted. Default is 4
;  nshow: scatterplot will show only every n'th point so that not more than nshow are displayed. Default is 10000. nshow=-1 forces display of all points
;   _extra : other keywords accepted by plot function.
;

; Example:
;   Try the main-level example program at the end of this file::
;
;     IDL> .run plot_par1_vs_params
;written by Hans J Deeg
; 
;version history
;8mar2019  first version, derived from utmparplot7mar19.pro
;10mar2019  added statvals and statname keywords
;13mar2019 added nshow keyword, similar to scattermatrix
;
; TODO: use more sophisticated panel positonining, like in scattermatrix.pro

  parr=parr_in            ;do it this way to avoid modification of input arrays
  pname=pname_in

  if keyword_set(statvals) then begin
     if isa(statname) eq 0 then begin
        print,'statname keyword must be given, since  statvals keyword is used'
        stop
     endif
     parr=[transpose(statvals),parr] ;'glue' the statvals into the first column
     pname=[statname,pname]          ;same with statname
  endif                              ;keyword_set(statvals)


  npar=n_elements(pname)-1 ;number of independent paramters (=number of panels)
                                ;number of panels to put side-by side
  if not(keyword_set(ncol)) then ncol=4 ;number of columns
  nrow=(npar-1)/ncol +1                 ;number of rows

;set non-standard defaults for plot (can be overriden through identical extra keywords in calling routine)
  font_size=7                   ;set default font size
  dimension=[700,600]           ;default window dimension 

  siz=size(parr)
  ndat=siz[2]                                      ;number of data points
  if not(keyword_set(nshow)) then nshow=10000L ; max number of pts that are plotted
  if nshow eq -1 or ndat le nshow then stride =1 else begin  ;plot all pts
     stride = ((ndat-1)/nshow) +1 ;stride in plotting the datarr
     print,"Sample of ",strcompress(string(ndat))," points. Every ",strcompress(string(stride)),"'th plotted"
  endelse

  for i=1,npar do begin                 
     if i eq 1 then current =0  else current=1 ;to start new win

     if i mod ncol eq 1 then begin ;identify leftmost plot
        ytitle=pname[0]   
        margin=[0.21,0.15,0.02,0.03] ; leftmost plot
        yshowtext=1
     endif  else begin 
        ytitle=""         
        margin=[0.15,0.15,0.02,0.03] ;other plots
        yshowtext=0
     endelse
     pl=plot(parr[i,0:*:stride],parr[0,0:*:stride],symbol='+',linestyle='none',layout=[ncol,nrow,i],xtitle=pname[i],ytitle=ytitle,margin=margin,yshowtext=yshowtext, font_size=font_size,dimension=dimension,current=current,_extra=ex) ;do explicit spec of symbol and linestyle to permit overide by _extra.
     pl['axis3'].showtext= 0 ;fudge to remove labeling from right y-axes if yshowtext =1 
                                ;  ax=pl.axes & ax[3].showtext=0 ;fudge to remove labeling from right y-axes if yshowtext =1 
  endfor

  return,pl
end


; main-level example program
;generates npar distributions, each one consisting of 2 gaussion distribs with ndat points and calls the plotting routine
npar=15  ;will generate npar-1 panels, since par0 will be y-axis on the others
npts=1000  ;number of points 

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

medianpar=0.

d1plot=plot_par1_vs_params(pname,datarr)  ;plot with default settings
;call below shows also that any further plot parameters might be passed as keywords
d2plot=plot_par1_vs_params(pname,datarr,dimension=[600,800],font_size=6,yminor=0,xminor=0,ncol=3,symbol='period')
savename='example_distributions.png'
d2plot.save,savename
print,'saved plot to file ',savename
;  for i=0,npar-1 do begin
;     print,pname[i],': best:',bestpars[i]
;  endfor
end


