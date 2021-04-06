function parevol_plot, fitparname, MC_OUT_PARS , statvals=MC_OUT_CHISQ, nshow=nshow,chimax=chimax
;function parevol_plot, fitparname, j_pars , chainid, nchain,statvals=j_chisq
;  plots temporal evolution of parameters

;keywords:   
;    nshow: the plot will show only every n'th point so that not more than nshow are displayed. Default is 1000. nshow=-1 forces display of all points
;    chimax : the chisq plot will be scaled to max of chimax
;;
;   HJD 3may2019, still in prelim stage
;       11may2019   added chimax  
;       22mar2020  addapted to work with single chain
;        7aug2020  mod to show true point-number in X-labels 

  tmp= size(mc_out_pars,/dimensions)
  
  npars = tmp[0]
  ndat=long(tmp[1])
  if n_elements(tmp) eq 3 then nchain=tmp[2] else nchain=1

;  ptnr=lindgen(niter)                ;point nr
;  chainid= lindgen(niter) mod nchain ; vector with ids of chains

  
  if not(keyword_set(nshow)) then nshow=1000L ; max number of pts that are plotted
  if nshow eq -1 or ndat le nshow then stride =1 else begin ;plot all pts
     stride = ((ndat-1)/nshow) +1 ;stride in plotting the datarr
     print,"Sample of ",strcompress(string(ndat))," points. Every ",strcompress(string(stride)),"'th is plotted"
  endelse

                                ;number of panels to put side-by side
  if not(keyword_set(ncol)) then ncol=2 ;number of columns
  nrow=(npars+2)/ncol                   ;number of rows
  if isa(font_size) eq 0 then font_size=8 ;set default font smaller than the usual 12  

  loadct,39,rgb_table=rgbtab,ncolor=nchain+1 ;use  nchain +1 to avoid use of the last color index (is white in table 39) 


  current=0                     ;plot first panel in new window
  for i=0, npars-1 do begin     
     for j=0, nchain-1 do begin
        rgbcol=transpose(rgbtab[j,*])
        if j eq 0 then over=0 else over =1 
        pare=plot(MC_OUT_PARS[i,0:*:stride,j], layout=[ncol,nrow,i+1],current=current,over=over,ytitle=fitparname[i],font_size=font_size,symbol='dot',sym_color=rgbcol,color=rgbcol)
     endfor
;     pare=plotcolsym(ptnr,j_pars[i,*],chainid,layout=[ncol,nrow,i+1],current=current,xtitle=fitparname[i],symbol='dot',font_size=font_size)
     ax=pare.axes
     ax[0].COORD_TRANSFORM =[0,stride] ;X-labels as point-nr in  MC_OUT_PARS
     current=1
  endfor
  
  for j=0,nchain-1 do begin     ;plot chains versus chisq
     rgbcol=transpose(rgbtab[j,*])
     if j eq 0 then over = 0 else over=1
     pare=plot(MC_OUT_CHISQ[0:*:stride,j], layout=[ncol,nrow,i+1],current=current,over=over,ytitle='chisq',font_size=font_size,symbol='dot',sym_color=rgbcol,color=rgbcol)
  endfor
  ax=pare.axes
  ax[0].COORD_TRANSFORM =[0,stride] ;X-labels as point-nr in  MC_OUT_PARS
  if keyword_set(chimax) then pare.yrange=[min(MC_OUT_CHISQ),chimax]
end


