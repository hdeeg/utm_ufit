@setupfile   ;compile setupfile library
pro utmparplot,filenametemplate ,nosave=nosave ;
;reads parameters from multiple .utm setupfiles, given by filenametemplate,
; and makes i) plots of params versus stat. value and ii) param-param correlation plots
;keywords
;  nosave  :if set, plots are not saved to png files


;Example: utmparplot,'trip9c.fout*.utm'
  
; NOTE:  These lines NEED MODIFICATION 
; pname=['sdev','2omega','2omegadot ...]:  for desired params. Keep sdev (or chisq ) first!
;  statlimit =0.003  ;limit plotted values to sdev (or chisq) below this value. ;  statlimit = 0   ; include all points


;HJD 24 Oct 2015
;26Feb 2019: changed plotting program to scattermatrix.pro
;8mar2019: externalized sdev versus parameter plots by plot_par1_vs_params.pro.
;         added histogram plot and nosave keyword


  common parrs,pvarr,pnarr
; pname=['sdev','thirdlight','radsuma','0inclin','0ecepoch','0ecc','srat']
  
  pname=['sdev','2omega','2omegadot','qmbin1','2period','2ecepoch','2ecc','2inclin','radAa2','radBaa1','2limbd','3pa','qlbin12','k1'] ;parameters to extract.  Keep sdev (or chisq ) first!


; statlimit =0.00011              ;limit plotted values to sdev (or chisq) below this value

  statlimit =0.00

  npar=n_elements(pname)
  filist=file_search(filenametemplate) ;list of files found
  nfil=n_elements(filist)

  print,'Selected files: '
  for i=0,nfil-1 do begin
     print,filist[i]
  endfor

  parr = dblarr(npar,nfil)      ; main array with extracted values
  for i=0,nfil-1 do begin
     readsetupfile,filist[i]

     for j=0,npar-1 do begin
        tmps=parasign(pname[j],'0')
        parr[j,i]=double(tmps)
     endfor
  endfor
  sdev=parr[0,*]                ;stat value (sdev or chisq)
;size=parr[1,*]

  if statlimit eq 0.0 then goodid = lindgen(nfil) else $
;     goodid=where(sdev le statlimit and size le 0.2) ;ids of good fits
             goodid=where(sdev le statlimit) ;ids of good fits

  parr=parr[*,goodid]           ;array with good values
;  sdevg=sdev[goodid]            ;  'good' stat values

  filistg = filist[goodid]      ;'good' param files

  print,'i: plots of params against stat. value'
  psdev=plot_par1_vs_params(pname,parr)


  print,'ii: parameter histograms'
  garr=parr[1:*,*]              ;parameter values without the sdev
  pname=pname[1:*]              ;same for parnames

  parhi=ufit_exofast_plotdist(garr, parnames=pname,nbins=10,dimension=[800,700])

  print,'iii: parameter correlation plots'
  
  sm=scattermatrix (garr,pname, /noticklabel,axis_title_size=6,xminor=0, yminor=0, dimension=[700,700])


;saving of plots
  if not(keyword_set(nosave)) then begin
     outrootn=strmid(filenametemplate,0,strlen(filenametemplate)-1)+'.png'
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
