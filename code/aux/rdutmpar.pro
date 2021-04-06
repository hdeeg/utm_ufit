@setupfile   ;compile setupfile library
pro rdutmpar,filenametemplate,outfile,valarr
;reads parameters from multiple .utm setupfiles, given by filenametemplate, and saves to .csv file, whose root-name is given by outfile
;valarr is optional output string-array with the parameters' values
;Example: rdutmpar,'trip9c.fout*.utm','trip9c_pars'


;NOTE: The following line NEEDS MODIFICATION 
;rdpar=[ 'par1','par2','par3'  ]  ; list of parameters to extract

;HJD 29Sept 2015

  common parrs,pvarr,pnarr
  rdpar=['2omega','2omegadot','qmbin1','2period','2ecepoch','2ecc','2inclin','radAa2','radBaa1','2limbd','3pa','qlbin12','k1','2oblate','sdev','chisq','amoeba_ncalls_done'] ;parameters to extract


  npar=n_elements(rdpar)
  filist=file_search(filenametemplate) ;list of files found
  nfil=n_elements(filist)

  for i=0,nfil-1 do begin
     print,filist[i]
  endfor

  tmps=''
  read,'continue ?(y/n) ',tmps
  if strupcase(tmps) ne 'Y' then stop
  openw,UN1,outfile+'.csv',/get_lun
     printf,UN1,format='(a,",",$)','setupfile'
  for i=0,npar-1 do begin
     printf,UN1,format='(a,",",$)',rdpar[i]
  endfor
  printf,UN1,' '                ;force line return

  valarr = strarr(nfil,npar)    ; main array with extracted values
  for i=0,nfil-1 do begin
     readsetupfile,filist[i]
     printf,UN1,format='(a,",",$)',filist[i]

     for j=0,npar-1 do begin
        tmps=parasign(rdpar[j],'0')
        valarr[i,j]=tmps
        printf,UN1,format='(a,",",$)',tmps
     endfor
     printf,UN1,' '             ;force line return
  endfor
  free_lun,un1
end
