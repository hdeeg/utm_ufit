@ufit
pro ufitrandomrep,rootname,runname,randominit,nloop,gauss=gauss
;repeats ufit several times, with random offsets to initial parameters.
 
;rootname: name (without .utm) of the 'root' ufit setupfile
;runname: some additonal string to include in output files-name
;randominit: spread of random offset to initial paramters
;    For uniform distribution (default), the offseet ranges is - to + (randominit *scale)
;    for Gaussian distrib., the offset has a standard deviation (sigma) given by (randominit *scale )
;gauss:  if keyword set, then the offset has a gaussian distributons, else uniform.

;Example:
;ufitrandomrep,'trip_t18','u2',2,50
;repeats UFIT  50 times with setup trip_t18.utm. Applies each time to the initial values of fitted parameters a random offset of 2x their scale values, with uniform distribution. Output files will be named 'trip_t18u2.fout0.utm' to 'trip_t18u2.fout49.utm'

; HJD 25Sep 2015

  common parrs,pvarr,pnarr
  !except=0                     ;don't report math errors
 ; rootname='trip12.f_in1'
;randominit=1. 
 ; nloop=40
 versarr=fltarr(nloop)  ;.fout version-nr
  sdevarr=fltarr(nloop)  ;sdev
  for i=1,nloop do begin
     pnarr =" "
     pvarr =" "                  ;reset arrays
     readsetupfile,rootname+'.utm'  ;read root setup
     paradd,'randominit',randominit      ;amount of randomization
if keyword_set(gauss) then uniform_flag =0 else uniform_flag =1

 paradd,'random_uniform',uniform_flag      ;uniform or gauss  randomization dist.
version_offset=int(parasign('version_nr',0))  ;offset to version-nr

;     runname='.r'   ;name extension for that run (shouldn't be empty to avoid overwrite of original setup

     paradd,'version_nr',i+version_offset      ;version nr
     paradd,'plcflag','0'       ;no final plot
     paradd,' savenextf_in','0' ;don't save f_in file with next version Nr.
 
     setfname=rootname+runname+'.utm'
     savesetupfile,setfname  ;save intermediat setup
     sdev=0.
     ufit,setfname,1,sdev
     versarr[i-1]=i+version_offset
     sdevarr[i-1]=sdev
     writelc,rootname+runname+'sdev.txt',versarr,sdevarr ;write entire array in case it breaks
 ;    spawn,'\rm '+setfname
  endfor
end
