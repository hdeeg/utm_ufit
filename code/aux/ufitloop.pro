@ufit
pro ufitloop
;loops UFIT with input setups varied by some initial value.
;in this case, 1radi from 0.05 to 0.4 
  common parrs,pvarr,pnarr
  !except=0                     ;don't report math errors
  rootname='2041radivar'
  rad0=0.05
  radmax=0.4

 for rad=rad0,radmax,0.005 do begin
  readsetupfile,rootname+'.utm'  ;note: if readsetupfile is done before loop, the next loop will use the final parameters of the current one.
  paradd,'savenextf_in','0'     ;no second output setup
  paradd,'showfinalfitflag','1'          ;no final plot of fit
     print,'LOOP: fitting for 1radi: ',rad
     paraput,'1radi',rad
     setfname=rootname+'.rad'+strcompress(string(fix(rad*1000)),/rem)+'.utm'
     savesetupfile,setfname
     ufit,setfname,1
  endfor
end
