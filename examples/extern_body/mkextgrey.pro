pro mkextgrey
;test/demo to generate an external absorbing body of type 'extocc'


 ;2D array size; fixed for now for tests
  gsextx=100  ;gr-size extern obj in x
  gsexty=100   ;gr-size extern obj in y
                               
transarr=lonarr(100,100)
for i=0,49 do transarr[i,*]=i*5
for i=0,49 do transarr[99-i,*]=i*5
save,transarr,filename='greystripe.sav'
tv,transarr
end
