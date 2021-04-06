pro mkextsq
;test/demo to generate an external absorbing body of type 'extocc
;makes a square


 ;2D array size; fixed for now for tests
  gsextx=100  ;gr-size extern obj in x
  gsexty=60   ;gr-size extern obj in y
                               
transarr=lonarr(100,100)+150
;transarr[30:70,10:50]=250
save,transarr,filename='sqarr.sav'
;stop
end
