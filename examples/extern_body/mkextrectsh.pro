pro mkextrectsh
;test/demo to generate an external absorbing body of type 'extocc
;makes an upside rectangle


 ;2D array size; fixed for now for tests
  gsextx=100  ;gr-size extern obj in x
  gsexty=100   ;gr-size extern obj in y
                               
transarr=lonarr(100,100)
transarr[40:60,*]=150
;transarr[30:70,10:50]=250
save,transarr,filename='rectsharr.sav'
tv,transarr
;stop
end
