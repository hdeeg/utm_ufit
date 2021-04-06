pro mkextest
;test/demo to generate an external absorbing body of type 'extocc'


 ;2D array size; fixed for now for tests
  gsextx=100  ;gr-size extern obj in x
  gsexty=100   ;gr-size extern obj in y
                               
transarr=(lindgen(gsextx,gsexty) mod (gsextx))*2
save,transarr,filename='transarr.sav'
tv,transarr
end
