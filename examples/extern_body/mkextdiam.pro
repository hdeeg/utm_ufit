pro mkexdiam
;test/demo to generate an external absorbing body of type 'extocc'
;makes a diamond shape


 ;2D array size; fixed for now for tests
  xs=100  ;gr-size extern obj in x
  ys=100   ;gr-size extern obj in y
                               
transarr=lonarr(xs,ys)
for i=1,49 do begin
transarr[i,49-i:49+i]=200
transarr[99-i,49-i:49+i]=200
endfor
save,transarr,filename='diamarr.sav'
end
