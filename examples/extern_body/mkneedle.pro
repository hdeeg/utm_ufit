pro mkneedle
;test/demo to generate an external absorbing body of type 'extocc'
;makes a needle shape


 ;2D array size; fixed for now for tests
  xs=100  ;gr-size extern obj in x
  ys=100   ;gr-size extern obj in y
                               
transarr=lonarr(xs,ys)
rad=49
for i=1,49 do begin
transarr[i,0+sqrt(rad^2-i^2):99-sqrt(rad^2-i^2)]=200
transarr[99-i,0+sqrt(rad^2-i^2):99-sqrt(rad^2-i^2)]=200

endfor
save,transarr,filename='needlearr.sav'
tv,transarr
end
