pro ufit_polynome,setfname,dflag,xdat,ymod ;
;modeller for polynomes
  readsetupfile,setfname
  order=int(parasign("order",1,/warn)) ;order of ploynome
  cof=dblarr(order+1)           ;coefficients of polynome
  for i=0, order do begin       ;read coefficents
     cofname='c'+strcompress(string(i),/rem)
     cof[i]=double(parasign(cofname))
  endfor
  ymod=poly(xdat,cof)           ;function that is being evaluated
end
