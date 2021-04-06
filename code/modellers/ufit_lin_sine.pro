pro ufit_lin_sine,setfname,dflag,xdat,ymod ;
;modeller for a superposition of a linear and a sinusoidal function, identical to the function used in Deeg et al. 2008 (CM Dra timing)
;with out-of bounds flag for negative values of kphi

  readsetupfile,setfname

; the coefficient namings from Deeg et al. 2008 (CM Dra timing) are used
  k0=double(parasign('k0'))
  k1=double(parasign('k1'))
  kd=double(parasign('kd'))
  kphi=double(parasign('kphi'))
  phi0=double(parasign('phi0'))
  ymod = k0 + xdat*k1 + kd * sin(xdat*kphi + phi0) ;function that is being evaluated
 ;raise the oobflag or reset to zero (if raised by previous invoke of model)
  if kphi lt 0 then paradd,"oobflag","1" else  paradd,"oobflag","0" 
end
