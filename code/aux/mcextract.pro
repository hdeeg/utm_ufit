;template code to 
;extract some specific values from MC chain
;(here, extracts some values which are high in 0mass and have good chisq)
restore,'M_R5000_no4.mcm.sav'
;fitparname:  krad 0mass 0inclin 0ecepoch
   tmp= size(MC_OUT_PARS,/dimensions)
  npar=tmp[0]
  nstep=long(tmp[1])
  nchain=tmp[2]
  niter=nchain*nstep   

 j_pars = reform(transpose(MC_OUT_PARS,[0,2,1]), npar, niter) ;join all chains and put dims in order par, step, chain
  j_chisq = reform(transpose(MC_OUT_CHISQ) , niter)            ; order step, chain
  chainid= lindgen(niter) mod nchain  

himid = where(j_pars[1,*] ge 0.787 and j_chisq le 15)  ;1 param-set remains
p=plot(j_pars[1,himid],j_chisq[himid],'+')

print,fitparname
print,j_pars[*,himid[0]]
print,'chisq=',j_chisq[himid[0]]
stop
end
