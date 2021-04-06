pro bin_prep_k2
  common parrs,pvarr,pnarr
;standard preprocessor for lightcurves using normalized fit-parameters
;the period might be 1 for phase-curves or the period in real units.
; uses radius ratio r1/r2 and (r1+r2)/a

;propagates some params from primary to secondary
;converts relative luminosity into  normalized ones for each component
;may also take a mass ratio if known(outcommented code)

  snum1='0'                     ;number of primary comp 1
  snum2='1'                     ;number of secondary comp 2


;params that are propagated from primary to secondary comp:
  period=double(parasign(snum1+"period"))
  paradd,snum2+"period",period  ;
  ecc=double(parasign(snum1+"ecc"))
  paradd,snum2+"ecc",ecc        ;
  omega=double(parasign(snum1+"omega"))
  paradd,snum2+"omega",omega+180. ;omega 2ndary = 180 + omega primary
  ecepoch1=double(parasign(snum1+"ecepoch"))
  paradd,snum2+"trepoch",ecepoch1 ;transit of 2ndary=eclipse of prim
;  limb1=double(parasign(snum1+"limbd"))   
;  paradd,snum2+"limbd",limb1    ;couple the limbd.


;convert k=r2/r1 and  (r1+r2)/a to radii
krad = double(parasign("krad"))
rovera = double(parasign("(r1+r2)/a")) ; (R1+R2)/a
r1=rovera / (1+krad)
r2=rovera / (1 + 1/krad)
 paradd,snum1+'radi',rovera / (1+krad)   ;radius 1
 paradd,snum2+'radi',rovera / (1+1/krad)   ;radius 2

;propagate inclination from primary to secondary comp, and mirror
; at 90deg if needed
  inclin=double(parasign(snum1+"inclin")) ;
  paradd,snum2+"inclin",inclin  ;

;luminosity
  lumtot = 1.                        ; normalized total luminosity
  contlum=double(parasign("contlum")) ;third light
  qlbin=double(parasign('qlbin'))     ;lum ratio comp2/comp1
  lumbin= lumtot - contlum            ;lumin of binary
  paradd,snum1+'lum',lumbin/(1+qlbin) ;lum 1
  paradd,snum2+'lum',lumbin*qlbin/(1+qlbin) ;lum 2

;outcomment below if a mass-ratio is known (will not change results of binary fit)
;qm=double(parasign("qmass"))  ;mass ratio q=m2/m1
;a_bin = 1  ;normalized relative semajor axes component 1 to comp 2
;a1=qm/(1+qm)*a_bin ;semim halfaxis 1 to Center of Mass
;a2=1/(1+qm)*a_bin
;print,'a1,a2:',a1,a2
;paradd,snum1+"rdist",a1
;paradd,snum2+"rdist",a2

;the normalized relative semajor halfaxes of component 1 to comp 2 is 1. For equal masses , then:
  paradd,snum1+'rdist',0.5      ;semi-maj. axis is 1 
  paradd,snum2+'rdist',0.5

end
