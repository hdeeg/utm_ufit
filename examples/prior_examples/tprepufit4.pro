pro tprepufit4
;test preporc  without built-in priors
;pro wd_fitR_lomass
  common parrs,pvarr,pnarr
;preproc for WD scenario with free R_WD=r1 M_WD=f(R_WD), for 5000K WD 
;from tracing of Parsons+17, Fig 9left.


m5000=[0.26,0.4,0.53,0.68,1.]
r5000=[.02,0.0159,0.0136,0.0116,0.008]
;p=plot(m5000,r5000,'-')

;input params are m1(free),m2(fix),r1(free),krad(free),inclin(free)

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
  inclin=double(parasign(snum1+"inclin")) ;
  paradd,snum2+"inclin",inclin            ;



  r1= double(parasign(snum1+"radi")) ;in solar radii
  m1=interpol(m5000,r5000,r1) ; interpolate on 5000K R(M) relation
  m1=m1*0.9   ;lower the mass of the WD
  m2= double(parasign(snum2+"mass"))
  mtot=m1+m2  

  paradd,snum1+'mass',m1        ;mass of prim., for record-keeping
   

  krad = double(parasign("krad"))
  r2=r1*krad
  pyr=period/365.24
  a_bin_AU=pyr^0.66666666*mtot^0.333333333 ;semimaj axes in AU
  a_bin = a_bin_AU*214.94                   ;semimaj axes in Rsol
  rovera=   (r1+r2)/a_bin 

  paradd,snum2+'radi',r2        ;radius 2
  paradd,"(r1+r2)/a",rovera     ;output to setup use only for diagnostics



;luminosity
;lumtot=double(parasign('lumtot'))  ;total luminosity entire system
  lumtot = 1.                         ; normalized total luminosity
  contlum=double(parasign("contlum")) ;third light
  qlbin=double(parasign('qlbin'))     ;lum ratio comp2/comp1
  lumbin= lumtot - contlum            ;lumin of binary
  paradd,snum1+'lum',lumbin/(1+qlbin) ;lum 1
  paradd,snum2+'lum',lumbin*qlbin/(1+qlbin) ;lum 2

;outcomment below if a mass-ratio is known (will not change results of binary fit)
;qm=double(parasign("qmass"))  ;mass ratio q=m2/m1
;a_bin = 1  ;normalized relative semajor axes component 1 to comp 2

  qm=m2/m1
  a1=qm/(1+qm)*a_bin            ;semim halfaxis 1 to Center of Mass
  a2=1/(1+qm)*a_bin
;print,'a1,a2:',a1,a2
  paradd,snum1+"rdist",a1
  paradd,snum2+"rdist",a2

;the normalized relative semajor halfaxes of component 1 to comp 2 is 1. For equal masses , then:
;  paradd,snum1+'rdist',0.5      ;semi-maj. axis is 1 
;  paradd,snum2+'rdist',0.5


r1pri=0.011   ;gaussian prior for WD radius; mail Roi 16mar20
r1sig=0.0002    ;+-0.1Rearth   

ld1=double(parasign("0limbd"))
ld1pri=0.315
ld1sig=0.15

;chisqadd=((r1-r1pri)/r1sig)^2 + ((ld1-ld1pri)/ld1sig)^2
;paradd,'chisqadd',chisqadd

end
