pro trip_prep3linld
common parrs,pvarr,pnarr
;preprocessing code for utm
;calculates individual component parameters for hierachical triple
;based on bin_prep.pro
;mass in Msol, period in day
;output rdist in Rsol

; inner binary Ba and Bb

snumA='2'  ;number of compA
snum1='3'  ;number of inner primary comp Ba
snum2='4'  ;number of inner secondary comp Bb

;luminosities
lumtot=double(parasign('lumtot'))  ;total luminos entire system
qlbin12=double(parasign('qlbin12'))  ;lum ratio  binary 1 / compA
qlbin1=double(parasign('qlbin1'))  ;lum ratio Bb/Ba in bin1
lumA=lumtot/(1+qlbin12) 
paradd,'2lum',lumA   ;lum A
lumbin1=lumtot*qlbin12/(1+qlbin12)  ;total lum of bin1
paradd,'3lum',lumbin1/(1+qlbin1)  ;lum Ba
paradd,'4lum',lumbin1*qlbin1/(1+qlbin1)  ;lum Ba


;masses
masstot=double(parasign('masstot'))  ;total mass entire system
qmbin12=double(parasign('qmbin12'))  ;mass ratio  binary 1 / compA
qmbin1=double(parasign('qmbin1'))  ;mas sratio Ba/Bb in bin1
massA=masstot/(1+qmbin12)
massbin1=masstot*qmbin12/(1+qmbin12)  ;total mass of bin1
massBa = massbin1/(1+qmbin1)  ;mass Ba
massBb = massbin1*qmbin1/(1+qmbin1)  ;mass Bb
paradd,'2mass',massA   ;mass A 
paradd,'3mass',massBa
paradd,'4mass',massBb

;semimajor axis of bin1
period1=double(parasign('3period'))
p1_y=period1/365.24
a1_AU=massbin1^(1D/3)*p1_y^(2D/3)  ;kepler eq. semimaj axes in AU
a1_bin=a1_AU * 214.95  ; relative semim-axis in Rsol
;print,'smaxis bin1:', a1_bin
aBa=a1_bin*qmbin1/(1+qmbin1) ;abs. axis Ba to Center of Mass
aBb=a1_bin*1/(1+qmbin1)
paradd,'3rdist',aBa
paradd,'4rdist',aBb

;semimajor axis of bin2
period2=double(parasign('2period'))
p2_y=period2/365.24
a2_AU=masstot^(1D/3)*p2_y^(2D/3)  ;kepler eq. semimaj axes in AU
a2_bin=a2_AU * 214.95  ; relative semim-axis bin2
;print,'smaxis bin2:', a2_bin
a2A=a2_bin*qmbin12/(1+qmbin12) ;abs. axis CM of all to CM of bin1
a21=a2_bin/(1+qmbin12); abs axis  CM of all of compA
paradd,'1rdist',a21 ;dist to reference point CM of bin1
paradd,'2rdist',a2A  ; to compA

;radii
radAa2 = double(parasign('radAa2')) ;radius A over smaxis bin2
radBaa1 = double(parasign('radBaa1')) ;radius A over smaxis bin2
;k2= double(parasign('k2'))  ;radius Ba / radius A
k1= double(parasign('k1'))  ;radius Bb / radius Ba
radA=radAa2*a2_bin
paradd,'2radi',radA
radBa=radBaa1*a1_bin
paradd,'3radi',radBa
radBb=k1*radBa
paradd,'4radi',radBb

;params taken from Ba and propagated to Bb
paradd,'4period',period1  ;
ecc1=double(parasign(snum1+'ecc'))
paradd,snum2+'ecc',ecc1  ;
inclin1=double(parasign(snum1+'inclin'))
paradd,snum2+'inclin',inclin1  ;
omega1=double(parasign(snum1+'omega'))
paradd,snum2+'omega',omega1+180.  ;
ecepoch1=double(parasign(snum1+'ecepoch'))
paradd,snum2+'trepoch',ecepoch1  ;transit of 2ndary=eclipse of prim
limbd1=double(parasign(snum1+'limbd'))  ;couple the limbd.
paradd,'4limbd',limbd1
quad=paratest('3limbd2') ;quad law used
if quad then begin
limbd21=double(parasign(snum1+'limbd2'))  ;couple the limbd.
paradd,'4limbd2',limbd21
endif
pa1=double(parasign('3pa'))
paradd,'4pa',pa1


; params from comp A and propagated to bin1 baryctr
paradd,'1period',period2  ;
ecc2=double(parasign('2ecc'))
paradd,'1ecc',ecc2  ;
inclin2=double(parasign('2inclin'))
paradd,'1inclin',inclin2  ;
omega2=double(parasign('2omega'))
paradd,'1omega',omega2+180.  ;
omegadot2=double(parasign('2omegadot'))
paradd,'1omegadot',omegadot2  ;

ecepoch2=double(parasign('2ecepoch'))
paradd,'1trepoch',ecepoch2  ;transit of 2ndary=eclipse of prim

end
