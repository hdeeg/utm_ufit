pro kepler_prep
 common parrs,pvarr,pnarr
;example of preprocessing code for utm
;puts distance of object 1 (1rdist) into setupfile from 0mass and 1period
;simple version for m_star >> m_orbiter
;mass in Msol, period in day
;output rdist in Rsol
snumctr='0'  ;number of center-obj (star)
snumorb='1'  ;number of planet
mass=double(parasign(snumctr+"mass"))
period=double(parasign(snumorb+"period"))
p_y=period/365.24
a_AU=mass^(1D/3)*p_y^(2D/3)  ;kepler eq. semimaj axes in AU
a_Rsol=a_AU * 214.95 ;convert to Rsol
paradd,snumorb+"rdist",a_Rsol
end
