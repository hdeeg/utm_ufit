pro rpl_apl_prep
common parrs,pvarr,pnarr
;preprocessing code for utm

;star is object 0
; permits use of input parameters 'rpl/rst' and 'a/rst' (planet/star radius-ratio and distance/star-radius) for a given star-radius 
;below example for two planets
;print,'preprocessor'
paradd,'1radi',double(parasign('1rpl/rst')) * double(parasign('0radi'))
paradd,'1rdist',double(parasign('1a/rst')) * double(parasign('0radi'))

;same for planet2
paradd,'2radi',double(parasign('2rpl/rst')) * double(parasign('0radi'))
paradd,'2rdist',double(parasign('2a/rst')) * double(parasign('0radi'))

end
