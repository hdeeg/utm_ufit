;utm_ancil.pro 
;collection of anxiliary routines that required for execution of UTM and UFIT;(does not contain routines that are specific to either of these programs) 
;24jul2020 first creation, by moving routines below from UTM.pro and UFIT.pro


function frac, fnum
;returns the fractional part of a number
;example:
; print frac(-3.1)
;     0.9000000
;print frac(3.1)
;     0.1000000
;
;7 mar 1997
;hans joerg deeg
return, (fnum mod 1) +1.0 - ((fnum mod 1) ge 0)
end

function digirem, num, ndig
;reduces a number 'num' to one with ndig digits before or after decimal point (for negative ndig). 
;The leading sign is kept. Return value is a double precision number. 
;examples: digirem(2456789.23D,3)=789.23
;          (digirem(-456.34,1) = -6.34
;return value is a double precision number. 
;15Oct2015 HJD
fac=10D^ndig
numd=double(num)  ;subtraction in next line requires double;
if min(numd) lt 0D then print,"WARN: don't set inlcdigi for time-stamps with negative values"
return ,numd-fac*long(numd/fac)
;return,(numd ge 0D) ? numd-fac*long(numd/fac) : numd+fac*long(-numd/fac)
end


function str,num
;converts number to string and removes white space.
return,strcompress(string(num),/remove_all)
end


function e_omega,sqrtecsin,sqrteccos
; calcs eccentricity and omega from sqrt(e)*sin(omega) and sqrt(e)*cos(omega)   
;returns 2 element vector [eccentricity,omega] with omega from 0..360 degrees   
omega_rad=atan(sqrtecsin,sqrteccos)
sinomega=sin(omega_rad)
if abs(sinomega) ge 0.1 then $
ecc = (sqrtecsin / sinomega)^2 $
else ecc = (sqrteccos / cos(omega_rad))^2  ;this avoids use of small sinomega
; values. Also, ecc will always be positve       
omega=omega_rad*180/!dpi
if omega lt 0 then omega+=360  ;to get output from 0..360                       
return,[ecc,omega]  
end



function mean_anomal,nue,ecc, EccA=EccA, frac=frac, verbose=verbose ; calculates mean anomaly M (and eccentric anomaly E) from true anomaly nue and eccentricity ecc

;equations from wikipedia on eccentric anomaly
; input:  any value of nue (including multiples of 2 pi), eccentricity
;output: mean anomaly in radians; including epoch-count of nue (nepoch * 2pi)

; keywords:
; EccA   if  supplied as a defined variable, receives the eccentric anomaly.
; frac    if set, output of mean and eccentric anomaly is limited to 0 - 2pi. 
; verbose     verbose output

;H Deeg 16mar 2011  first version based on wikipedia 'eccentric anomaly'
;13apr11 inserted hook to set near-zero values of Ecca0 to 0
;30apr15 inserted frac keyword to limit output of M and Ecca to values 0 to 2pi

pi2 = 2*3.1415926535897932384626433832795028841971D 
nepoch=floor(nue/pi2)
sinE = sqrt(1-ecc^2)*sin(nue)/(1+ecc*cos(nue))
cosE = (ecc+cos(nue))/(1+ecc*cos(nue))
Ecca0=atan(sinE,cosE)
if abs(Ecca0) le 1e-11 then Ecca0=0D  ; avoid a num.error that leads to bad use of 'mod'  in the following line. Probably caused by limits in calculation of sin
Ecca=(Ecca0 +pi2) mod pi2 ;Ecc. Anom. -pi to pi is shifted to 0 to 2pi
M=Ecca  - ecc * sin(Ecca) ;M from 0 to 2pi
if not(keyword_set(frac)) then begin  ;add 2pi * epoch-count
M +=  nepoch*pi2
Ecca +=nepoch*pi2
endif
if keyword_set(verb) then print,'nue,nepoch,sinE,cosE,Ecca0,Ecca,M:',nue,nepoch,sinE,cosE,Ecca0,Ecca,M
return,M
end


function ecc_anomal,  M_in, ecc_in, prec=prec
; calculates the eccentric anomaly as function of the mean anomaly M_in 
;and eccentricity
;by inversion of M = E - ecc *sin (E) 
;prec is precision to achieve, in units of rad
;returns E in double-precision
;gets inefficient for ecc > 0.95
;HJD 30Apr2015
if not(keyword_set(prec)) then prec= 1e-7
M=double(M_in)
ecc=double(ecc_in)
En=M  
diff=1.
while abs(diff) ge prec do begin  ;Newtons Method,  see wikipedia Kepler's eq
E = En - (En - ecc*sin(En)-M)/(1-ecc*cos(En))  ;E is E_{n+1}
diff=  E - ecc *sin (E) - M
;print,M, En, E, diff
En = E  ;for next loop
endwhile
return, E
end

function true_anomal, E_in,ecc_in, mean_an=mean_an
; calculates the true anomaly as function of eccentric anomaly E and eccentricity. 
;keyword mean_an  : calculates from mean anomaly M. Input for E is then ignored
;works for any positive and negative E
;HJD 30Apr2015
E=double(E_in)
ecc=double(ecc_in)
pi2 = 2*3.1415926535897932384626433832795028841971D 
if keyword_set(mean_an) then begin
E=ecc_anomal(mean_an, ecc_in)
endif 
nepoch=floor(E/pi2)
E= E mod pi2  ;limit to 0 .. 2pi
if E lt 0 then E += pi2  ; to avoid negative E
;print,'nepoch, E, Efrac',nepoch, E, Efrac
polarx=sqrt(1-ecc)*cos(E/2.)
polary=sqrt(1+ecc)*sin(E/2.)
nue=2*atan(polary,polarx)  ;this eq is correct for -2pi < E < 2pi
nue+=nepoch*pi2
return,nue
end


function rdnumtab,filename,ncol
;reads numerical data into double-array, with as many colums as indicated by ncol
;with ignorance of headerlines (very short or starting with #)
datstr=strarr(5000000)           ;max number of data-pts   
comstr=strarr(2000)               ;max lines with comments
ndat=0l & ncom=0 & tstr= ' '
openr,UN1,/get_lun,filename     ;do that here so file errors show up NOW
while(not eof(UN1)) do begin
    readf,UN1,tstr
    if (strmid(tstr,0,1) eq '#') or strlen(tstr) le 2 then begin
        comstr(ncom)=tstr
        ncom=ncom+1
   endif else begin
        datstr(ndat)=tstr
        ndat=ndat+1
    endelse
endwhile
if ndat ge 1 then datstr=datstr(0:ndat-1)
if ncom ge 1 then comstr=comstr(0:ncom-1)
;print,ndat,' parameter values read from file ',setfname
;print,ncom,' comment lines read'
free_lun,UN1
;define data arrays
tmarr=dblarr(ndat,ncol)
;read in the main data block
for k=0,ndat-1 do begin
    tok=str_sep(strtrim(strcompress(datstr(k)),2)," ")
    if n_elements(tok) lt ncol then print,"error in file ",filename," with line: ",tok
    tmarr(k,*)=double(tok[0:ncol-1])
endfor
return,tmarr
end



function subshift, arr,dx,dy, padv=padv, padorg=padorg, padbil=padbil
;shifts the 2D-array 'arr' by subpixel offsets dx and dy, while maintaining flux. 
; Does not wrap, but by default pads with zero
; values. For conversation of flux or total count, the input array should have an outer 1-pix wide frame
; of zero values. Else, the keywords provide some other options for padding

;input:
; arr       a 2D array
; dx,dy   amount of offset in x and y direction. The absolute value of dx and dy must be less or equal to 1 
;
;output
;   the shifted array. it will be of float or double type. 
;
;keywords 
;     padv     value that is used for padding, instead of the default of 0
;     padorg  if set, padding is done with the orginal array values 
;     padbil   if set, padding is done using a bilinear interpolation within the row or colum to be padded.
;               The one corner-point from which the array is shifted away from is being filled with the 
;               original value
;
;  HJD 5 jun2019  first version
;      12jun 2019 added padorg keyword
;      25 jun2019 added padbil keyword

  if ~isa(padv) then padv = 0

  sx=signum(dx)                 ;sign of x, y -shift
  sy=signum(dy)
  dx=abs(float(dx))
  dy=abs(float(dy))

;  arrtyp=size(arr,/type)

  sharr_x=shift(arr,sx,0)          ;array shifted a pix in x
  if sx eq 1 then sharr_x[0,*]=padv ;fill wrapped rows/cols
  if sx eq -1 then  sharr_x[-1,*]=padv  

  sharr_y=shift(arr,0,sy)       ;array shifted a pix in y
  if sy eq 1 then sharr_y[*,0]=padv 
  if sy eq -1 then sharr_y[*,-1]=padv  

  sharr_xy=shift(arr,sx,sy)     ;array shifted a pix in x and in y
  if sx eq 1 then sharr_xy[0,*]=padv 
  if sx eq -1 then sharr_xy[-1,*]=padv  
  if sy eq 1 then sharr_xy[*,0]=padv 
  if sy eq -1 then sharr_xy[*,-1]=padv  
  
;add all shifted arrays for final convoluton
  sharr=(1.-dx)*(1.-dy)*arr+dx*(1.-dy)*sharr_x+(1.-dx)*dy*sharr_y+dx*dy*sharr_xy
;print,dx,dy,sx,sy,total(sharr)

if keyword_set(padorg) then begin  ;pad with original values
if sx eq 1 then sharr[0,*]=arr[0,*]
if sx eq -1 then sharr[-1,*]=arr[-1,*]
if sy eq 1 then sharr[*,0]=arr[*,0]
if sy eq -1 then sharr[*,-1]=arr[*,-1]
endif

  if keyword_set(padbil) then begin ;pad with original values using a bilinear interpol in the 
;corresponding row or colum
     if sx eq 1 then sharr[0,*]=dy*shift(arr[0,*],sy)+(1.-dy)*(arr[0,*]) ;shift left-most column up or down
     if sx eq -1 then sharr[-1,*]=dy*shift(arr[-1,*],sy)+(1.-dy)*(arr[-1,*]) ;right-most col
     if sy eq 1 then sharr[*,0]=dx*shift(arr[*,0],sx)+(1.-dx)*arr[*,0]       ;top row
     if sy eq -1 then sharr[*,-1]=dx*shift(arr[*,-1],sx)+(1.-dx)*arr[*,-1]    ; bottom row
     if sx ne 0 and sy ne 0 then begin ;we need to pad a corner point
        sx = sx < 0                    ;to set it to -1 or 0
        sy = sy < 0  
;the one corner-point from which the array is shifted away from is being filled with the original value
        sharr[sx,sy]=arr[sx,sy]
     endif
  endif                         ;if keyword_set(padbi)

  return, sharr
end
 
