;setupfile.pro is a collection of programs used to deal with setup-files
;
;version history
;3/5/1999 first collection of setupfile routines into setupfile.pro
;7/4/2000 added 'nowarn' keyword to parasign
;21/10/2002 added comment keyword to savesetupfile
;23/08/2005 added paradd routine
;20/10/2010 case insensitivity put as default to parasign; added casesens and warn kw
;11/09/2015  as comment, a pardelete that permits deletion (with wildcard) of parameters should be added. Useful for the saving of setup-files.
;29/11/2017 in readsetupfile, changed datstr and comstr from 200 to 400 elements
;09/03/2019  minor fixes only in comments lines, no code change
;26/03/2020  added parid kw to paratest
;10/4/2020  added foundflag to parasign

;the calling program needs to contain the following statements before
;accessing any of the programs in this file:
;   '@setupfile' 
;   common parrs,pvarr,pnarr


pro readsetupfile,setfname
;reads the setup-file 'setfname'
;and returns two arrays into the common block 'parrs':
;pvarr contains values of the parmeters
;pnarr contains names of the parameters
;HJD 1/1/1998
;HJD 29/11/2017  ;changed datstr and comstr from 200 to 400 elements
;we read the whole file into two arrays of strings
common parrs,pvarr,pnarr
datstr=strarr(400)               ;can deal with 400 parameters
comstr=strarr(400)               ;lines with comments
ndat=0l & ncom=0 & tstr= ' '
openr,UN1,/get_lun,setfname     ;do that here so file errors show up NOW
while(not eof(UN1)) do begin
    readf,UN1,tstr
                                ;in following, lines with a length of <= 2 characters or
                                ;lines beginning with '#' are considered comment lines.
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
pnarr=strarr(ndat)
pvarr=strarr(ndat)
;read in the main data block
for k=0,ndat-1 do begin
    tok=str_sep(strtrim(strcompress(datstr(k)),2)," ")
    if n_elements(tok) lt 2 then print,"error in setupfile with line: ",tok

    pnarr(k)=tok(0)
    pvarr(k)=tok(1)
;        print,pnarr(k)," : ",pvarr(k)
endfor
;end of reading stuff
end


function parasign,parname,default,nowarn=nowarn,casesens=casesens, warn=warn, foundflag=foundflag
;assigns a paramaeter 'parname' to its value. (extracts value of 'parname')
;used in combination with readsetupfile.pro
;The value that belongs to parname from 
;arrays pvarr and pnarr
;is returned in a variable, with the type of 
;pvarr, (which is normally of type string)
;
;the second parameter is optional and assigns a 
;default value if parname is absent. 
;Without a supplied default value,  
;parameters that aren't found will cause stopping.
;keywords
;   nowarn supresses printing of warning messages
;   warn allows printing of a non-generic warning message
;   casesens turns on up/low case-sensitivity to match parameter names; default is that case is ignored 
;   foundflag indicates if parameter was found (1) or if default was used (0)
;HJD 7/4/2000
;HJD 20/10/2010  Ignores case by default. With casesens kw set, will respect case. kw warn allows to give specific warning strings.
;HJD 13/5/2015  inserted warning if warn keyword is not a string
;HJD 10/4/2020  added foundflag
common parrs,pvarr,pnarr
foundflag = 0
if keyword_set(casesens) then begin ;case sensitive

for i=0,n_elements(pnarr)-1 do begin
    if pnarr(i) eq parname then begin
        parvalue = pvarr(i)
        foundflag = 1
    endif
endfor
endif else begin  ;ignore case (default)
for i=0,n_elements(pnarr)-1 do begin 
    if strlowcase(pnarr(i)) eq strlowcase(parname) then begin
        parvalue = pvarr(i)
        foundflag = 1
    endif
endfor

endelse

if keyword_set(nowarn) then wflag=0 else wflag=1
if keyword_set(warn) then if strlen(string(warn)) le 3 then print,"parasign.pro: If keyword warn is set, it should contain a warning message. If none wanted, use /nowarn"

if foundflag eq 0 then begin   ;no default specified. If parameter missing, stop

    if n_params() eq 1 then begin ; parameter not found, no default given
if keyword_set(warn) then print,warn else print,"STOPPED: parameter ",parname," not found in parameter file"
        stop                    
        endif
    if n_params() eq 2 then begin   ;parameter not found, use default
        if wflag then begin
if keyword_set(warn) then print,warn else $
        print,"WARNING: ",parname," not found in parameter file, using ", parname," = ",default
        endif
        parvalue=default
    endif
endif
return,parvalue
end


function paratest,parname,warn=warn,parid=parid
;tests, if parameter 'parname' exists in pnarr
;yes: returns 1, no: returns 0
;keywords:
;    warn: if set, returns warning if parname is absent
;    parid: if paramter is found, returns its index in pnarr. 
;used in combination with readsetupfile.pro
;HJD 27/4/1999
;    26/3/2020   added parid kw
common parrs,pvarr,pnarr
parid=[]
foundflag = 0
for i=0,n_elements(pnarr)-1 do begin
    if pnarr(i) eq parname then begin
    foundflag = 1 
    parid=i
    endif
endfor
if foundflag eq 0 and keyword_set(warn) then print,"WARNING: parameter ",parname," not found in setup-file"
return,foundflag
end


pro paraput,parname,value
;assigns a new value to existing parameter 'parname'
;
;Intended for use with program-modified setup files saved by 
;subsequent call to 'savesetupfile' 
;
;To add a new parmater, see paradd
;
;HJD 30/4/1999
common parrs,pvarr,pnarr
foundflag = 0
for i=0,n_elements(pnarr)-1 do begin
    if pnarr(i) eq parname then begin
        pvarr(i)=value
        foundflag = 1
    endif
endfor
if foundflag eq 0 then print,"Error: parameter: ",parname," not in setup file"
end



pro paradd,parname,value
;adds a new parameter-value pair to pvarr,pnarr
;if parameter already exits, modifies only the value and
; works identical to paraput
;
;value may be a variable of any kind, but take care to avoid
;conversion errors from floats to string in the builtin string procedure. If
;given number-format is needed, convert to string before calling paradd
;
;Intended for use with program-modified setup files saved by 
;subsequent call to 'savesetupfile' 
;
;HJD 23/8/2005
common parrs,pvarr,pnarr
foundflag = 0
for i=0,n_elements(pnarr)-1 do begin
    if pnarr(i) eq parname then begin
        pvarr(i)=value
        foundflag = 1
    endif
endfor
if foundflag eq 0 then begin
    pnarr=[pnarr,parname]
    pvarr=[pvarr,string(value)]
endif
end




pro savesetupfile,setfname,comment=comment
;saves setupfile from pvarr,pnarr
;used to save modified setupfiles
;comment lines are not saved in this version
;if comment keyword is set, setupfile is prepended by '#' comment symbols
;HJD 30/4/1999
;    21/10/2002 added comment keyword
common parrs,pvarr,pnarr
if keyword_set(comment) then comchar="# " else comchar=""
openw,unout,setfname,/get_lun
printf,unout,"# auto-saved setup file"
for i=0,n_elements(pvarr)-1 do begin
    printf,unout,comchar,pnarr(i),"  ",pvarr(i)
endfor
free_lun,unout
end
