pro ufit_mccomb,mclist,outname_root,burn=burn
  mclist='lp141fitM_R5000savlist'
;combines mc.sav files from several runs of UFIT. By default, all links within the burn-in will be removed. Output file will be name root+'mcm.sav'
;the numbers of parralel chains and the number of fit-parameters in each mc.sav needs to be identical (There is no checking if the names of the fit-params agree)


;paramters: 
;mclist  : text-file containing file-names of mc.sav files, one per line.
;outname_root  : name of merged .mcm.sav file. Give name without .mcm.sav extension
;
;keywords:
;burn  : if set, the full chain including burn-in will be saved into merger. The burn-index (burnndx) of the first chain will then be transmitted to the merged chain.

  rdtab,mclist,filename,vtype=[7],ndat=ndat

  print,'Merging files:'
  MC_OUT_PARS_C=[]
  MC_OUT_CHISQ_C=[]
  for i=0,ndat-1 do begin
     print,filename[i]
     burnndx=[]                 ;reset in case it is missing in restored .sav
     restore,filename[i]
;restores these vars:
; BURNNDX         LONG      =          149
; FITPARNAME      STRING    = Array[4]
; MC_OUT_CHISQ    DOUBLE    = Array[20000, 8]
; MC_OUT_PARS     DOUBLE    = Array[4, 20000, 8]

     lastid=n_elements(where(MC_OUT_CHISQ[*,0] gt 0))-1-1 ;id of last filled element. Discard last step, since it might be incomplete (if generated from stop during demc.pro execution)
     if ~isa(burnndx) then begin
        print,"WARN: ",filename[i]," did not reach burn-in. All links are kept in merger"
        burnndx =0
     endif
     if keyword_set(burn) then begin          ;remove only unused last part of chain
        MC_OUT_PARS = MC_OUT_PARS[*,0:lastid,*] ;array of [npar,niter,nchain]
        MC_OUT_CHISQ = MC_OUT_CHISQ[0:lastid,*]
        if i eq 0 then burnndx_first=burnndx   
     endif else begin           ;remove burn-in and unused last part
        MC_OUT_PARS = MC_OUT_PARS[*,burnndx:lastid,*] 
        MC_OUT_CHISQ = MC_OUT_CHISQ[burnndx:lastid,*] 
     endelse
     tmp= size(MC_OUT_PARS,/dimensions)
     npar=tmp[0]
     nstep=long(tmp[1])
     nchain=tmp[2]
     
     MC_OUT_PARS_C=[[MC_OUT_PARS_C],[MC_OUT_PARS]] ;concatenate along 2nd dimension
     MC_OUT_CHISQ_C=[MC_OUT_CHISQ_C,MC_OUT_CHISQ]
  endfor
  MC_OUT_PARS=MC_OUT_PARS_C     ;rename back to original
  MC_OUT_CHISQ=MC_OUT_CHISQ_C  
  if keyword_set(burn) then burnndx=burnndx_first else burnndx =0  
  outmcname=outname_root+'.mcm.sav'
  save,mc_out_pars,mc_out_chisq,burnndx,fitparname,filename=outmcname 
  print,"MCMC chains merged: analyze with ufit_mcshow,'"+outmcname+"'"
end
