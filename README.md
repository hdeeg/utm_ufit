README for utm/ufit distribution    
Last update: 6Apr2021
(see further historial at end)

Short description of UTM to follow here





This directory and its subdirectories contain the full UTM/UFIT distribution. 
Its content is: 

/code   All required IDL code  

/documentation      Documentation of UTM/UFIT 

/examples       Files requires for a variety of example runs of UTM or UFIT

For more details, see _README files in the subdirectories

--
Citing UTM or UFIT: The preferred way is through its entry in the Astrophysics source code library 
at http://ascl.net/1412.003 
In ADS:  Deeg 2014ascl.soft12003D  


----
Versions of UTM/UFIT until 20210327 (2021 March 27) are available at
ftp://tep:fu9dufa5@ftp.iac.es/idl_hans_lib/utm/
The UTM/UFIT distro on that ftp server will not longer be updated and UTM versions posterior to 2021 March 27 will only be available at github.

For reference, below a list with the principal updates prior to 2021 March 27:
(For more detail, see _README files in the subdirectories and the headers of utm.pro, ufit.pro and other codes)

The 20210327 release contains several major updates:
- UTM uses now by default a hybrid calculation of the transit flux that is pixel-based or analytical, depending on the circumstances.
- UFIT permits the specification of uniform or gaussian priors to any parameter
- the directory-structure of the UTM/UFIT distributions has been re-organized. All code is now in the /code subdirectory

The 20151201 release increases the speed of UTM by several times and provides improved precision of the transit calculations.

The 20150917 release is a major update, for UTM.PRO, UFIT.PRO, documentation and examples. The principal improvements are:  
- All bodies can now be placed on arbitrarily orientated eccentric orbits; 
- All types of bodies can be orbiting other ones; 
- Interactive display of final model.
- UFIT uses AMOEBA fit as default and permits generation of MCMC chains. 
- The documention (utm_doc.txt and ufit_doc.txt) has been revised profoundly
- New example setup-files have been added, all old ones have been revised.





