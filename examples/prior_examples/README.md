Examples of UFIT setups with priors are given.\
All of them fit a model to a lightcurve of an eclipse of a white dwarf 

A gaussion prior of 0radi=0.0108±0.002 and an uniformed prior that limits the linear limbdarkening 0limbd between 0.165 and 0.465  have been set:

IDL> ufit,'prior_in_setup.utm',2,/nofit\
Has the priors defined at the end of the setupfile

IDL> ufit,'prior_in_preproc.utm',2,/nofit\
defines priors and calculates the corresponding chisq in the preprocessor tprepufit1.pro  (in last few code lines; else it is identical to the preprocessor tprepufit4.pro)

In both cases, the output should be:
```
------------ufit: initial parameters--------------
krad :       12.843492
0radi :     0.010603270
0inclin :       88.345625
0limbd :      0.29887926
0ecepoch :      0.72061889
chisq (all pts)    :       27.466706  meas_error:    0.00500000
chisq_red (all pts):      0.45027386
AICc               :       37.466706  BIC       :        48.490170
stdev (all pts)    :    0.0032255317  npt_all   :           67
stdev (off-eclipse):    0.0033244054  npt_off   :           35
stdev (on-eclipse) :    0.0031669198  npt_on    :           32
chisq (on-eclipse) :       12.436433
chisq_red (on-ecl) :      0.47832433
chisq (priors)     :      0.96756732
chisq (total)      :       28.434273
------------------------ufit finished-------------------------
````
IDL> ufit,'noprior.utm',2,/nofit  
Is like above, but without priors. The last two lines in above output will be missing (the total chisq is now the value indicated in 'chisq (all pts)'

Without the /nofit keyword, a very short MCMC run will be initiated.
For a useful MCMC run, give larger values to these setup parameters:  
mc_nchains  (e.g 10 resp. twice the number of fit-parameters)  
and  
mc_maxsteps  (e.g 1000, better 5000)






