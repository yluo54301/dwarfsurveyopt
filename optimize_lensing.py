This is a python code that optimized the lensing

FIX in a configuration file:
- lensing survey parameters
- target emission line
- completeness C (how many galaxies at fixed M* you want to get). One number for the bin.
- M* range [mmin, mmax]. Default 10^8 and 10^9
- Nnights fixed number of nights (let's say 50 nights)
- assume size of imager (let's say 1 deg^2)
- hours per night: assume 8
- SN_detect: what we are willing to call a detection of an emission line, Let's call this SN_detect and assume SN_detect=4 as our minimum.

For a z array from 0 to 0.2.
Vary LW.
Select exposure time (just assume one long exposure) that satisfies completeness > C.
=> run_monte_carlo.py

If exposure time > number of nights then discard this lw value!

Compute area covered by survey in Nnights
    
For all values of LW that satisfy this criteria, compute the volume probed by the survey (volume is combination of area and depth of the survey as given by lw)

Compute the signal to noise of the lensing signal: SN_lens
(=> the code we wrote with Sukhdeep)
 
Find the maximum value of SN_lens.

Write out:
- the maximum value of SN_lens
- lw
- targetted redshift range
- survey area








