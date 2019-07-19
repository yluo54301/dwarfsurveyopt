This is a python code that optimized the lensing

FIX in a configuration file
- lensing survey parameters
- target emission line
- completeness C (how many galaxies at fixed M* you want to get). One number for the bin.
- M* range [mmin, mmax]. Default 10^8 and 10^9
- Nnights fixed number of nights (let's say 50 nights)
- assume size of imager (for HSC this is 1.7 deg^2)
- A telescope configuration (let's do HSC first, 8.2)                         
- hours per night: assume 8
- SN_detect: what we are willing to call a detection of an emission line, Let's call this SN_detect and assume SN_detect=4 as our minimum.

A for loop over three variables:                        
1. For a z array from 0.01 to 0.2 by deltaz of 0.01 (this correponds to the center of the filter)
2. Vary LW (between a typical narrrow band and a medium band)
3. Vary total (all dithers summed together) exposure time from low to high staritng with 5 min to 8h by steps of 10 min.
Select the minimum exposure time (just assume one long exposure) that satisfies completeness > C.
=> run_monte_carlo.py
(figure out how many times you need to draw in monte carlo to get a robust estimate of C)

If exposure time > number of nights then discard this lw value!

Compute area covered by survey in Nnights
    
For all values of LW that satisfy this criteria, compute the volume probed by the survey (volume is combination of area and depth of the survey as given by lw)

Compute the signal to noise of the lensing signal: SN_lens
For this we just need to fix a radial range, Let's use 0.1 to 5 Mpc
(the bining scheme doesn't matter here)                       
                             
(=> the code we wrote with Sukhdeep)
 
Make a 2d figure of Z (=central wavelength) versus Lw color coded by signal to noise
                         
Find the maximum value of SN_lens.

Write out:
- the maximum value of SN_lens
- lw
- targetted redshift range
- survey area








