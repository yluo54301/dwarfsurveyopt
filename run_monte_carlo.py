Python code to run a monte carlo test

This monte carlo tests figures out the exposure time that satisfies completeness > C.

It monte carlos over:
- galaxy sizes (as drawn from COSMOS)
- galaxy model parameters (distribution of age, metallicity, etc ....)

This calls compute_sn.py

If compute_sn > SN_detect then we have detected the emission line.

Compute completeness as how many things we have detected over what went in.




