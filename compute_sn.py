
A python function that computes the expected SN of a detection. Core function to be called in the monte carlo routine.

Input :
z : redshift (central wavelength is computed from z)
lw : filter width
psf: which PSF to assume
texp : exposure time
galaxy model : M*, age, metallicity, etc ...
size : galaxy size
pizelsize : fix for now to 0.18"

option: which line to use (Halpha, OIII)
option: which profile to assume for dwarfs (e.g. exponential)

Assumes: a top hat filter of width lw.

STEPS OF THE CODE:

1) compute the background level for a given exposure time. To get the code up and running for now. Just fix the background level here for now.
I asked Kevin about computing the background in a given medium band filter. He did not think it would be too difficult. You compute the number of background photons given lw. Then assume poisson noise for the background. That gives the noise. There is one aspect to be careful here is that to make the postage stamp you want to know the number of background photons per pixel.
Kevin can help us more here.

2) compute the fluz of the galaxy in this mediaum band (given by z=> central wavelength and by lw)

3) generate a postage stamp for this galaxy,
Make sure that the size of the galaxy scales with the redshift assumed
Make a postage stamp that scales with the size of the galaxy (use 5 times the size of the galaxy?)

4) We assume that we know the centroid of the galaxy (forced photometry)
Compute aperture flux. Run SEP.
question: what aperture to use here? Use one that is matched to the size of the galaxy?

5) now we want to compute the error on the detected flux.

6) compute the predicted flux in this medium band filter WITHOUT the emission line.

7) Compute the SN (flux-broadbandflux)/sigma of the "detection" where "detection" is the fluz in 4) ABOVE the background level.

7) return the SN of the detection
(not of the detected flux, this is of the detected flux above he broadband level)

NOTE: keep note that in the future, we should include the error on the broadband flux


RETURNS : detected signal to noise of the emission line







