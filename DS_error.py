import numpy as np
from astropy import units, constants
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import cosmology
plt.rcParams["font.family"] = "serif"

# H0=70  Omega_M = 0.3, physical units for R and DS
cosmo=cosmology.Cosmo(H0=70, omega_m = 0.3)

#DS profile
file = open("BolPlanck_r_DS_50.txt","r")
matrix=[[float(i) for i in line.split()] for line in file.readlines()]
lines=(np.array(matrix)).T
radius_mpc=lines[0]
DS=lines[1]
err_DS=[]
#Lens redshift
z_lens_min = 0.05
z_lens_max = 0.2
z_lens_mean = (z_lens_min+z_lens_max)/2
#Sources redshift
zs = 0.8
#cosmos = 1.64; HSC = 300
area_survey = 300
#Stellar mass range
Ms_min = 8
Ms_max = 9
#Source density in unit of galaxies/arcmin**2
source_d = 18
shape_noise = 0.28

#Stellar mass function from COSMOS2015
file = open("cosmos2015_dic2017_smf_z01-04_STY0.dat","r")
matrix=[[float(i) for i in line.split()] for line in file.readlines()]
lines=(np.array(matrix)).T
log_m=lines[0]
log_phi=lines[1]
log_phi_inf=lines[2]
log_phi_sup=lines[3]

SMF = interp1d(log_m, log_phi, kind=3)
x_SMF = np.linspace(0,13,10000)
Phi_interp1d = SMF(x_SMF)
SMF_new = interp1d(x_SMF, pow(10,Phi_interp1d), kind=3)

CSQUARE_OVER_4PIG = (constants.c ** 2/(4 * np.pi * constants.G)).to(units.Msun / units.pc).value

def sigma_crit(zl, zs, cosmo, comoving=False):
    """Calculate the distance term in Sigma_crit.
        
        Convert it into pc^1
        
        Parameter
        ---------
        zl : float or numpy array
        Redshift of lens.
        zs : float or numpy array
        Redshift of source.
        cosmos : cosmology.Cosmology object
        Cosmology object from `cosmology` package by Erin Sheldon.
        comoving : boolen, optional
        Flag for using comoving instead of physical unit. Default: False
        
        Return
        ------
        Critical surface density measurements
        
        """
    dist_term = ((1e-6 * cosmo.Da(0, zs) /
                  (cosmo.Da(zl, zs) * cosmo.Da(0, zl))))
        
    if comoving:
        return CSQUARE_OVER_4PIG * dist_term * (1.0 / (1. + zl)) ** 2

    return CSQUARE_OVER_4PIG * dist_term

sigma_crit_1 = 1/sigma_crit((z_lens_min+z_lens_max)/2, zs, cosmo, comoving=False)


def DS_error(z_lens_min, z_lens_max, Ms_min, Ms_max, source_d, shape_noise):
   i=0
   while i < len(DS)-1 :     
        # Bins in arcmin :
        am_1 = radius_mpc[i]/(cosmo.Da(0,z_lens_mean))*10800/np.pi
        am_2 = radius_mpc[i+1]/(cosmo.Da(0,z_lens_mean))*10800/np.pi
        # area
        area = np.pi*((am_2**2)-(am_1**2))
        # Total number of sources in this bin: 
        phiintegral= integrate.quad(lambda x: SMF_new(x), Ms_min, Ms_max)[0]
        #volume = integrate the voume givein zmin zmax and area
        volume = cosmo.V(z_lens_min,z_lens_max) * area_survey/(4*np.pi*(180/np.pi)**2)
        nlens = phiintegral * volume
        n_this_bin = area*source_d*nlens
        # Error on gamma
        err_gamma  = shape_noise/np.sqrt(n_this_bin)
        # Error on Delta Sigma
        err_ds = err_gamma/sigma_crit_1
        err_DS.append(err_ds)
        i=i+1
        #print(err_ds)

DS_error(z_lens_min, z_lens_max, Ms_min, Ms_max, source_d, shape_noise)
err_DS.append(0)
err_DS=np.array(err_DS)

plt.figure(figsize=(7.5,6))
x=np.array([-1000,1000])
y=np.array([0,0])
plt.plot(x,y,'k')
plt.errorbar(radius_mpc, DS,yerr=err_DS,fmt='o',markersize=5, capsize=5,ecolor='b',color='k')
plt.xscale('log')
plt.xlim(0.01,20)
plt.xlabel('R [Mpc]',fontsize=14)
plt.ylabel('$\Delta\Sigma$ [M$_\odot$ pc$^{-2}$]',fontsize=14)

