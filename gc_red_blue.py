from astropy.table import Table
from math import *
import matplotlib.pyplot as plt 
import scipy.stats
import numpy as np
from gc_stat import *

cat = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat = cat[cat['p_gc']>0.95] # sort out the GCs
#cat = cat[cat['gmag']<24.5] 
cat_gc_b = cat[cat['gi0']<0.8] 
cat_gc_r = cat[cat['gi0']>0.8] 
M87 = (187.70583, 12.39111)
DIS = 17.21*1000 # in kpc (luminosity distance from NED)
radius_max = 500. #kpc   max radius for bin statistics (gc and ucd)

#Binned Statistics
gc_density_b, gc_density_err_b, bin_mean_gc_b = gc_density_stat(cat_gc_b)
gc_density_r, gc_density_err_r, bin_mean_gc_r = gc_density_stat(cat_gc_r)
'''
r_maj = np.exp(np.arange(0.955,7.,0.05))
r_maj = r_maj[r_maj<radius_max] 
e = 0.5
gc_density_b, gc_density_err_b = elliptical_stat(cat_gc_b, r_maj,e=e)
gc_density_r, gc_density_err_r = elliptical_stat(cat_gc_r, r_maj,e=e)
'''

'''==================fitting=========================='''
#fit range
fit_min = 25   # kpc
fit_max = 210  # kpc

slope_b, intercept_b, r_value_b = gc_fitting(gc_density_b, bin_mean_gc_b, fit_min, fit_max)
abline_values_b = [radius**slope_b*exp(intercept_b) for radius in np.arange(1,1000)] # parameter for drawing the fitted line

slope_r, intercept_r, r_value_r = gc_fitting(gc_density_r, bin_mean_gc_r, fit_min, fit_max)
abline_values_r = [radius**slope_r*exp(intercept_r) for radius in np.arange(1,1000)] # parameter for drawing the fitted line

print 'blue slope:',slope_b,'blue intercept:',intercept_b,'r-squared:', r_value_b**2
print 'red slope:',slope_r,'red intercept:',intercept_r,'r-squared:', r_value_r**2

'''===================plot============================'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(bin_mean_gc_b,gc_density_b,yerr=gc_density_err_b,fmt='.b')
ax.errorbar(bin_mean_gc_r,gc_density_r,yerr=gc_density_err_r,fmt='.r')
ax.plot(np.arange(1,1000), abline_values_b, 'b',label = r'$n_{\rm blue}=$'+str(round(slope_b,2)))
ax.plot(np.arange(1,1000), abline_values_r, 'r',label = r'$n_{\rm red}=$'+str(round(slope_r,2)))
ax.axvline(x=fit_min, ymin=0, ymax = 10, linewidth=1, linestyle='dashed', color='k')
ax.axvline(x=fit_max, ymin=0, ymax = 10, linewidth=1, linestyle='dashed',color='k')
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlabel(r'$R_{\rm maj} {\rm [kpc]}$',fontsize=15)
ax.set_ylabel(r'Surface density [N/kpc$^2$]')
ax.set_xlim([1,500])
ax.set_ylim([0.0001,10])
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(numpoints=1,frameon=False,loc='lower left')
plt.title('fit range:'+str(fit_min)+'~'+str(fit_max))
plt.savefig('pics/gc_red_blue.png')
plt.show()