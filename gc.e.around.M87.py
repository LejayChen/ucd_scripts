from astropy.table import Table
from math import *
import matplotlib.pyplot as plt 
import scipy.stats
import numpy as np
from gc_stat import elliptical_stat
np.set_printoptions(precision=2)

cat = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat = cat[cat['p_gc']>0.95] # sort out the GCs
#cat = cat[cat['gmag']<24.5] # sort out the GCs
cat_ucd = Table.read('NGVS.pilot.92ucds.fits')
M87 = (187.70583, 12.39111)
DIS = 17.21*1000 # in Kpc (luminosity distance from NED)
radius_max = 500. #kpc   max radius for bin statistics (gc and ucd)

def gc_fitting(gc_density, bin_mean_gc, fit_min, fit_max):

	gc_density = gc_density[~np.isnan(bin_mean_gc)]
	bin_mean_gc = bin_mean_gc[~np.isnan(bin_mean_gc)]

	gc_density_fit = gc_density[bin_mean_gc>fit_min]
	bin_mean_gc_fit = bin_mean_gc[bin_mean_gc>fit_min]

	gc_density_fit = gc_density_fit[bin_mean_gc_fit<fit_max]
	bin_mean_gc_fit = bin_mean_gc_fit[bin_mean_gc_fit<fit_max]

	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))

	return slope, intercept, r_value

'''================Bin Statistics for gc====================='''
r_maj = np.exp(np.arange(0.955,7.,0.05))
r_maj = r_maj[r_maj<radius_max] 
e = 0.5
gc_density, gc_density_err = elliptical_stat(cat, r_maj,e=e)

'''================== Bin Statistics for ucd ======================'''
dis_list_ucd = []
for i in range(len(cat_ucd)):
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	distance_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)  # in degree
	distance_ucd = distance_ucd/180.*pi*DIS #in kpc
	dis_list_ucd.append(distance_ucd)

bs_ucd = 12  # bin size for ucd (kpc)
count_ucd_binned = scipy.stats.binned_statistic(dis_list_ucd,dis_list_ucd,statistic='count',bins=radius_max/bs_ucd,range=(0.,radius_max))
bin_mean_ucd = scipy.stats.binned_statistic(dis_list_ucd,dis_list_ucd,statistic='mean',bins=radius_max/bs_ucd,range=(0.,radius_max))[0] #mean value for distance in each bin
bin_left_edges_ucd = np.array(count_ucd_binned[1][:-1],dtype='f8')  #0,6,12 ...
bin_mean_ucd[np.isnan(bin_mean_ucd)] = bin_left_edges_ucd[np.isnan(bin_mean_ucd)]+bs_ucd/2   # fill the NaN entries with middle value of each bin
ucd_counts = np.array(count_ucd_binned[0],dtype='f8')                    #counts
ucd_counts_err  = np.sqrt(ucd_counts)                                             #noise (assume poisson)
areas_ucd = ((bin_left_edges_ucd+bs_ucd)**2 - bin_left_edges_ucd**2)*pi
ucd_density = ucd_counts/areas_ucd
ucd_density_err = ucd_counts_err/areas_ucd

'''==================fitting=========================='''
#fit range
fit_min = 25   # kpc
fit_max = 210  # kpc

#slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))
slope, intercept, r_value = gc_fitting(gc_density, r_maj[:-1], fit_min, fit_max)
abline_values = [radius**slope*exp(intercept) for radius in np.arange(1,1000)] # parameter for drawing the fitted line
print 'gc slope:',slope,'gc intercept:',intercept,'r-squared:', r_value**2

#mask the bins that contains no ucd
ucd_density_fit = ucd_density[fit_min/bs_ucd:fit_max/bs_ucd] # fit range [15,190]
bin_mean_ucd_fit = bin_mean_ucd[fit_min/bs_ucd:fit_max/bs_ucd] # fit range [15,190]
mask  = ucd_density_fit != 0.0   #mask selection
ucd_density_fit = ucd_density_fit[mask]
bin_mean_ucd_fit = bin_mean_ucd_fit[mask]

#linear fit for ucd surface density
slope_ucd, intercept_ucd, r_value_ucd, p_value_ucd, std_err_ucd = scipy.stats.linregress(np.log(bin_mean_ucd_fit),np.log(ucd_density_fit))
abline_values_ucd = [radius**slope_ucd*exp(intercept_ucd) for radius in np.arange(1,1000)]

'''===================plot============================'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(r_maj[:-1],gc_density,yerr=gc_density_err,fmt='.k')
ax.errorbar(bin_mean_ucd,ucd_density,yerr=ucd_density_err,fmt='.g')
ax.plot(np.arange(1,1000), abline_values, 'b',label = r'$n_{\rm GC}=$'+str(round(slope,3)))
ax.plot(np.arange(1,1000), abline_values_ucd, 'r',label = r'$n_{\rm UCD}=$'+str(round(slope_ucd,3)))
ax.axvline(x=fit_min, ymin=0, ymax = 10, linewidth=1, linestyle='dashed', color='k')
ax.axvline(x=fit_max, ymin=0, ymax = 10, linewidth=1, linestyle='dashed',color='k')
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlabel(r'$R_{\rm maj} {\rm [kpc]}$',fontsize=15)
ax.set_ylabel(r'Surface density [N/kpc$^2$]')
ax.set_xlim([1,1000])
ax.set_ylim([0.0001,10])
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(numpoints=1,frameon=False,loc='lower left')
plt.title('fit range:'+str(fit_min)+'~'+str(fit_max)+', e= '+str(e))
plt.savefig('pics/gc.e.around.M87.png')
plt.show()