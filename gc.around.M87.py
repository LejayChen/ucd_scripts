from astropy.table import Table
from math import *
import matplotlib.pyplot as plt 
import scipy.stats
import numpy as np

cat = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat = cat[cat['p_gc']>0.95] # sort out the GCs
cat_ucd = Table.read('NGVS.pilot.92ucds.fits')
M87 = (187.70583, 12.39111)
DM = 22.2*1000 # in Kpc
radius_max = 420. #kpc   max radius for bin statistics (gc and ucd)

# Bin Statistics for gc
dis_list = []
for i in range(len(cat)):
	ra = cat[i]['ra']
	dec = cat[i]['dec']
	distance = sqrt((ra-M87[0])**2+(dec-M87[1])**2)  # in degree
	#distance = acos(sin(dec)*sin(M87[1])+cos(dec)*cos(M87[1])*cos(abs(ra-M87[0]))) # in degree
	distance = distance/180.*pi*DM
	dis_list.append(distance)

#bins_gc = bin edges, set larger bin size for larger distance from M87 
bins_gc = np.arange(0.,400.,1.)**1.43
bins_gc = bins_gc[bins_gc<radius_max]
count_gc_binned = scipy.stats.binned_statistic(dis_list,dis_list,statistic='count',bins=bins_gc,range=(0.,radius_max))
bin_mean_gc = scipy.stats.binned_statistic(dis_list,dis_list,statistic='mean',bins=bins_gc,range=(0.,radius_max))[0] #mean value for distance in each bin
bin_mean_gc[0]=0
gc_counts = np.array(count_gc_binned[0],dtype='f8')
gc_counts_err = np.sqrt(gc_counts)
bin_left_edges = np.array(count_gc_binned[1][:-1],dtype='f8')  #0,2,4,6 ...
areas = np.array([],dtype='f4')
for i in range(len(bins_gc[:-1])):
	areas = np.append(areas, (bins_gc[i+1]**2 - bins_gc[i]**2)*pi)
gc_density = gc_counts/areas
gc_density_err = gc_counts_err/areas

# Bin Statistics for ucd
dis_list_ucd = []
for i in range(len(cat_ucd)):
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	distance_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)  # in degree
	#distance = acos(sin(dec)*sin(M87[1])+cos(dec)*cos(M87[1])*cos(abs(ra-M87[0]))) # in degree
	distance_ucd = distance_ucd/180.*pi*DM
	dis_list_ucd.append(distance_ucd)


bs_ucd = 6  # bin size for ucd (kpc)
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
fit_min = 20   # kpc
fit_max = 250  # kpc

#linear fit for gc surface density, mask data out of fitting range
gc_density_fit = gc_density[bin_mean_gc>20.]
bin_mean_gc_fit = bin_mean_gc[bin_mean_gc>20.]
gc_density_fit = gc_density_fit[bin_mean_gc_fit<250.]
bin_mean_gc_fit = bin_mean_gc_fit[bin_mean_gc_fit<250.]
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))
abline_values = [radius**slope*exp(intercept) for radius in np.arange(1,1000)]
# print slope,intercept

#mask the bins that contains no ucd
ucd_density_fit = ucd_density[fit_min/bs_ucd:fit_max/bs_ucd] # fit range [20.,250.]
bin_mean_ucd_fit = bin_mean_ucd[fit_min/bs_ucd:fit_max/bs_ucd]
mask  = ucd_density_fit != 0.0
ucd_density_fit = ucd_density_fit[mask]
bin_mean_ucd_fit = bin_mean_ucd_fit[mask]
#linear fit for ucd surface density
slope_ucd, intercept_ucd, r_value_ucd, p_value_ucd, std_err_ucd = scipy.stats.linregress(np.log(bin_mean_ucd_fit),np.log(ucd_density_fit))
abline_values_ucd = [radius**slope_ucd*exp(intercept_ucd) for radius in np.arange(1,1000)]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(bin_mean_gc,gc_density,yerr=gc_density_err,fmt='.k')
ax.errorbar(bin_mean_ucd,ucd_density,yerr=ucd_density_err,fmt='.g')
ax.plot(np.arange(1,1000), abline_values, 'b',label = r'$n_{GC}=$'+str(round(slope,2)))
ax.plot(np.arange(1,1000), abline_values_ucd, 'r',label = r'$n_{UCD}=$'+str(round(slope_ucd,2)))
ax.axvline(x=fit_min, ymin=0, ymax = 10, linewidth=1, linestyle='dashed', color='k')
ax.axvline(x=fit_max, ymin=0, ymax = 10, linewidth=1, linestyle='dashed',color='k')
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlabel('Radius [kpc]')
ax.set_ylabel(r'Surface density [N/kpc$^2$]')
ax.set_xlim([1,500])
ax.set_ylim([0.0001,10])
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(numpoints=1,frameon=False,loc='lower left')
plt.show()