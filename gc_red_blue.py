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

def mask(bins_gc, gc_counts, dis_list, theta_list):

	def slice_list(dis_list,theta_list,i):
		theta_list_slice = theta_list[dis_list>float(bins_gc[i])]
		dis_list_slice = dis_list[dis_list>float(bins_gc[i])]

		theta_list_slice = theta_list_slice[dis_list_slice<float(bins_gc[i+1])]	
		return theta_list_slice

	cat_ucd = Table.read('NGVS.pilot.92ucds.fits')
	cat_dEN = Table.read('pp.gal.nuc2.s.master.new.fits')

	dis_ucd_list = np.sqrt((cat_ucd['RA']-M87[0])**2+(cat_ucd['DEC']-M87[1])**2) /180.*pi*DIS
	dis_dEN_list = np.sqrt((cat_dEN['RA']-M87[0])**2+(cat_dEN['DEC']-M87[1])**2) /180.*pi*DIS

	theta_ucd_list = np.arctan((cat_ucd['RA']-M87[0])/(cat_ucd['DEC']-M87[1]))
	theta_dEN_list = np.arctan((cat_dEN['RA']-M87[0])/(cat_dEN['DEC']-M87[1]))

	gc_densities = []
	areas = []
	for i in range(len(bins_gc[:-1])):
		theta_list_slice = slice_list(dis_list, theta_list,i)
		theta_ucd_list_slice = slice_list(dis_ucd_list, theta_ucd_list,i)
		theta_dEN_list_slice = slice_list(dis_dEN_list, theta_dEN_list,i)

		theta_list_slice_byucd = []
		theta_list_slice_bydEN = []
		for theta in theta_ucd_list_slice:
			theta_list_slice_byucd = theta_list_slice[theta_list_slice<theta+0.05*pi]
			theta_list_slice_byucd = theta_list_slice_byucd[theta_list_slice_byucd>theta-0.05*pi]
			gc_counts[i] = gc_counts[i] - len(theta_list_slice_byucd)

		for theta in theta_dEN_list_slice:
			theta_list_slice_bydEN = theta_list_slice[theta_list_slice<theta+0.05*pi]
			theta_list_slice_bydEN = theta_list_slice_bydEN[theta_list_slice_bydEN>theta-0.05*pi]
			gc_counts[i] = gc_counts[i] - len(theta_list_slice_bydEN)

		area = (bins_gc[i+1]**2 - bins_gc[i]**2)*pi
		area = area*(1- 0.1/2*(len(theta_list_slice_byucd)+len(theta_list_slice_bydEN)))
		areas.append(area)
		gc_densities.append(gc_counts[i]/area)

	return np.array(gc_densities),np.array(areas)

def gc_density_stat(cat_GC):
	dis_list = []
	theta_list = []
	for i in range(len(cat_GC)):
		ra = cat_GC[i]['ra']
		dec = cat_GC[i]['dec']
		distance = sqrt((ra-M87[0])**2+(dec-M87[1])**2) /180.*pi*DIS # in kpc
		theta = np.arctan((ra-M87[0])/(dec-M87[1]))
		dis_list.append(distance)
		theta_list.append(theta)
	bins_gc = np.exp(np.arange(0.955,7.,0.05)) #bins_gc = bin edges, set larger bin size for larger distance from M87 
	bins_gc = bins_gc[bins_gc<radius_max]        #exclued outer part that GC distribution is flattened out

	count_gc_binned = scipy.stats.binned_statistic(dis_list,dis_list,statistic='count',bins=bins_gc,range=(0.,radius_max))
	bin_mean_gc = scipy.stats.binned_statistic(dis_list,dis_list,statistic='mean',bins=bins_gc,range=(0.,radius_max))[0] #mean value for distance in each bin

	gc_counts = np.array(count_gc_binned[0],dtype='f8')
	gc_counts_err = np.sqrt(gc_counts)
		
	gc_density, areas = mask(bins_gc, gc_counts, dis_list, theta_list)
	gc_density_err = gc_counts_err/areas

	return gc_density, gc_density_err,bin_mean_gc

def gc_fitting(gc_density, bin_mean_gc, fit_min, fit_max):

	gc_density = gc_density[~np.isnan(bin_mean_gc)]
	bin_mean_gc = bin_mean_gc[~np.isnan(bin_mean_gc)]

	gc_density_fit = gc_density[bin_mean_gc>fit_min]
	bin_mean_gc_fit = bin_mean_gc[bin_mean_gc>fit_min]

	gc_density_fit = gc_density_fit[bin_mean_gc_fit<fit_max]
	bin_mean_gc_fit = bin_mean_gc_fit[bin_mean_gc_fit<fit_max]

	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))

	return slope, intercept, r_value
	

#Binned Statistics
gc_density_b, gc_density_err_b, bin_mean_gc_b = gc_density_stat(cat_gc_b)
gc_density_r, gc_density_err_r, bin_mean_gc_r = gc_density_stat(cat_gc_r)

#r_maj = np.exp(np.arange(0.955,7.,0.05))
#r_maj = r_maj[r_maj<radius_max] 
#e = 0.5
#gc_density_b, gc_density_err_b = elliptical_stat(cat_gc_b, r_maj,e=e)
#gc_density_r, gc_density_err_r = elliptical_stat(cat_gc_r, r_maj,e=e)

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