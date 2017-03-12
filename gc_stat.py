import scipy.stats
import numpy as np
from math import *
from astropy.table import Table
np.set_printoptions(precision=2)

def gc_fitting(gc_density, bin_mean_gc, fit_min, fit_max):

	gc_density = gc_density[~np.isnan(bin_mean_gc)]
	bin_mean_gc = bin_mean_gc[~np.isnan(bin_mean_gc)]

	gc_density_fit = gc_density[bin_mean_gc>fit_min]
	bin_mean_gc_fit = bin_mean_gc[bin_mean_gc>fit_min]

	gc_density_fit = gc_density_fit[bin_mean_gc_fit<fit_max]
	bin_mean_gc_fit = bin_mean_gc_fit[bin_mean_gc_fit<fit_max]

	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))

	return slope, intercept, r_value

def gc_density_stat(cat_GC):
	dis_list = []
	for i in range(len(cat_GC)):
		ra = cat_GC[i]['ra']
		dec = cat_GC[i]['dec']
		distance = sqrt((ra-M87[0])**2+(dec-M87[1])**2) /180.*pi*DIS # in kpc
		dis_list.append(distance)
	bins_gc = np.exp(np.arange(0.955,7.,0.05)) #bins_gc = bin edges, set larger bin size for larger distance from M87 
	bins_gc = bins_gc[bins_gc<radius_max]        #exclued outer part that GC distribution is flattened out

	count_gc_binned = scipy.stats.binned_statistic(dis_list,dis_list,statistic='count',bins=bins_gc,range=(0.,radius_max))
	bin_mean_gc = scipy.stats.binned_statistic(dis_list,dis_list,statistic='mean',bins=bins_gc,range=(0.,radius_max))[0] #mean value for distance in each bin

	gc_counts = np.array(count_gc_binned[0],dtype='f8')
	gc_counts_err = np.sqrt(gc_counts)
	areas = np.array([],dtype='f4')
	for i in range(len(bins_gc[:-1])):
		areas = np.append(areas, (bins_gc[i+1]**2 - bins_gc[i]**2)*pi)
	gc_density = gc_counts/areas
	gc_density_err = gc_counts_err/areas

	return gc_density, gc_density_err,bin_mean_gc

def elliptical_stat(cat, r_maj, PA=149./180*pi,e=0.5):
	M87 = (187.70583, 12.39111)
	DIS = 17.21*1000
	r_min = e*r_maj
	gc_counts = np.zeros(len(r_maj) -1) 

	for i in range(len(cat)):   # for GC in GC catalog
		ra_dis = (cat[i]['ra'] - M87[0])/180.*pi*DIS
		dec_dis = (cat[i]['dec'] - M87[1])/180.*pi*DIS
		dis = sqrt(ra_dis**2+dec_dis**2)
		theta = np.arctan(dec_dis/ra_dis) - PA
		for j in range(len(r_maj)):   # for ring in elliptical rings
			if dis > r_maj[j]:
				continue
			else:
				if dis>sqrt((r_maj[j]*cos(theta))**2+(r_min[j]*sin(theta))**2):
					continue
				else:
					if j>=1: gc_counts[j-1] += 1
					break

	gc_counts_err = np.sqrt(gc_counts)
	areas = np.array([],dtype='f4')
	for i in range(len(r_maj[:-1])):
		areas = np.append(areas, pi*e*(r_maj[i+1]**2-r_maj[i]**2))
	gc_density = gc_counts/areas
	gc_density_err = gc_counts_err/areas

	return gc_density, gc_density_err

def  which_ring(ra,dec,r_maj, PA=149./180*pi, e=0.5):
	M87 = (187.70583, 12.39111)
	DIS = 17.21*1000
	r_min = e*r_maj
	ra_dis = (ra - M87[0])/180.*pi*DIS
	dec_dis = (dec - M87[1])/180.*pi*DIS
	dis = sqrt(ra_dis**2+dec_dis**2)
	theta = np.arctan(dec_dis/ra_dis) - PA

	gc_counts = np.zeros(len(r_maj) -1) 
	for j in range(len(r_maj)):   # for ring in elliptical rings
		if dis > r_maj[j]:
			continue
		else:
			if dis>sqrt((r_maj[j]*cos(theta))**2+(r_min[j]*sin(theta))**2):
				continue
			else:
				if j>=1: gc_counts[j-1] += 1
				return j-1

def bkg_e(ra,dec,r_maj,intercept,slope):
	exp_density = exp(intercept)*r_maj[which_ring(ra,dec,r_maj)]**(slope)  # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(intercept)*230**(slope)    # uniform bakground in the field (possibly non-GC objects)	
	return exp_density,bkg_unif
