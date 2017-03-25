import scipy.stats
import numpy as np
from math import *
from astropy.table import Table

M87 = (187.70583, 12.39111)
PA_M87 = 149./180*pi
DIS = 17.21*1000 # in kpc (luminosity distance from NED)
radius_max = 500. #kpc   max radius for bin statistics (gc and ucd)
ell_default = 0.5

'''
    this script is for estimating the background spatial distribution of objects around M87
    1D circle rings and elliptical rings statistics
    circle rings funcs:
    	mask_c()
    	gc_density_stat()
    elliptical rings funcs:
              mask_e()
              elliptical_stat()
    fitting func: gc_fitting()
'''
def load_ucd_dEN_cat():
	'''load the catalogs (ra, dec info) for mask'''
	cat_ucd = Table.read('NGVS.pilot.92ucds.fits')
	cat_dEN = Table.read('pp.gal.nuc2.s.master.new.fits')
	ra_ucd_list = cat_ucd['RA']
	dec_ucd_list = cat_ucd['DEC']
	ra_dEN_list = cat_dEN['RA']
	dec_dEN_list = cat_dEN['DEC']
	return ra_ucd_list, dec_ucd_list, ra_dEN_list, dec_dEN_list

def mask_c(bins_gc, gc_counts, dis_list, theta_list):
	'''mask the UCDs and dENs in circle shape ring stats'''

	def dis_slice_list(dis_list,theta_list,index):
		'''slice the UCD and dEN list between certain rings'''

		theta_list_slice = theta_list[dis_list>bins_gc[index]]
		dis_list_slice = dis_list[dis_list>bins_gc[index]]
		theta_list_slice = theta_list_slice[dis_list_slice<bins_gc[index+1]]
		dis_list_slice = dis_list_slice[dis_list_slice<bins_gc[index+1]]
		return dis_list_slice, theta_list_slice

	#load RA,DEC data from UCD and dE,N catalogs 
	ra_ucd_list, dec_ucd_list, ra_dEN_list, dec_dEN_list = load_ucd_dEN_cat()
	dis_ucd_list = np.sqrt((ra_ucd_list-M87[0])**2+(dec_ucd_list-M87[1])**2) /180.*pi*DIS
	dis_dEN_list = np.sqrt((ra_dEN_list-M87[0])**2+(dec_dEN_list - M87[1])**2) /180.*pi*DIS
	theta_ucd_list = np.arctan((ra_ucd_list - M87[0]) / (dec_ucd_list - M87[1]))
	theta_dEN_list = np.arctan((ra_dEN_list - M87[0]) / (dec_dEN_list - M87[1]))	

	gc_counts_cor = []
	areas = []  
	for i in range(len(bins_gc[:-1])):
		#print '=======',bins_gc[i],'============='
		dis_list_slice, theta_list_slice = dis_slice_list(dis_list, theta_list,i)
		dis_ucd_list_slice, theta_ucd_list_slice = dis_slice_list(dis_ucd_list, theta_ucd_list,i)
		dis_dEN_list_slice, theta_dEN_list_slice = dis_slice_list(dis_dEN_list, theta_dEN_list,i)

		theta_list_slice_byucd = []
		theta_list_slice_bydEN = []
		for k in range(len(theta_ucd_list_slice)):
			theta = theta_ucd_list_slice[k]
			theta_list_slice_byucd = theta_list_slice[abs(dis_list_slice - dis_ucd_list_slice[k]) < 2]
			theta_list_slice_byucd = theta_list_slice_byucd[abs(theta_list_slice_byucd - theta) < 0.01*pi]
			gc_counts[i] = gc_counts[i] - len(theta_list_slice_byucd)
			#print 'GC cut by ucd',k+1,':', len(theta_list_slice_byucd)

		for k in range(len(theta_dEN_list_slice)):
			theta = theta_dEN_list_slice[k]
			theta_list_slice_bydEN = theta_list_slice[abs(dis_list_slice - dis_dEN_list_slice[k]) < 4]				
			theta_list_slice_bydEN = theta_list_slice_bydEN[abs(theta_list_slice_bydEN - theta) < 0.01*pi]
			gc_counts[i] = gc_counts[i] - len(theta_list_slice_bydEN)	
			#print 'GC cut by dEN',k+1,':', len(theta_list_slice_bydEN)		

		area = (bins_gc[i+1]**2 - bins_gc[i]**2)*pi
		area = area*(1- 0.02/2*(len(theta_list_slice_byucd)+len(theta_list_slice_bydEN)))
		areas.append(area)
		gc_counts_cor.append(gc_counts[i])

	return np.array(gc_counts_cor), np.array(areas)
	
def gc_fitting(gc_density, bin_mean_gc, fit_min, fit_max):
	'''1D power law fitting (linear in log-log scale) '''

	gc_density = gc_density[~np.isnan(bin_mean_gc)]
	bin_mean_gc = bin_mean_gc[~np.isnan(bin_mean_gc)]

	gc_density_fit = gc_density[bin_mean_gc>fit_min]
	bin_mean_gc_fit = bin_mean_gc[bin_mean_gc>fit_min]

	gc_density_fit = gc_density_fit[bin_mean_gc_fit<fit_max]
	bin_mean_gc_fit = bin_mean_gc_fit[bin_mean_gc_fit<fit_max]

	slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(bin_mean_gc_fit),np.log(gc_density_fit))

	return slope, intercept, r_value

def gc_density_stat(cat_GC, mask='on'):
	dis_list = []
	theta_list = []
	for i in range(len(cat_GC)):
		ra = cat_GC[i]['ra']
		dec = cat_GC[i]['dec']
		distance = np.sqrt((ra-M87[0])**2+(dec-M87[1])**2) /180.*pi*DIS # in kpc
		theta = np.arctan((ra-M87[0])/(dec-M87[1]))
		dis_list.append(distance)
		theta_list.append(theta)

	bins_gc = np.exp(np.arange(0.955,7.,0.05)) #bins_gc = bin edges, set larger bin size for larger distance from M87 
	bins_gc = bins_gc[bins_gc<radius_max]        #exclued outer part that GC distribution is flattened out

	gc_counts = scipy.stats.binned_statistic(dis_list, dis_list,statistic='count',bins=bins_gc,range=(0.,radius_max))[0]
	bin_mean_gc = scipy.stats.binned_statistic(dis_list, dis_list,statistic='mean',bins=bins_gc,range=(0.,radius_max))[0] #mean value for distance in each bin

	#gc_counts = np.array(count_gc_binned[0],dtype='f8')
	if mask=='on': 
		gc_counts, areas = mask_c(bins_gc, gc_counts, np.array(dis_list), np.array(theta_list)) 
	elif mask=='off':
		areas = np.array([])
		for i in range(len(bins_gc[:-1])):
			areas = np.append(areas, (bins_gc[i+1]**2 - bins_gc[i]**2)*pi)
	else:
		raise KeyError('Argument mask only have options "on" and "off" ')

	gc_density = gc_counts/areas
	gc_counts_err = np.sqrt(gc_counts)
	gc_density_err = gc_counts_err/areas	

	return gc_density, gc_density_err, bin_mean_gc

def  which_ring(ra,dec,r_maj, PA=PA_M87, ell=ell_default):
	ra_dis = (ra - M87[0])/180.*pi*DIS
	dec_dis = (dec - M87[1])/180.*pi*DIS
	dis = sqrt(ra_dis**2+dec_dis**2)
	theta = np.arctan(dec_dis/ra_dis) - PA

	for j in range(len(r_maj)):   # for ring in elliptical rings
		if dis > r_maj[j]:
			continue
		else:
			if dis>sqrt((r_maj[j]*cos(theta))**2+(ell*r_maj[j]*sin(theta))**2):
				continue
			else:
				return j-1
	return None #doesn't belong to anly ring

def mask_e(r_maj, gc_counts, areas, gc_ra_list, gc_dec_list, mask_r=2.0, ell=ell_default):
	'''mask the UCDs and dENs in elliptical shape ring stats'''

	def count_gc_mask_r(ra, dec, gc_ra_list, gc_dec_list, index):
		'''count the GCs to be substracted from bkg estimation and recalculate the area'''
		mask = np.sqrt((gc_ra_list - ra)**2 + (gc_dec_list - dec)**2)/180.*pi*DIS<mask_r
		gc_ra_list = gc_ra_list[mask]
		gc_dec_list = gc_dec_list[mask]
		gc_in_mask_r = len(gc_ra_list)

		num_gc_out = num_gc_in = num_gc_mid =  0 # num_gc_out: number of GCs in index+1 ring, num_gc_in: ... in index-1 ring
		for i in range(len(gc_ra_list)):
			ra_gc = gc_ra_list[i]
			dec_gc = gc_dec_list[i]
			index_gc = which_ring(ra_gc, dec_gc, r_maj,ell=ell)
			if index_gc == index+1 or which_ring(ra_gc, dec_gc, r_maj)==None:  
				num_gc_out += 1
			elif index_gc == index-1:
				num_gc_in += 1
			elif index_gc == index:
				num_gc_mid +=1

		return gc_in_mask_r, num_gc_mid, num_gc_in, num_gc_out

	ra_ucd_list, dec_ucd_list, ra_dEN_list, dec_dEN_list = load_ucd_dEN_cat()
	for i in range(len(ra_ucd_list)):
		ra = ra_ucd_list[i]
		dec = dec_ucd_list[i]
		index = which_ring(ra,dec,r_maj)
		if index == None: continue

		gc_in_mask_r, num_gc_mid, num_gc_in, num_gc_out = count_gc_mask_r(ra,dec,np.copy(gc_ra_list), np.copy(gc_dec_list),index)
		if gc_in_mask_r>0:
			areas[index] = areas[index] - (pi*mask_r**2)*ell*num_gc_mid/float(gc_in_mask_r)
			areas[index+1] = areas[index+1] - (pi*mask_r**2)*ell*num_gc_out/float(gc_in_mask_r)
			areas[index-1] = areas[index-1] - (pi*mask_r**2)*ell*num_gc_in/float(gc_in_mask_r)

		gc_counts[index] = gc_counts[index] - num_gc_mid
		gc_counts[index+1] = gc_counts[index+1] - num_gc_out
		gc_counts[index-1] = gc_counts[index-1] - num_gc_in

	return np.array(gc_counts),np.array(areas)

def elliptical_stat(cat, r_maj, PA=PA_M87,ell=ell_default, mask='on'):
	r_min = ell*r_maj
	gc_counts = np.zeros(len(r_maj) -1) 
	gc_ra_list = cat['ra']
	gc_dec_list  =cat['dec']

	'''binned stat in elliptical rings (bin edges = r_maj)'''
	for i in range(len(cat)):   # for GC in GC catalog
		ra_dis = (gc_ra_list[i] - M87[0])/180.*pi*DIS
		dec_dis = (gc_dec_list[i] - M87[1])/180.*pi*DIS
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
	
	areas = np.array([],dtype='f4')
	for i in range(len(r_maj[:-1])):
		areas = np.append(areas, pi*ell*(r_maj[i+1]**2-r_maj[i]**2))

	if mask=='on': gc_counts, areas = mask_e(r_maj, gc_counts, areas, gc_ra_list, gc_dec_list, mask_r=5.0, ell=ell)

	gc_density = gc_counts/areas 
	gc_counts_err = np.sqrt(gc_counts)
	gc_density_err = gc_counts_err/areas

	return gc_density, gc_density_err

def  cal_r_maj(ra,dec, PA= PA_M87, ell=ell_default):
	'''calculate corresponding r_maj for expected GC density around UCD '''

	ra_dis = (ra - M87[0])/180.*pi*DIS
	dec_dis = (dec - M87[1])/180.*pi*DIS
	dis = sqrt(ra_dis**2+dec_dis**2)

	theta = np.arctan(dec_dis/ra_dis) - PA
	r_maj_ucd = dis/sqrt(cos(theta)**2+(ell*sin(theta))**2)

	return r_maj_ucd