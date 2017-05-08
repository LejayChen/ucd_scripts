from astropy.table import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import random

from gc_stat import *
from robust_mean import *
from func import *

DIS = 17.21*1000 #kpc
M87 = (187.70583, 12.39111)
fit_min = 25
fit_max = 210
bkg_level_radius  = 250
step = 0.25  # in kpc
start = 0.5
radius = np.arange(start, start+4, step)   # in kpc 
cat_ucd = Table.read('new.NGVS.pilot.92ucds.fits')  #catalog of UCDs

cat_gc = Table.read('gc_UCD_dEN_excluded.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra_gc = np.array(cat_gc['ra'])      #RA of GCs
dec_gc = np.array(cat_gc['dec'])  #Dec of GCs

C = np.array([])   # C=clustering signal     (measured/bkg)
gc_counts = np.array([])  # count of GCs 
dis_ucd_gc_list = np.array([])  # list of each gc's distance to UCD

def bkg_1d(dis_M87_ucd,slope, intercept):
	exp_density = exp(intercept)*dis_M87_ucd**(slope)  # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(intercept)*bkg_level_radius**(slope)    # uniform bakground in the field (possibly non-GC objects)
	return exp_density,bkg_unif

def gc_stat(r, dis_ucd_gc, exp_density, bkg_unif):
	'''count GC around each UCD in differently sized bins'''

	dis_ucd_gc_mask = np.array([])  #mask out the area for counting  (ring shaped mask from r to r+step)
	dis_ucd_gc_mask = dis_ucd_gc[dis_ucd_gc<(r+step)] #outer bound
	dis_ucd_gc_mask = dis_ucd_gc_mask[dis_ucd_gc_mask>r]  #inner bound
	area = pi*((r+step)**2 - r**2)  # area of that ring  (in kpc^2)	

	# expected count and measued count of GCs
	gc_expected = (exp_density - bkg_unif)*area
	gc_count = len(dis_ucd_gc_mask)
	gc_count_cor = gc_count - bkg_unif*area
	return gc_count, gc_count_cor, gc_expected, dis_ucd_gc_mask

def C_ave(gc_counts, dis_ucd_gc_list, C):
	dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(dis_ucd_gc_list)/len(radius), len(radius))   #shape: no_ucd*len(radius)
	dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0) 

	C = C.reshape(len(C)/len(radius),len(radius))  #original clustering signals from individual UCDs
	gc_counts = gc_counts.reshape(len(gc_counts)/len(radius),len(radius)) 
	gc_counts_mean = np.nanmean(gc_counts,0) 

	C_mean = np.array([])
	C_std = np.array([])
	for i in range(len(radius)):
		C_cloumn, mean, std,robust_mean_log = robust_mean(C[:,i],iter=0,show_step=True)
		C_mean = np.append(C_mean,mean)   #robust averge C signal over all UCDs
		C_std = np.append(C_std,std)  #robust averge C signal over all UCDs

	return C_mean,C_std,dis_ucd_gc_list_mean,gc_counts_mean

def mean_dis(dis_ucd_gc_list, gc_count, dis_ucd_gc_mask):
	'''GC mean distance to host UCD for individual UCD'''

	if gc_count == 0:
		dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
	else:
		dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))	
	return dis_ucd_gc_list

def sample_selection(cat_gc):
	names = ['ra','dec','dis']
	dtypes = ['f8','f8','f8']
	cat_ucd_dEN = Table(names=(names),dtype=(dtypes))

	for i in range(len(cat_ucd)):
		dis_ucd_M87 = np.sqrt((cat_ucd[i]['RA'] - M87[0])**2+(cat_ucd[i]['DEC'] - M87[1])**2)/180.*pi*DIS
		mask = abs(np.sqrt((cat_gc['ra'] - M87[0])**2+(cat_gc['dec'] - M87[1])**2)/180.*pi*DIS - dis_ucd_M87) < 0.5
		cat_gc_slice = cat_gc[mask]

		index = random.randint(1,len(cat_gc_slice)-1)
		ra = cat_gc_slice[index]['ra']
		dec = cat_gc_slice[index]['dec']
		dis = np.sqrt((ra - M87[0])**2+(dec - M87[1])**2)/180.*pi*DIS
	
		cat_ucd_dEN.add_row([ra,dec,dis])

	cat_ucd_dEN = cat_ucd_dEN[cat_ucd_dEN['dis']>fit_min]
	cat_ucd_dEN = cat_ucd_dEN[cat_ucd_dEN['dis']<fit_max]

	#cat_ucd_dEN = cat_ucd_dEN[:65]
	return cat_ucd_dEN

for test_run in range(50):
	print 'GC autoclustering test run No.'+str(test_run+1)
	cat_ucd_dEN = sample_selection(cat_gc)
	gc_count_total = 0
	for i in range(len(cat_ucd_dEN)): 
		ra_ucd = cat_ucd_dEN[i]['ra']
		dec_ucd = cat_ucd_dEN[i]['dec']
		dis_M87_ucd = cat_ucd_dEN[i]['dis']

		exp_density, bkg_unif = bkg_1d(dis_M87_ucd,slope = -1.88364525707, intercept =  5.20946384186)
		dis_ucd_gc = np.sqrt((ra_gc - ra_ucd)**2 + (dec_gc - dec_ucd)**2)/180.*pi*DIS  # in kpc (list)
		for r in radius:
			gc_count, gc_count_cor, gc_expected, dis_ucd_gc_mask = gc_stat(r, dis_ucd_gc, exp_density, bkg_unif)
			gc_count_total += gc_count
			C = np.append(C, gc_count_cor/gc_expected) 

			gc_counts = np.append(gc_counts,gc_count)
			dis_ucd_gc_list = mean_dis(dis_ucd_gc_list, gc_count, dis_ucd_gc_mask)

C_mean, C_std, dis_ucd_gc_list_mean, gc_counts_mean = C_ave(gc_counts, dis_ucd_gc_list, C)
#==============print statistics on the screen=================
np.set_printoptions(precision=2)
no_ucd = len(cat_ucd_dEN)
print 'Number of UCDs:', no_ucd
print 'Number of GCs:',gc_count_total
print 'Radius bins:',radius
print '==============================='
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean C:',C_mean
print 'GC counts Mean:',gc_counts_mean
print 'Std:',C_std/sqrt(no_ucd) 

#============plot the final figure=================================================
fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/sqrt(no_ucd),fmt='ko',label = 'All GCs')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]',fontsize=14)
ax.set_ylabel(r'C=<$\Sigma_{\rm measured}/\Sigma_{\rm expected}$>',fontsize=16)
ax.get_xaxis().set_tick_params(direction='in', which='major', width=1)
ax.get_xaxis().set_tick_params(direction='in', which='minor', width=1)
ax.get_yaxis().set_tick_params(direction='in', which='major', width=1)
ax.get_yaxis().set_tick_params(direction='in', which='minor', width=1)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0.25,4.5)
ax.set_ylim(0,int(C_mean.max())+2)
ax.legend(numpoints=1,frameon=False,loc='upper right')
plt.savefig('pics/clustering.test.png')

plt.show()