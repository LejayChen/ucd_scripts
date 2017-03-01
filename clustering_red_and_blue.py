from astropy.table import *
from math import *
from robust_mean import *
from func import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

DIS = 17.21*1000 #kpc
M87 = (187.70583, 12.39111)
logfile = open('clustering_red_and_blue.log','w')
logfile.write(str(datetime.now())+'\n\n')
logfile.write('r, gc_count, gc_count_b, gc_count_r, gc_expected, C \n\n')

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra_gc = np.array(cat_gc['ra'])      #RA of GCs
dec_gc = np.array(cat_gc['dec'])  #Dec of GCs
#cat_gc = cat_gc[cat_gc['gmag']<24] 

cat_gc_b = cat_gc[cat_gc['gmag'] - cat_gc['imag']<0.87] 
ra_gc_b = np.array(cat_gc_b['ra'])      #RA of GCs
dec_gc_b = np.array(cat_gc_b['dec'])  #Dec of GCs

cat_gc_r = cat_gc[cat_gc['gmag'] - cat_gc['imag']>0.87] 
ra_gc_r = np.array(cat_gc_r['ra'])      #RA of GCs
dec_gc_r = np.array(cat_gc_r['dec'])  #Dec of GCs

step = 0.5  # in kpc
radius = np.arange(0.5,4.5,step)   # in kpc  from 0.25 to 3.75 step 0.5
str_radius = str(radius).strip('[]').split()  #for table name
tab_names = ['ID','ra','dec','DIS','r_h','g-i']
for str_r in str_radius:
	tab_names.append(str_r+'_count')
	tab_names.append(str_r+'_count_cor')
	tab_names.append(str_r+'_C')
tab_dtype = ['a8','f8','f8','f4','f4','f4']+['i2','f2','f4']*len(radius)
ucd_data = Table(names=(tab_names),dtype=(tab_dtype))

gc_count_total = 0
C = np.array([])   # C=clustering signal     (measured/bkg)
gc_counts = np.array([])  # count of GCs 
dis_ucd_gc_list = np.array([])  # list of each gc's distance to UCD

gc_count_total = 0
C_b = np.array([])   
gc_counts_b = np.array([])  
dis_ucd_gc_list_b = np.array([])  

C_r = np.array([])   
gc_counts_r = np.array([])  
dis_ucd_gc_list_r = np.array([]) 

def bkg_2d(ra,dec):
	dis_ra = (ra - M87[0])/180.*pi*DIS  
	dis_dec = (dec - M87[1])/180.*pi*DIS  
	n = -1.5716456
	PA =  0.658513273
	e = 0.740146285
	H = 85.1889069
	Const = -0.0123046213

	exp_density = func((dis_ra,dis_dec), n, PA, e, H, Const)
	bkg_unif = 0
	return exp_density,bkg_unif

def bkg_1d(dis_M87_ucd,slope, intercept):
	exp_density = exp(intercept)*dis_M87_ucd**(slope)  # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(intercept)*200**(slope)    # uniform bakground in the field (possibly non-GC objects)	
	return exp_density,bkg_unif

def gc_stat(dis_ucd_gc, exp_density, bkg_unif):
	'''
	count GC around each UCD
	'''
	dis_ucd_gc_mask = np.array([])
	#mask out the area for counting  (ring shaped mask from r to r+step)
	dis_ucd_gc_mask = dis_ucd_gc[dis_ucd_gc<(r+step)] #outer bound
	dis_ucd_gc_mask = dis_ucd_gc_mask[dis_ucd_gc_mask>r]  #inner bound
	area = pi*((r+step)**2-r**2)  # area of that ring  (in kpc^2)

	# expected count and measued count of GCs
	gc_expected = (exp_density-bkg_unif)*area
	gc_count = len(dis_ucd_gc_mask)

	gc_count_cor =gc_count -bkg_unif*area

	return gc_count, gc_count_cor, gc_expected, dis_ucd_gc_mask

def C_ave(gc_counts, dis_ucd_gc_list, C):
	'''
	calculate averaged C value among UCDs
	'''
	dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(dis_ucd_gc_list)/len(radius),len(radius)) 
	C = C.reshape(len(C)/len(radius),len(radius)) 
	gc_counts = gc_counts.reshape(len(gc_counts)/len(radius),len(radius)) 
	#C = C[C[:,0]+C[:,1]>0]
	C_mean = np.array([])
	C_std = np.array([])
	for i in range(len(radius)):
		logfile.write('========= C statistics for r ='+str(radius[i])+'============\n')
		C_cloumn = C[:,i]
		C_cloumn, mean, std,log_robust_mean = robust_mean(C[:,i],iter=0,show_step=True)
		if log_robust_mean != ' ': logfile.write(log_robust_mean)
		C_mean = np.append(C_mean,mean)   #robust averge C signal over all UCDs
		C_std = np.append(C_std,std)  #robust averge C signal over all UCDs
	C_max = np.max(C,0)  #max
	C_min = np.min(C,0)    #min
	dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0)
	return C_mean,C_std,dis_ucd_gc_list_mean

def mean_dis(dis_ucd_gc_list, gc_count, dis_ucd_gc_mask):
	'''
	GC distance to host UCD
	'''
	if gc_count == 0:
		dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
	else:
		dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))	

	return dis_ucd_gc_list

for i in range(len(cat_ucd)): 
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	r_h = cat_ucd[i]['RH']
	g_i = cat_ucd[i]['MAGCOR_AP8'][1] - cat_ucd[i]['MAGCOR_AP8'][3] 

	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180.*pi*DIS
	if dis_M87_ucd>190 or dis_M87_ucd<15:
		continue

	exp_density,bkg_unif = bkg_1d(dis_M87_ucd,slope = -1.84539280985, intercept = 5.07330399332)
	exp_density_b,bkg_unif_b = bkg_1d(dis_M87_ucd,slope = -1.65574548112, intercept = 3.98509677054)
	exp_density_r,bkg_unif_r = bkg_1d(dis_M87_ucd,slope = -1.87975286933, intercept = 4.07040541314)

	#exp_density,bkg_unif = bkg_2d(ra_ucd,dec_ucd)                         

	col_value = [ID,ra_ucd,dec_ucd,round(dis_M87_ucd,2),r_h, g_i]  #prepare for the info table
	dis_ucd_gc = np.sqrt((ra_gc-ra_ucd)**2+(dec_gc-dec_ucd)**2)/180.*pi*DIS  # in kpc (list)
	dis_ucd_gc_b = np.sqrt((ra_gc_b-ra_ucd)**2+(dec_gc_b-dec_ucd)**2)/180.*pi*DIS  # in kpc (list)
	dis_ucd_gc_r = np.sqrt((ra_gc_r-ra_ucd)**2+(dec_gc_r-dec_ucd)**2)/180.*pi*DIS  # in kpc (list)

	logfile.write('============'+str(ID)+' '+str(ra_ucd)+' '+str(dec_ucd)+' '+str(round(dis_M87_ucd,4))+'=================\n') 
	for r in radius:
		'''
		  gc_count is a scalar (number of GCs in each bin)
		  gc_counts is a list (number of GCs for every bin)
		'''
		gc_count, gc_count_cor, gc_expected, dis_ucd_gc_mask = gc_stat(dis_ucd_gc, exp_density, bkg_unif)
		gc_count_b, gc_count_cor_b, gc_expected_b, dis_ucd_gc_mask_b = gc_stat(dis_ucd_gc_b, exp_density_b, bkg_unif_b)
		gc_count_r, gc_count_cor_r, gc_expected_r, dis_ucd_gc_mask_r = gc_stat(dis_ucd_gc_r, exp_density_r, bkg_unif_r)

		#Clustering Signal
		C = np.append(C,gc_count_cor/gc_expected) 
		C_b = np.append(C_b,gc_count_cor_b/gc_expected_b) 
		C_r = np.append(C_r,gc_count_cor_r/gc_expected_r) 

		gc_counts = np.append(gc_counts,gc_count)
		gc_counts_b = np.append(gc_counts_b,gc_count_b)
		gc_counts_r = np.append(gc_counts_r,gc_count_r)

		dis_ucd_gc_list = mean_dis(dis_ucd_gc_list, gc_count, dis_ucd_gc_mask)
		dis_ucd_gc_list_b = mean_dis(dis_ucd_gc_list_b, gc_count_b, dis_ucd_gc_mask_b)
		dis_ucd_gc_list_r = mean_dis(dis_ucd_gc_list_r, gc_count_r, dis_ucd_gc_mask_r)

		col_value.append(gc_count)
		col_value.append(gc_count_cor)
		col_value.append(gc_count_cor/gc_expected)
		'''
		col_value.append(gc_count_b)
		col_value.append(gc_count_cor_b)
		col_value.append(gc_count_cor_b/gc_expected_b)
		col_value.append(gc_count_r)
		col_value.append(gc_count_cor_r)
		col_value.append(gc_count_cor_r/gc_expected_r)
		'''
		write_to_log = np.array( [r, gc_count, gc_count_b, gc_count_r, gc_expected, gc_count_cor/gc_expected],dtype=np.float16)
		logfile.write(np.array_str(write_to_log)+'\n')

	ucd_data.add_row(col_value)

logfile.write('=================C clustering signal======================\n')
C_mean, C_std, dis_ucd_gc_list_mean = C_ave(gc_counts, dis_ucd_gc_list, C)
C_mean_b, C_std_b, dis_ucd_gc_list_mean_b = C_ave(gc_counts_b, dis_ucd_gc_list_b, C_b)
C_mean_r, C_std_r, dis_ucd_gc_list_mean_r = C_ave(gc_counts_r, dis_ucd_gc_list_r, C_r)

np.set_printoptions(precision=2)
print 'Number of UCDs:',C.shape[0]
print 'Radius bins:',radius
print '==============================='
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean:',C_mean
print 'Std:',C_std/sqrt(len(cat_ucd)) 
print '==============================='
print 'Blue Mean DIS:',dis_ucd_gc_list_mean_b
print 'Blue Mean:',C_mean_b
print 'Blue Std:',C_std_b/sqrt(len(cat_ucd)) 
print '==============================='
print 'Red Mean DIS:',dis_ucd_gc_list_mean_r
print 'Red Mean:',C_mean_r
print 'Red Std:',C_std_r/sqrt(len(cat_ucd)) 
#print 'Max:',C_max
#print 'Min:', C_min
ucd_data.write('ucd_gc_clustering.fits',overwrite=True)

#plot the final figure
fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/sqrt(len(cat_ucd)),fmt='ko')
ax.errorbar(dis_ucd_gc_list_mean_b,C_mean_b,yerr=C_std_b/sqrt(len(cat_ucd)),fmt='bs')
ax.errorbar(dis_ucd_gc_list_mean_r,C_mean_r,yerr=C_std_r/sqrt(len(cat_ucd)),fmt='r^')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]',fontsize=16)
ax.set_ylabel(r'C=<$\Sigma_{\rm measured}/\Sigma_{\rm expected}$>',fontsize=17)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0.25,4.5)
ax.set_ylim(0,int(C_mean.max())+2)
plt.savefig('clustering.new.blueandred.png')
plt.savefig('./pics/clustering_blueandred'+str_radius[0]+'_'+str(step)+'.eps')

logfile.close()