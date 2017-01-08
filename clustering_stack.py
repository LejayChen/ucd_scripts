from astropy.table import *
from math import *
from robust_mean import *
import numpy as np
import matplotlib.pyplot as plt
DIS = 22.2*1000 #kpc
M87 = (187.70583, 12.39111)

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra_gc = np.array(cat_gc['ra'])      #RA of GCs
dec_gc = np.array(cat_gc['dec'])  #Dec of GCs

step = 0.5  # in kpc
radius = np.arange(0.25,4.25,step)   # in kpc  from 0.25 to 3.75 step 0.5
str_radius = str(radius).strip('[]').split()  #for table name

gc_count_total = 0
gc_exps = np.array([])   # C=clustering signal     (measured/bkg)
gc_counts = np.array([])  # count of GCs 
dis_ucd_gc_list = np.array([])  # list of each gc's distance to UCD

for i in range(len(cat_ucd)): 
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	r_h = cat_ucd[i]['RH']

	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180.*pi*DIS  # in kpc (scalar)
	dis_ucd_gc = np.sqrt((ra_gc-ra_ucd)**2+(dec_gc-dec_ucd)**2)/180.*pi*DIS  # in kpc (list)
	if dis_M87_ucd>250. or dis_M87_ucd<20.:
		continue
	exp_density = exp(5.09201844251)*dis_M87_ucd**(-1.85787561291)   # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(5.09201844251)*270**(-1.85787561291)    # uniform bakground in the field (possibly non-GC objects)

	dis_ucd_gc_mask = np.array([])
	print '============',str(ID),ra_ucd,dec_ucd,round(dis_M87_ucd,4),'================='
	for r in radius:
                          #mask out the area for counting  (ring shaped mask from r to r+step)
		dis_ucd_gc_mask = dis_ucd_gc[dis_ucd_gc<(r+step)] #outer bound
		dis_ucd_gc_mask = dis_ucd_gc_mask[dis_ucd_gc_mask>r]  #inner bound
		area = pi*((r+step)**2-r**2)  # area of that ring  (in kpc^2)

		# expected count and measued count of GCs
		gc_expected = (exp_density-bkg_unif)*area
		gc_count = len(dis_ucd_gc_mask)
		gc_count_total += gc_count
		gc_count_cor =gc_count -bkg_unif*area
		#Clustering Signal

		gc_counts = np.append(gc_counts,gc_count_cor)
		gc_exps = np.append(gc_exps,gc_expected)

		# gc distance to host UCD
		if gc_count == 0:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
		else:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))

		print ID,r,round(gc_count,4),round(gc_expected,4),round(gc_count_cor/gc_expected,4)

'''
   This part is statistics for the clustering signal
    Mean Value, Standard Error, etc...
    implemented  robust statistics
'''

dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(dis_ucd_gc_list)/len(radius),len(radius)) 
gc_exps = gc_exps.reshape(len(gc_exps)/len(radius),len(radius)) 
gc_counts = gc_counts.reshape(len(gc_counts)/len(radius),len(radius)) 

gc_exps_mean = np.array([])
gc_exps_std = np.array([])
gc_counts_mean = np.array([])
gc_counts_std = np.array([])
for i in range(len(radius)):# histagram for C signal distribution 
	print '========= C statistics for r =',radius[i],'============'
	gc_exps_r = gc_exps[:,i]
	gc_exps_r,exps_mean_r,exps_std_r = robust_mean(gc_exps[:,i],iter=0,show_step=True)
	gc_counts_r = gc_counts[:,i]
	gc_counts_r,counts_mean_r,counts_std_r = robust_mean(gc_counts[:,i],iter=0,show_step=True)

	gc_exps_mean = np.append(gc_exps_mean,exps_mean_r)   #robust averge gc_exps signal over all UCDs
	gc_exps_std = np.append(gc_exps_std,exps_std_r)  #robust averge C signal over all UCDs
	gc_counts_mean = np.append(gc_counts_mean,counts_mean_r)   #robust averge gc_counts signal over all UCDs
	gc_counts_std = np.append(gc_counts_std,counts_std_r)  #robust averge C signal over all UCDs

C_mean = gc_counts_mean/gc_exps_mean
C_std = gc_exps_std/gc_exps_std
dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0)

print '=================C clustering signal======================'
np.set_printoptions(precision=2)
print 'Number of UCDs:',gc_counts.shape[0],', Number of GCs:',gc_count_total
print 'Radius bins:',radius
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean:',C_mean
print 'Std:',C_std/sqrt(len(cat_ucd)) 

#plot the final figure
fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/sqrt(len(cat_ucd)),fmt='o')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]')
ax.set_ylabel(r'C=<$\Sigma_{measured}/\Sigma_{expected}$>',fontsize=16)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0,4.25)
ax.set_ylim(0,int(C_mean.max())+2)
plt.savefig('clustering.new.stack.png')