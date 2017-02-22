from astropy.table import *
from math import *
from robust_mean import *
from func import func
import numpy as np
import matplotlib.pyplot as plt
DIS = 17.21*1000 #kpc
M87 = (187.70583, 12.39111)

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra_gc = np.array(cat_gc['ra'])      #RA of GCs
dec_gc = np.array(cat_gc['dec'])  #Dec of GCs

step = 0.5  # in kpc
radius = np.arange(0.25,4.25,step)   # in kpc  from 0.25 to 3.75 step 0.5
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

for i in range(len(cat_ucd)): 
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	r_h = cat_ucd[i]['RH']
	g_i = cat_ucd[i]['MAG_BEST'][1] - cat_ucd[i]['MAG_BEST'][3] 

	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180.*pi*DIS  # in kpc (scalar)
	dis_ucd_gc = np.sqrt((ra_gc-ra_ucd)**2+(dec_gc-dec_ucd)**2)/180.*pi*DIS  # in kpc (list)
	if dis_M87_ucd>190 or dis_M87_ucd<15:
		continue
	
	exp_density = exp(5.10118386605)*dis_M87_ucd**(-1.85088530618 )   # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(5.10118386605)*200**(-1.85088530618)    # uniform bakground in the field (possibly non-GC objects)	
	#exp_density = func(((ra_ucd-M87[0])/180*pi*DIS,(dec_ucd-M87[1])/180*pi*DIS), -1.732794,  123.30887701)  
	#bkg_unif = 0                                                                        
	col_value = [ID,ra_ucd,dec_ucd,round(dis_M87_ucd,2),r_h, g_i]

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
		C = np.append(C,gc_count_cor/gc_expected) 
		gc_counts = np.append(gc_counts,gc_count)
		col_value.append(gc_count)
		col_value.append(gc_count_cor)
		col_value.append(gc_count_cor/gc_expected)

		# gc distance to host UCD
		if gc_count == 0:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
		else:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))

		print ID,r,round(gc_count,4),round(gc_expected,4),round(gc_count_cor/gc_expected,4)

	ucd_data.add_row(col_value)

'''
   This part is statistics for the clustering signal
    Mean Value, Standard Error, etc...
    implemented  robust statistics
'''

dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(dis_ucd_gc_list)/len(radius),len(radius)) 
C = C.reshape(len(C)/len(radius),len(radius)) 
gc_counts = gc_counts.reshape(len(gc_counts)/len(radius),len(radius)) 
#C = C[C[:,0]+C[:,1]>0]
C_mean = np.array([])
C_std = np.array([])
for i in range(len(radius)):# histogram for C signal distribution 
	print '========= C statistics for r =',radius[i],'============'
	C_r = C[:,i]
	C_r,mean_r,std_r = robust_mean(C[:,i],iter=0,show_step=True)

	fig, ax = plt.subplots()
	ax.hist(C_r,bins=10)
	ax.set_title(str(radius[i])+'~'+str(radius[i]+step))
	plt.savefig('C.hist.'+str(radius[i])+'.png')

	fig2, ax2 = plt.subplots()
	ax2.hist(gc_counts[:,i],bins=range(0, 8, 1))
	ax2.set_title(str(radius[i])+'~'+str(radius[i]+step))
	plt.savefig('gc.count.hist.'+str(radius[i])+'.png')

	C_mean = np.append(C_mean,mean_r)   #robust averge C signal over all UCDs
	C_std = np.append(C_std,std_r)  #robust averge C signal over all UCDs
C_max = np.max(C,0)  #max
C_min = np.min(C,0)    #min
dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0)

print '=================C clustering signal======================'
np.set_printoptions(precision=2)
print 'Number of UCDs:',C.shape[0],', Number of GCs:',gc_count_total
print 'Radius bins:',radius
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean:',C_mean
print 'Std:',C_std/sqrt(len(cat_ucd)) 
print 'Max:',C_max
print 'Min:', C_min
ucd_data.write('ucd_gc_clustering.fits',overwrite=True)

#plot the final figure
fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/sqrt(len(cat_ucd)),fmt='o')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]',fontsize=16)
ax.set_ylabel(r'C=<$\Sigma_{\rm measured}/\Sigma_{\rm expected}$>',fontsize=17)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0,4.25)
ax.set_ylim(0,int(C_mean.max())+2)
plt.savefig('clustering.new.png')
plt.savefig('./pics/clustering_'+str_radius[0]+'_'+str(step)+'.eps')