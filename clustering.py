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
ra = np.array(cat_gc['ra'])      #RA of GCs
dec = np.array(cat_gc['dec'])  #Dec of GCs

step = 0.5  # in kpc
radius = np.arange(0.25,4.25,step)   # in kpc  from 0.25 to 3.75 step 0.5
str_radius = str(radius).strip('[]').split()  #for table name

tab_names = ['ID','ra','dec','DIS']+ str_radius 
tab_dtype = ['a8','f8','f8','f4']+['f4']*len(radius)
ucd_data = Table(names=(tab_names),dtype=(tab_dtype))

C = np.array([])   # C=clustering signal     (measured/bkg)
dis_ucd_gc_list = np.array([])  # list of each ucd's distance to M87
for i in range(len(cat_ucd)): 
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']

	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180*pi*DIS  # in kpc (scalar)
	dis_ucd_gc = np.sqrt((ra-ra_ucd)**2+(dec-dec_ucd)**2)/180*pi*DIS  # in kpc (list)
	#if dis_M87_ucd>250 or dis_M87_ucd<20:
	#	continue
	exp_density = exp(5.1902)*dis_M87_ucd**(-1.88)    # expected density (assume uniform in the vicinity of a UCD)
	bkg_unif = exp(5.1902)*250**(-1.88)                       # uniform bakground in the field (possibly non-GC objects)
	col_value = [ID,ra_ucd,dec_ucd,round(dis_M87_ucd,2)]

	dis_ucd_gc_mask = np.array([])
	print '============',str(ID),ra_ucd,dec_ucd,round(dis_M87_ucd,4),'================='
	for r in radius:
                          #mask out the area for counting  (ring shaped mask from r to r+step)
		dis_ucd_gc_mask = dis_ucd_gc[dis_ucd_gc<(r+step)] #outer bound
		dis_ucd_gc_mask = dis_ucd_gc_mask[dis_ucd_gc_mask>r]  #inner bound
		area = pi*((r+step)**2-r**2)  # area of that ring  (in kpc^2)

		# expected count and measued count of GCs
		gc_expected = (exp_density-bkg_unif)*area
		gc_count = len(dis_ucd_gc_mask)-bkg_unif*area
		C = np.append(C,gc_count/gc_expected) 

		# gc distance to host UCD
		if len(dis_ucd_gc_mask) == 0:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
		else:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))

		col_value.append(gc_count/gc_expected)
		print ID,r,round(len(dis_ucd_gc_mask),4),round(gc_expected,4),round(gc_count/gc_expected,4)

	ucd_data.add_row(col_value)

'''
   This part is statistics for the clustering signal
    Mean Value, Standard Error, etc...
    implemented  robust statistics
'''
C = C.reshape(len(C)/len(radius),len(radius)) 
#C = C[C[:,0]+C[:,1]>0]
C_mean = np.array([])
C_std = np.array([])
for i in range(len(radius)):# histagram for C signal distribution 
	print '========= C statistics for r =',radius[i],'============'
	fig, ax = plt.subplots()
	ax.hist(C[:,i],bins=50)
	ax.set_title(str(radius[i])+'~'+str(radius[i]+step))
	plt.savefig('C.hist.'+str(radius[i])+'.png')
	mean_r,std_r = robust_mean(C[:,i],iter=10,show_step=True)
	C_mean = np.append(C_mean,mean_r)   #robust averge C signal over all UCDs
	C_std = np.append(C_std,std_r)   #robust averge C signal over all UCDs
C_max = np.max(C,0)  #max
C_min = np.min(C,0)  #min
dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(dis_ucd_gc_list)/len(radius),len(radius)) 
dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0)

print '=================C clustering signal======================'
np.set_printoptions(precision=2)
print 'radius:',radius
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean:',C_mean
print 'Std:',C_std
print 'Max:',C_max
print 'Min:', C_min
ucd_data.write('ucd_gc_clustering.fits',overwrite=True)

#plot the final figure
fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/sqrt(len(cat_ucd)),fmt='o')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]')
ax.set_ylabel(r'C=<$\Sigma_{measured}/\Sigma_{expected}$>',fontsize=16)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0,4.25)
ax.set_ylim(0,10)
plt.savefig('clustering.new.png')