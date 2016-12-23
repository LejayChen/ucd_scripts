from astropy.table import *
from math import *
from robust_mean import *
import numpy as np
import matplotlib.pyplot as plt
DM = 22.2*1000 #kpc
M87 = (187.70583, 12.39111)

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra = np.array(cat_gc['ra'])  #RA of GCs
dec = np.array(cat_gc['dec'])  #Dec of GCs

C = np.array([])   # C=clustering signal     measured/bkg
dis_ucd_gc_list = np.array([])  # list of ucd's distance to M87

step = 0.5
radius = np.arange(0.25,4.25,step)   # in kpc  from 0.25 to 3.75 step 0.5

str_radius = str(radius).strip('[]').split()  #for table name
tab_names = ['ID','ra','dec','DIS']+ str_radius 
tab_dtype = ['a8','f8','f8','f4']+['i2']*len(radius)
ucd_data = Table(names=(tab_names),dtype=(tab_dtype))
for i in range(len(cat_ucd)): 
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180*pi*DM  # in kpc
	dis_ucd_gc = np.sqrt((ra-ra_ucd)**2+(dec-dec_ucd)**2)/180*pi*DM  # in kpc
	bkg_density = exp(5.1902)*dis_M87_ucd**(-1.88)    # assume uniform in the vicinity of a UCD
	col_value = [ID,ra_ucd,dec_ucd,round(dis_M87_ucd,2)]
	print '================',str(ID),ra_ucd,dec_ucd,dis_M87_ucd,'========================'
	for r in radius:
		bkg = bkg_density*pi*((r+step)**2-r**2)
		mask = dis_ucd_gc<(r+step)  
		dis_ucd_gc_mask = dis_ucd_gc[mask]
		mask2 = dis_ucd_gc_mask>r
		dis_ucd_gc_mask = dis_ucd_gc_mask[mask2]

		gc_count = len(dis_ucd_gc_mask)
		C = np.append(C,gc_count/bkg)  
		if gc_count == 0:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.NaN)
		else:
			dis_ucd_gc_list = np.append(dis_ucd_gc_list, np.mean(dis_ucd_gc_mask))

		col_value.append(gc_count)
		print ID,r,gc_count,bkg,gc_count/bkg

	ucd_data.add_row(col_value)

C = C.reshape(len(cat_ucd),len(radius)) 
#C = C[C[:,0]+C[:,1]>0]
C_mean = np.array([])
for i in range(len(radius)): # histagram for C signal distribution 
	fig, ax = plt.subplots()
	ax.hist(C[:,i],bins=50)
	ax.set_title(str(radius[i])+'~'+str(radius[i]+step))
	plt.savefig('C.hist.'+str(radius[i])+'.png')
	C_mean = np.append(C_mean,robust_mean(C[:,i]))   #robust averge C signal over all UCDs

C_max = np.max(C,0)  #max
C_min = np.min(C,0)  #min
C_std = np.std(C,0)  # standard error
dis_ucd_gc_list = dis_ucd_gc_list.reshape(len(cat_ucd),len(radius)) 
dis_ucd_gc_list_mean = np.nanmean(dis_ucd_gc_list,0)

print '=================C======================'
print 'radius:',radius
print 'Mean DIS:',dis_ucd_gc_list_mean
print 'Mean:',C_mean
print 'Std:',C_std
print 'Max:',C_max
print 'Min:', C_min
ucd_data.write('ucd_gc_clustering.fits',overwrite=True)

fig, ax = plt.subplots()
ax.errorbar(dis_ucd_gc_list_mean,C_mean,yerr=C_std/10,fmt='o')
ax.axhline(y=1,xmin=0,xmax=4,color='k')
ax.set_xlabel('Aperture around UCD [kpc]')
ax.set_ylabel(r'C=<$\Sigma_{measured}/\Sigma_{expected}$>',fontsize=16)
ax.tick_params(which='major',length=12) 
ax.tick_params(which='minor',length=5) 
ax.set_xlim(0,4.25)
ax.set_ylim(0,10)
plt.savefig('clustering.new.png')