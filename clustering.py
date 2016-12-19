from astropy.table import *
import numpy as np
from math import *
import matplotlib.pyplot as plt
DM = 22.2*1000 #kpc
M87 = (187.70583, 12.39111)

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
ra = np.array(cat_gc['ra'])  #RA of GCs
dec = np.array(cat_gc['dec'])  #Dec of GCs

C = np.array([])   # C=clustering signal     measured/bkg
radius = np.arange(0.25,4.25,0.5)   # in kpc  from 0.25 to 3.75 step 0.5

str_radius = str(radius).strip('[]').split()  #for table name
tab_names = ['ID','ra','dec']+ str_radius 
tab_dtype = ['a8','f8','f8']+['i2']*len(radius)
ucd_data = Table(names=(tab_names),dtype=(tab_dtype))
for i in range(len(cat_ucd)):
	ID = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']
	dis_M87_ucd = sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)/180*pi*DM
	bkg_density = exp(5.192)*dis_M87_ucd**(-1.89)    # assume uniform in the vicinity of a UCD
	col_value = [ID,ra_ucd,dec_ucd]
	for r in radius:
		bkg = bkg_density*pi*r**2
		mask = np.sqrt((ra-ra_ucd)**2+(dec-dec_ucd)**2)<(r+0.5)/DM/pi*180  # radius in angular seperation (degree)
		mask2 = np.sqrt((ra-ra_ucd)**2+(dec-dec_ucd)**2)<r/DM/pi*180
		C = np.append(C,(len(cat_gc[mask])-len(cat_gc[mask2]))/bkg)  
		col_value.append(len(cat_gc[mask])-len(cat_gc[mask2]))
		#print ID,r,bkg,(len(cat_gc[mask])-len(cat_gc[mask2]))/bkg
	ucd_data.add_row(col_value)
		
C = C.reshape(len(cat_ucd),len(radius))  #  C in reserve order
C_mean = np.mean(C,0)   #averge C signal over all UCDs
C_max = np.max(C,0)
C_min = np.min(C,0)
C_std = np.std(C,0)  # standard error

print '=================C======================'
print 'Mean:',C_mean
print 'Std:',C_std
print 'Max:',C_max
print 'Min:', C_min
ucd_data.write('ucd_gc_clustering.fits',overwrite=True)

plt.errorbar(radius+0.25,C_mean,yerr=C_std/10,fmt='o')
plt.axhline(y=1,xmin=0,xmax=4,color='k')
plt.xlabel('Aperture around UCD [kpc]')
plt.ylabel(r'C=<$\Sigma_{measured}/\Sigma_{expected}$>',fontsize=16)
plt.tick_params(which='major',length=12) 
plt.tick_params(which='minor',length=5) 
plt.xlim(0,4.25)
plt.ylim(0,10)
plt.show()