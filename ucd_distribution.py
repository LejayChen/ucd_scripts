import matplotlib.pyplot as plt 
from astropy.table import Table 
from math import pi

cat = Table.read('NGVS.pilot.92ucds.fits')
cat2 = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat2 = cat2[cat2['p_gc']>0.95]
cat3 = Table.read('../dEdata.fits')

M87 = (187.7058,12.3911)
DM = 22.2*1000 #kpc
#major axis: 499.10 arcsec
#minor axis: 394.29 arcsec

fig, ax = plt.subplots()
ax.set_aspect('equal')
circle1 = plt.Circle((M87[0],M87[1]),radius=20./DM/pi*180,fill=False)
circle1.set_zorder(100)
circle2 = plt.Circle((M87[0],M87[1]),radius=250./DM/pi*180,fill=False)
ax.add_artist(circle1)
ax.add_artist(circle2)

ra2 = cat2['ra']
dec2 = cat2['dec']
plt.plot(ra2,dec2,'.k',color='0.7',markersize=2,label='GCs')
ra3 = cat3['ra']
dec3 = cat3['dec']
plt.plot(ra3,dec3,'.b',label='dE,Ns')
ra = cat['RA']
dec = cat['DEC']
plt.plot(ra,dec,'.r',label='UCDs')
plt.plot(M87[0],M87[1],'+b',markersize=12,markeredgewidth=2,label='M87')


plt.xlim([M87[0]-1.12,M87[0]+1.12])
plt.ylim([M87[1]-1,M87[1]+1])
plt.xlabel('RA [Degree]')
plt.ylabel('Dec [Degree]')
plt.grid(True)
plt.legend(numpoints=1,frameon=True,loc='upper right')
plt.show()