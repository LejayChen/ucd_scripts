import aplpy 
from astropy.table import Table
from astropy.io import fits
from math import pi
cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs
aperture =12./3600   # in degree
DM = 22.2*1000 #kpc

def pixel_max(filename):
	im = fits.open('ucd.'+filename+'.fits')
	data = im[0].data
	pixel_max = data.max()
	x,y = data.shape
	pixel_value_central = data[x/2-2:x/2+2,y/2-3:y/2+2].max()
	return pixel_value_central*1.1

for i in range(len(cat_ucd)):
	INDEX = cat_ucd[i]['INDEX']
	ra_ucd = cat_ucd[i]['RA']
	dec_ucd = cat_ucd[i]['DEC']

	mask1 =  abs(cat_gc['ra'] -ra_ucd)<aperture
	cat_gc_around_ucd = cat_gc[mask1]
	mask2 =  abs(cat_gc_around_ucd['dec'] -dec_ucd)<aperture 
	cat_gc_around_ucd = cat_gc_around_ucd[mask2]
	ra_gc = cat_gc_around_ucd['ra']
	dec_gc = cat_gc_around_ucd['dec']


	filename = str(INDEX)
	fig = aplpy.FITSFigure('ucd_newcut.'+filename+'.fits')
	fig.show_grayscale(stretch='log',vmin=0,vmid=-0.3,vmax=pixel_max(filename),pmin=30,invert='y')
	fig.show_circles(ra_ucd,dec_ucd,0.25/DM/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,0.75/DM/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,1.25/DM/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,1.75/DM/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,2.25/DM/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,2.75/DM/pi*180,color='green')
	if len(cat_gc_around_ucd)>0: fig.show_circles(ra_gc,dec_gc,0.8/3600,color='red')
	fig.set_title('UCD '+filename)
	fig.save('ucd_gctag.'+filename+'.png')