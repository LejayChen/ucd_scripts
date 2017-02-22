import aplpy 
from astropy.table import Table
from astropy.io import fits
from math import pi,sqrt

cat_ucd = Table.read('NGVS.pilot.92ucds.fits')  #catalog of UCDs

cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]   #masked (p_gc>0.95) catalog of GCs

cat_v_ucd = Table.read('match_ucd_v.fits')
cat_v_gc = Table.read('match_v_radial_gc.fits')
cat_v_gc = cat_v_gc[cat_v_gc['p_gc']>0.95]

aperture = 50./3600   # in degree
DIS = 17.21*1000 #kpc
M87 = (187.70583,  12.39111)

def pixel_max(filename):
	im = fits.open('ucd.'+filename+'.fits')
	data = im[0].data
	pixel_max = data.max()
	x,y = data.shape
	pixel_value_central = data[x/2-2:x/2+2,y/2-3:y/2+2].max()
	return pixel_value_central*1.1

def mask(cat):
	mask1 =  abs(cat['ra'] -ra_ucd)<aperture
	cat_gc_around_ucd = cat[mask1]
	mask2 =  abs(cat_gc_around_ucd['dec'] -dec_ucd)<aperture 
	cat_gc_around_ucd = cat_gc_around_ucd[mask2]
	return cat_gc_around_ucd

for i in range(len(cat_v_ucd)):
	INDEX = cat_v_ucd[i]['INDEX']
	ra_ucd = cat_v_ucd[i]['RA']
	dec_ucd = cat_v_ucd[i]['DEC']
	dis = round(sqrt((ra_ucd-M87[0])**2+(dec_ucd-M87[1])**2)*pi/180*DIS,1)
	v = round(cat_v_ucd[i]['v_msc'],2)
	g_i_ucd = cat_v_ucd[i]['MAG_BEST'][1] - cat_ucd[i]['MAG_BEST'][3]  

             #mask the GC catalog around UCD
	cat_gc_around_ucd = mask(cat_gc)
	ra_gc = cat_gc_around_ucd['ra']
	dec_gc = cat_gc_around_ucd['dec']
	g_i_gc = cat_gc_around_ucd['gmag'] - cat_gc_around_ucd['imag']

             #mask the GC with radial velocity catalog around UCD
	cat_gc_ucd_v = mask(cat_v_gc)
	ra_gc_v = cat_gc_ucd_v['ra']
	dec_gc_v = cat_gc_ucd_v['dec']
	gc_v = cat_gc_ucd_v['v_msc']

	filename = str(INDEX)
	fig = aplpy.FITSFigure('ucd_newcut.'+filename+'.fits')
	fig.tick_labels.set_xformat('dd.ddd')
	fig.tick_labels.set_yformat('dd.ddd')

	for i in range(len(ra_gc_v)):
		fig.add_label(ra_gc_v[i],dec_gc_v[i] + 0.0004,text=str(gc_v[i]),color='red')
	for i in range(len(ra_gc)):
		fig.add_label(ra_gc[i],dec_gc[i] - 0.0004,text=str(g_i_gc[i]),color='green')
	fig.add_label(ra_ucd,dec_ucd+0.0004,text=str(v))
	fig.add_label(ra_ucd,dec_ucd-0.0004,text=str(g_i_ucd),color='green')

	fig.show_grayscale(stretch='log',vmin=0,vmid=-0.3,vmax=pixel_max(filename),pmin=30,invert='y')
	fig.show_circles(ra_ucd,dec_ucd,0.25/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,0.75/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,1.25/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,1.75/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,2.25/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,2.75/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,3.25/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,3.75/DIS/pi*180,color='green')
	fig.show_circles(ra_ucd,dec_ucd,4.25/DIS/pi*180,color='green')
	if len(cat_gc_around_ucd)>0: fig.show_circles(ra_gc,dec_gc,0.8/3600,color='red')
	fig.set_title('UCD '+filename+'  '+str(dis))
	fig.save('ucd_gctag.'+filename+'.png')