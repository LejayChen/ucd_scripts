import os
from astropy.table import *
from astropy.io import fits
import numpy as np 
from math import *

cat_dEN = Table.read('../pp.gal.nuc2.s.master.new.fits')
cat_ucd = Table.read('../NGVS.pilot.92ucds.fits')  #catalog of UCDs

def run_sex(FileName,cat_type):
	if cat_type == 'ucd':
		os.system('sextractor  ../ucd_newcut.'+FileName+'.fits'+' -c ngvs.sex -CATALOG_NAME '+'./sex_fits/'+cat_type+'.'+FileName+'.fits '+'-CHECKIMAGE_NAME '+'./sex_fits/'+ cat_type+'.'+FileName+'.check.fits')
	elif cat_type == 'dEN':
		os.system('sextractor  ../dEN_pics/dE.N.'+FileName+'.fits'+' -c ngvs.sex -CATALOG_NAME '+'./sex_fits/'+cat_type+'.'+FileName+'.fits '+'-CHECKIMAGE_NAME '+'./sex_fits/'+ cat_type+'.'+FileName+'.check.fits')
	else:
		raise KeyError('cat_type argument only takes ucd or dEN')

def img_scale(img_name):
	header = fits.open(img_name)[0].header
	x_width = header['NAXIS1']
	y_width = header['NAXIS2']
	return x_width, y_width

def find_central_obj(catalog, ra, dec,cat_type):
	x,y = catalog['X_IMAGE'], catalog['Y_IMAGE']
	if cat_type == 'ucd':
		img_width_x, img_width_y = img_ucd_width_x, img_ucd_width_y
	elif cat_type == 'dEN':
		img_width_x, img_width_y = img_dEN_width_x, img_dEN_width_y

	limit = 20
	x_slice = np.copy(x)
	while len(x_slice)>1:
		x_slice = x[abs(x - img_width_x/2.)<limit]
		y_slice = y[abs(x - img_width_x/2.)<limit]
		x_slice = x_slice[abs(y_slice - img_width_y/2.)<limit]		
		y_slice = y_slice[abs(y_slice - img_width_y/2.)<limit]
		limit = limit/2. + 0.5

	return np.where(x==x_slice[0])[0][0]

def read_tables(ID, ra, dec, cat_type):
	catalog = Table.read('./sex_fits/'+cat_type+'.'+str(ID)+'.fits')
	obj_id =  find_central_obj(catalog, ra, dec,cat_type)

	r1 = 32*0.187/2
	r2 = 64*0.187/2

	mag_r1 = catalog[obj_id]['MAG_APER'][8]
	mag_r2 = catalog[obj_id]['MAG_APER'][9]

	flux_r1 = 10**((mag_r1 - 30)/(-2.5))
	flux_r2 = 10**((mag_r2 - 30)/(-2.5))

	sb = -2.5*np.log10((flux_r2-flux_r1)/(pi*(r2**2-r1**2)))+30         
	return sb, flux_r1, flux_r2 

def write_table(sb_column, cat_type):
	sb_column = Column(data = np.array(sb_column), name='sb')
	if cat_type =='ucd':
		cat_ucd.add_column(sb_column)
		cat_ucd.write('../new.NGVS.pilot.92ucds.fits',overwrite=True)
	if cat_type== 'dEN':
		cat_dEN.add_column(sb_column)
		cat_dEN.write('../new.pp.gal.nuc2.s.master.new.fits',overwrite=True)

img_ucd_width_x, img_ucd_width_y = img_scale('../ucd_newcut.31819.fits')
img_dEN_width_x, img_dEN_width_y = img_scale('../dEN_pics/dE.N.20157.fits')
for cat in [(cat_ucd,'ucd'), (cat_dEN,'dEN')]:
	sb_column = []
	for obj in cat[0]:
		cat_type = cat[1]
		ID = obj['INDEX']
		#run_sex(str(ID), str(cat_type))
		ra,dec = obj['RA'],obj['DEC']
		sb,f1,f2 = read_tables(ID, ra, dec, cat_type)
		print cat_type, ID, sb,f1,f2
		sb_column.append(sb)
	write_table(sb_column, cat_type)