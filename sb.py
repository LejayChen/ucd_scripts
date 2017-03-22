import os
from astropy.table import *
import numpy as np 

cat_dEN = Table.read('../pp.gal.nuc2.s.master.new.fits')
cat_ucd = Table.read('../NGVS.pilot.92ucds.fits')  #catalog of UCDs

def run_sex(FileName,cat_type):
	if cat_type == 'ucd':
		os.system('sextractor  ../ucd_newcut.'+FileName+'.fits'+' -c ngvs.sex -CATALOG_NAME '+'./sex_fits/'+cat_type+'.'+FileName+'.fits')
	elif cat_type == 'dEN':
		os.system('sextractor  ../dEN_pics/dE.N.'+FileName+'.fits'+' -c ngvs.sex -CATALOG_NAME '+'./sex_fits/'+cat_type+'.'+FileName+'.fits')
	else:
		raise KeyError('cat_type argument only takes ucd or dEN')

def find_central_obj(ra,dec):
	
	return obj_id

def read_tables(ID, ra, dec, cat_type):
	catalog = Table.read('./sex_fits/'+cat_type+'.'+ID+'.fits')
	obj_id =  find_central_obj(ra,dec)

	mag_r1 = catalog[obj_id]['mag_APER'][8]
	mag_r2 = catalog[obj_id]['mag_APER'][9]
	sb = cal_sb(mag_r1, mag_r2)
	sb_column = np.append(sb_column, sb)

	return sb_column

def write_table(sb_column, cat_type):
	return 0

for cat in [(cat_ucd,'ucd'), (cat_dEN,'dEN')]:
	for obj in cat[0]:
		cat_type = cat[1]
		ID = obj['INDEX']
		ra,dec = obj['RA'],obj['DEC']
		run_sex(str(ID), str(cat_type))
		sb_cloumn = read_tables(ID, ra, dec, cat_type)
	write_table(sb_column, cat_type)