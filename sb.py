import os
from astropy.table import *
from astropy.io import fits
import numpy as np 
from math import *
import matplotlib.pyplot as plt 

cat_dEN = Table.read('../pp.gal.nuc2.s.master.new.fits')
cat_ucd = Table.read('../NGVS.pilot.92ucds.fits')  

def run_sex(FileName,cat_type):
	assoc_name = 'position'
	cat_name = './sex_fits/'+cat_type+'.'+FileName+'.fits '
	check_img = './sex_fits/'+ cat_type+'.'+FileName+'.check.fits'
	if cat_type == 'ucd':
		os.system('sextractor  ../ucd_newcut.'+FileName+'.fits'
			+' -c ngvs.sex -CATALOG_NAME ' + cat_name
			+' -ASSOC_NAME ' + assoc_name
			+' -CHECKIMAGE_NAME '+check_img)
	elif cat_type == 'dEN':
		os.system('sextractor  ../dEN_pics/dE.N.'+FileName+'.fits'
			+' -c ngvs.sex -CATALOG_NAME ' + cat_name
			+' -ASSOC_NAME ' + assoc_name
			+' -CHECKIMAGE_NAME '+check_img)
	else:
		raise KeyError('cat_type argument only takes ucd or dEN')

def growth_curve(ID, cat_type, radii,mags):
	radii = radii[mags !=99.0]
	mags = mags[mags != 99.0]

	r_edge = radii[-1]
	for i in range(len(mags[:-1])):
		if abs(mags[i+1] - mags[i])/(radii[i+1]-radii[i])<0.012: 
			r_edge = radii[i+1]
			plt.plot(radii[i+1],mags[i+1],'^g')
			break

	for i in range(len(mags[:-1])):
		if abs((mags[i+1] - mags[i])/(radii[i+1] - radii[i]))<0.1:
			r1 = radii[i]
			break

	r2 = 1.01*r_edge
	index1  = len(radii[radii<r1]) - 1
	index2  = len(radii[radii<r2]) - 1	

	plt.plot(radii,mags,'.-')
	plt.plot(radii[index1],mags[index1],'.r')
	plt.plot(radii[index2],mags[index2],'.r')
	plt.title(cat_type+str(ID))
	plt.gca().invert_yaxis()
	plt.savefig('growth_curves/'+cat_type+'/'+ID+'.growth_curve.png')
	plt.clf()

	return index1, index2

def cal_sb(ID, ra, dec, cat_type):
	catalog = Table.read('sex_fits/'+cat_type+'.'+str(ID)+'.fits')
	obj_id =  0

	radii = np.array([2,4,6, 8, 10,12,14,16,18, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96,110,130,150])*0.187/2  # in arcsec
	mags = catalog[obj_id]['MAG_APER']
	fluxes = 10**((mags - 30)/(-2.5))

	index1, index2 = growth_curve(str(ID), cat_type, radii,mags)
	sb = - 2.5*np.log10((fluxes[index2] - fluxes[index1])/(pi*(radii[index2]**2 - radii[index1]**2)))+30

	return sb

def write_table(sb_column, cat_type):
	sb_column = Column(data = np.array(sb_column), name='sb')
	if cat_type =='ucd':
		cat_ucd.add_column(sb_column)
		cat_ucd.write('../new.NGVS.pilot.92ucds.fits',overwrite=True)
	if cat_type == 'dEN':
		cat_dEN.add_column(sb_column)
		cat_dEN.write('../new.pp.gal.nuc2.s.master.new.fits',overwrite=True)

def plot_hist(sbs,cat_type):
	sbs = sbs[~np.isnan(sbs)]
	plt.hist(sbs)
	plt.xlabel(r'surface brightness ${\rm mag/arcsec^2}$')
	plt.title(cat_type + ' surface brightness histogram')
	plt.savefig('../pics/'+cat_type+'_sb_hist.png')
	plt.clf()
	print cat_type+' fig saved'

for cat in [(cat_ucd,'ucd'), (cat_dEN,'dEN')]:
	sb_column = []
	for obj in cat[0]:
		cat_type = cat[1]
		ID = obj['INDEX']
		#run_sex(str(ID), str(cat_type))
		ra,dec = obj['RA'],obj['DEC']
		sb = cal_sb(ID, ra, dec, cat_type)
		sb_column.append(sb)
		print cat_type, ID,sb

	write_table(sb_column, cat_type)
	plot_hist(np.array(sb_column),cat_type)