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

def plot_growth_curve(index1,index2,radii, mags, mag_auto, fluxes, sbs):

	fig, ax1 = plt.subplots()
	plt.gca().invert_yaxis()

	ax1.plot(radii,mags,'.-k')
	#ax1.plot(radii[:-1],sbs,'.-')
	ax1.plot(radii[index1],mags[index1],'.r')
	ax1.plot(radii[index2],mags[index2],'.r')
	ax1.axhline(y=mag_auto, xmin=0, xmax = 10, linewidth=1, linestyle='dashed', color='k')
	ax1.set_title(cat_type+str(ID))
	ax1.set_xlabel('aperture radius [arcsec]',fontsize=14)
	ax1.set_ylabel(r'${\rm m}_g$',fontsize=14)
	ax1.tick_params(direction='in', width=1)
	
	ax2 = ax1.twinx()
	ax2.plot(radii,fluxes,'.-')
	ax2.set_ylabel('flux [unit]')
	ax2.tick_params(direction='in', width=1)

	fig.tight_layout()
	fig.savefig('growth_curves/'+cat_type+'/'+str(ID)+'.growth_curve.png')
	fig.clf()

def growth_curve(ID, cat_type, radii, mags, fluxes, mag_auto):
	radii = radii[mags !=99.0]
	mags = mags[mags != 99.0]
	fluxes = fluxes[mags != 99.0]

	radii1 = radii[:-1]
	radii2 = radii[1:]
	fluxes1 = fluxes[:-1]
	fluxes2 = fluxes[1:]
	sbs = - 2.5*np.log10((fluxes2 - fluxes1)/(pi*(radii2**2 - radii1**2)))+30
	for i in range(len(sbs)): 
		if np.isnan(sbs[i]): sbs[i] = sbs[i-1] 

	r_edge = radii[-1]
	for i in range(len(mags[:-1])):
		if abs(mags[i+1] - mags[i])/(radii[i+1]-radii[i])<0.012: 
			r_edge = radii[i+1]
			break

	for i in range(len(mags[:-1])):
		if abs((mags[i+1] - mags[i])/(radii[i+1] - radii[i]))<0.1:
			r1 = radii[i]
			break

	r2 = 1.01*r_edge
	index1  = len(radii[radii<r1]) - 1
	index2  = len(radii[radii<r2]) - 1	
	plot_growth_curve(index1,index2,radii, mags,mag_auto, fluxes, sbs)

	return index1, index2

def cal_sb(ID, ra, dec, cat_type):
	catalog = Table.read('sex_fits/'+cat_type+'.'+str(ID)+'.fits')

	radii = np.array([2,4,6, 8, 10,12,14,16,18, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96,110,130,150])*0.187/2  # in arcsec
	mags = catalog[0]['MAG_APER']
	mag_auto = catalog[0]['MAG_AUTO']
	fluxes = catalog[0]['FLUX_APER']
	
	index1, index2 = growth_curve(str(ID), cat_type, radii, mags, fluxes,  mag_auto)
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
		run_sex(str(ID), str(cat_type))
		ra,dec = obj['RA'],obj['DEC']
		sb = cal_sb(ID, ra, dec, cat_type)
		sb_column.append(sb)
		print cat_type, ID,sb

	#write_table(sb_column, cat_type)
	plot_hist(np.array(sb_column),cat_type)