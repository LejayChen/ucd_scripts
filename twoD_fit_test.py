import numpy as np
from scipy.optimize import curve_fit
from astropy.table import *
from math import *
from scipy.stats import binned_statistic_2d as stats_2d
import matplotlib.pyplot as plt 
from matplotlib.colors import Colormap
from func import *

spacing = 4
fit_min =15
fit_max = 210
np.set_printoptions(precision=2)

def mask(xx, yy, ra_dis_dEN, dec_dis_dEN, original_data,method='mask'):
	xx_m = np.copy(xx)
	yy_m = np.copy(yy)
	original_data_m = np.copy(original_data)
	for i in range(len(ra_dis_dEN)):
		if abs(ra_dis_dEN[i])>fit_max or  abs(dec_dis_dEN[i])>fit_max:
			continue
		index_ra = int((ra_dis_dEN[i] + fit_max - spacing/2)/spacing)
		index_dec = int((dec_dis_dEN[i] + fit_max - spacing/2)/spacing)

		if method == 'mask':
			xx_m[index_ra-1:index_ra+2,index_dec-1:index_dec+2], yy_m[index_ra-1:index_ra+2,index_dec-1:index_dec+2] = np.NaN, np.NaN
			original_data_m[index_ra-1:index_ra+2,index_dec-1:index_dec+2] = np.NaN
		elif method == 'smooth':
			sf = 8 #smooth factor
			original_data_m[index_ra-1:index_ra+2,index_dec-1:index_dec+2] = np.mean(original_data[max(index_ra-sf,0):index_ra+sf+1,max(index_dec-sf,0):index_dec+sf+1])
		else:
			raise KeyError(str(method)+' is not a valid argument!')

	return xx_m, yy_m, original_data_m

def fit(bin_stats,xx,yy):
	for i in range(len(x)):
		for j in range(len(y)):
			if sqrt(x[i]**2+y[j]**2)<fit_min or sqrt(x[i]**2+y[j]**2)>fit_max:
				xx[i,j],yy[i,j],bin_stats[i,j] =  np.NaN,np.NaN,np.NaN
	xx_ravel=xx.ravel()[~np.isnan(xx.ravel())]
	yy_ravel=yy.ravel()[~np.isnan(yy.ravel())]
	bin_stats_ravel=bin_stats.ravel()[~np.isnan(bin_stats.ravel())]

	# Guess intial parameters
	n = -2
	PA = 1.2/4*pi
	e = 0.751
	H = 100.
	C = -0.2
	initial_guess = [n, PA, e, H, C]

	lower = [-np.inf,-np.inf, 0, 0,-10]
	upper = [0,np.inf, 1, 30000,10]
	bounds = [lower, upper]

	pred_params, uncert_cov = curve_fit(func, (xx_ravel, yy_ravel), bin_stats_ravel, p0=initial_guess, bounds=bounds, method='trf')

	perr = np.sqrt(np.diag(uncert_cov))
	err =err_func((xx_ravel, yy_ravel), bin_stats_ravel ,pred_params[0],pred_params[1],pred_params[2],pred_params[3],pred_params[4])

	n1 = -1.85088530618
	H1 = exp(5.10118386605)
	err_1D = err_func_1D((xx_ravel,yy_ravel),bin_stats_ravel,n1,H1)

	n2 = -1.83080934889
	H2 = exp(5.01768186697)
	err_1D_2 = err_func_1D((xx_ravel,yy_ravel),bin_stats_ravel,n2,H2)

	n3 = -1.84539280985
	H3 = exp(5.07330399332)
	err_1D_3 = err_func_1D((xx_ravel,yy_ravel),bin_stats_ravel,n3,H3)

	print '2D fitting error:',err,
	print '1D fitting error',err_1D,err_1D_2,err_1D_3
	print 'parameter standard error 2D:',perr

	return pred_params

M87 = (187.70583, 12.39111)
DIS = 17.21*1000 # in Kpc (luminosity distance from Mei et.al. 2011)
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]
cat_gc = cat_gc[np.sqrt((cat_gc['ra']-M87[0])**2+(cat_gc['dec']-M87[1])**2) /180.*pi*DIS<fit_max]
cat_gc = cat_gc[np.sqrt((cat_gc['ra']-M87[0])**2+(cat_gc['dec']-M87[1])**2) /180.*pi*DIS>fit_min]

cat_dEN = Table.read('pp.gal.nuc2.s.master.new.fits')
ra_dEN = cat_dEN['RA']
dec_dEN = cat_dEN['DEC'] 
ra_dis_dEN = np.array((ra_dEN-M87[0]) /180.*pi*DIS) # in kpc
dec_dis_dEN = np.array((dec_dEN-M87[1]) /180.*pi*DIS) # in kpc

ra = cat_gc['ra']
dec = cat_gc['dec']
dis_ra = np.array((ra-M87[0]) /180.*pi*DIS) # in kpc
dis_dec = np.array((dec-M87[1]) /180.*pi*DIS) # in kpc

#===============data===================
bin_ra = np.arange(-fit_max,fit_max+.1,spacing)
bin_dec = np.arange(-fit_max,fit_max+.1,spacing)
x = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
y = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
xx, yy = np.meshgrid(x, y)

ret = stats_2d(dis_ra,dis_dec,None,'count',bins=[bin_ra,bin_dec])
values = ret.statistic/(spacing**2)
values[np.isnan(values)]=0.0

#===========mask dE,N regions=============================
xx_m, yy_m, values_m = mask(xx, yy, ra_dis_dEN, dec_dis_dEN, values,method='mask')  

#==========plot data======================================
plt.hist2d(dis_ra,dis_dec,bins=[bin_ra,bin_dec])
plt.colorbar()
plt.show()

#======fit data======
print 'fitting function: I = H*np.sqrt(x2**2+(y2/e)**2)**n + C'
pred_params = fit(values_m,xx_m,yy_m)
print 'Parameters 2D fit [n,PA,e,H,C]:',pred_params

#==============plot fit result (residue)=======================
xx[xx==0] = 0.001
yy[yy==0] = 0.001

func_value = func((xx,yy),pred_params[0],pred_params[1],pred_params[2],pred_params[3],pred_params[4])

residual = func_value - values
for i in range(len(x)):
	for j in range(len(y)):
		if sqrt(x[i]**2+y[j]**2)<fit_min or sqrt(x[i]**2+y[j]**2)>fit_max:
			residual[i,j]=  0.0

plt.pcolor(xx, yy, residual)
plt.colorbar()
plt.xlim([-fit_max,fit_max])
plt.ylim([-fit_max,fit_max])
plt.show()