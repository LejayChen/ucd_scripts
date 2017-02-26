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
fit_max = 190

def fit(bin_stats):
	x = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
	y = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
	xx, yy = np.meshgrid(x, y)
	#print xx.shape,bin_stats.shape

	for i in range(len(x)):
		for j in range(len(y)):
			#print i,j
			if sqrt(x[i]**2+y[j]**2)<fit_min or sqrt(x[i]**2+y[j]**2)>fit_max:
				xx[i,j],yy[i,j],bin_stats[i,j]=  np.NaN,np.NaN,np.NaN
	xx_ravel=xx.ravel()[~np.isnan(xx.ravel())]
	yy_ravel=yy.ravel()[~np.isnan(yy.ravel())]
	bin_stats_ravel=bin_stats.ravel()[~np.isnan(bin_stats.ravel())]

	'''
	PA = 1.2/4*pi
	sigmax = 11
	sigmay = 11
	H = 300
	initial_guess = [PA, sigmax, sigmay, H]

	lower = [-pi,0, 0, 0]
	upper = [pi, np.inf,np.inf, 3000]
	bounds = [lower, upper]
	'''
	# Guess intial parameters
	n = -2
	PA = 1.2/4*pi
	e = 0.751
	H = 2800.
	C = -0.2
	initial_guess = [n, PA, e, H, C]

	lower = [-np.inf,-np.inf, 0, 0,-10]
	upper = [0,np.inf, 1, 30000,10]
	bounds = [lower, upper]

	pred_params, uncert_cov = curve_fit(func, (xx_ravel, yy_ravel), bin_stats_ravel, p0=initial_guess, bounds=bounds, method='trf')
	#pred_params, uncert_cov = curve_fit(func_gaussian, (xx_ravel, yy_ravel), bin_stats_ravel,p0=initial_guess,bounds=bounds,method='trf')

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

	print '2D fitting error:',err,'1D fitting',err_1D,err_1D_2,err_1D_3
	print 'parameter standard error 2D:',perr

	return pred_params

def fit_test(data_gen_value):
	x = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
	y = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
	xx, yy = np.meshgrid(x, y)
	for i in range(len(x)):
		for j in range(len(y)):
			if sqrt(x[i]**2+y[j]**2)<fit_min or sqrt(x[i]**2+y[j]**2)>fit_max:
				xx[i,j],yy[i,j],data_gen_value[i,j]=  np.NaN,np.NaN,np.NaN
	xx_ravel=xx.ravel()[~np.isnan(xx.ravel())]
	yy_ravel=yy.ravel()[~np.isnan(yy.ravel())]
	data_gen_ravel=data_gen_value.ravel()[~np.isnan(data_gen_value.ravel())]
	# Guess intial parameters
	n = -1.1
	PA = 1.2/4*pi
	e = 0.795
	H = 600
	initial_guess = [n, PA, e, H]

	lower = [-np.inf, 0.1, 0.1, 1]
	upper = [-0.1, np.inf, 1, 30000]
	bounds = [lower, upper]

	pred_params, uncert_cov = curve_fit(func, (xx_ravel, yy_ravel), data_gen_ravel,p0=initial_guess,bounds=bounds)
	err =err_func((xx_ravel, yy_ravel), data_gen_ravel ,pred_params[0],pred_params[1],pred_params[2],pred_params[3],pred_params[4])
	print uncert_cov,err
	return pred_params

M87 = (187.70583, 12.39111)
DIS = 17.21*1000 # in Kpc (luminosity distance from Mei et.al. 2011)
cat_gc = Table.read('ngvs_pilot_xdclass1.0_g18.0-25.0.fits')
cat_gc = cat_gc[cat_gc['p_gc']>0.95]
cat_gc = cat_gc[np.sqrt((cat_gc['ra']-M87[0])**2+(cat_gc['dec']-M87[1])**2) /180.*pi*DIS<fit_max+spacing]
#cat_gc = cat_gc[np.sqrt((cat_gc['ra']-M87[0])**2+(cat_gc['dec']-M87[1])**2) /180.*pi*DIS>fit_min]

ra = cat_gc['ra']
dec = cat_gc['dec']
dis_ra = np.array((ra-M87[0]) /180.*pi*DIS) # in kpc
dis_dec = np.array((dec-M87[1]) /180.*pi*DIS) # in kpc

#===============data===================
bin_ra = np.arange(-fit_max,fit_max+.1,spacing)
bin_dec = np.arange(-fit_max,fit_max+.1,spacing)
ret = stats_2d(dis_ra,dis_dec,None,'count',bins=[bin_ra,bin_dec])
values = ret.statistic/(spacing**2)  #2D statistic
values[np.isnan(values)]=0.0

#==========plot data======================================
plt.hist2d(dis_ra,dis_dec,bins=[bin_ra,bin_dec])
plt.colorbar()
plt.show()
#======fit data======
pred_params = fit(values)
print pred_params

#======================simulated data=======================
x = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
y = np.arange(-fit_max+spacing/2.,fit_max-spacing/2.+0.1,spacing)
xx, yy = np.meshgrid(x, y)
#data_gen_value = data_gen(func,xx,yy)
#pred_params = fit_test(data_gen_value)

#=================plot residual==============================
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