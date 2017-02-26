import numpy as np
from math import *

def func(xy, n, PA, e, H, C):
	x, y = xy
	x2 = x*cos(PA)-y*sin(PA)
	y2 = y*cos(PA)+x*sin(PA)
	I = H*np.sqrt(x2**2+(y2/e)**2)**n + C
	return I

def func_1D(xy,n,H):
	x, y = xy
	I = H*np.sqrt(x**2+y**2)**n
	return I

def err_func_1D(xy,values,n, H):
	err = 0
	x,y=xy
	for i in range(len(x)):
		err = err + (func_1D((x[i],y[i]), n, H) - values[i])**2
	return sqrt(err)

def func_gaussian(xy,PA, sigmax, sigmay, H):
	x, y = xy
	theta = PA
	x2 = x*cos(theta)-y*sin(theta)
	y2 = y*cos(theta)+x*sin(theta)
	I = H*np.exp(-x2**2/sigmax)*np.exp(-y2**2/sigmay)

	return I

def err_func(xy,values,n, PA, e, H, C):
	err = 0
	x,y=xy
	for i in range(len(x)):
		err = err + (func((x[i],y[i]), n, PA, e, H, C) - values[i])**2
	return sqrt(err)

def data_gen(func_name,xx,yy):
	n = -1.0
	PA = 1.2/4*pi
	e = 0.795
	H = 160.0
	xy = (xx,yy)

	I = func_name(xy,n,PA,e,H) 
	if len(I.shape)==1:
		data_gen_values = I + np.random.random(I.shape)/10.
	if len(I.shape)==2:
		data_gen_values = I + np.random.random([I.shape[0],I.shape[1]])/10.
	return data_gen_values