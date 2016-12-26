import  numpy as np 

def robust_mean(x,iter=1,show_step=False):
	'''
	Robust Statistics Function, cast away data point out of 3 sigmas each iteration
	Based on numpy

	INPUT: 
		x: 1D numpy array
		   the array to do statistics on
		iter: int 
		    number of iterations (iter=0 means normal np.mean())
		show_step: boolean
		    if Ture, will print mean and std of each step

	'''
	for i in range(iter):
		mean = np.mean(x)
		std = np.std(x)
		x = x[x<mean+3*std]
		x = x[x>mean-3*std]
		if show_step == True:
			print 'Iteration:',i+1,' Mean:',round(mean,4),' Std:',round(std,4)
	return np.mean(x),np.std(x)