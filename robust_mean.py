import  numpy as np 

'''
  INPUT: 
      x: 1D numpy array
'''
def robust_mean(x,iter=1,show_step=False):
	for i in range(iter):
		mean = np.mean(x)
		std = np.std(x)
		x = x[x<mean+3*std]
		x = x[x>mean-3*std]
		if show_step == True:
			print 'Iteration:',i+1,' Mean:',round(mean,4),' Std:',round(std,4)
	return np.mean(x),np.std(x)