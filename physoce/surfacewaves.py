""" 
Functions for analysis of surface gravity waves
Tom Connolly 2014 
"""

import numpy as np

def wavedisp(wavper,h):
	""" (omega,k,Cph,Cg) = wavedisp(wavper,h)
	------------------
	Returns [omega,k,Cph,Cg]
	
	Inputs (can use arrays): 	
				wavper	- wave period
				h	- water depth
	Outputs: 	
				omega	- wave frequency
				k	- wave number	
				Cph	- phase speed
				Cg	- group velocity 
			
	T Connolly 2014
	based on Matlab function wavedisp.m from S Lentz """
	
	wavper=np.array([wavper])
	h=np.array([h])
	
	omega = 2*np.pi/wavper
	g = 9.81
	c = omega**2*h/g
	
	x = np.sqrt(c)	
	
	d = 100*np.ones(np.shape(wavper))
	tol = 5e-16
	while any([d>tol]):
		f1=x*np.tanh(x)-c
		f2=x*(1/np.cosh(x))**2+np.tanh(x)
		x=x-f1/f2
		d=np.abs((c-x*np.tanh(x))/c)
	k=x/h
	Cph=omega/k
	Cg=(g*k*(1/np.cosh(x))**2+g*np.tanh(x))/(2*np.sqrt(g*k*np.tanh(x)))
	
	return (omega,k,Cph,Cg)

def ustokes(Hsig,wavper,wdepth,zst=()):
	""" (ust,zst) = ustokes(Hsig,wavper,wdepth,zst=()):
	Stokes drift velocity
	---------------------
	Inputs: 	Hsig	- Significant wave height [m]
				wavper	- wave period [s]
				wdepth	- bottom depth [m]
				zst	- depths of calculation [m] 
					(optional, every 1m if not given)
	Outputs:	ust	- Stokes drift [m/s]
				zst	- depths of calculation
	
	T Connolly 2014
	based on Matlab function ustokes.m from S Lentz"""
	
	if not any(zst):
		zst = np.arange(-wdepth,0)
		
	(omega,wavnum,Cph,Cg)=wavedisp(wavper,wdepth)
	A=Hsig**2*omega*wavnum/(16*np.sinh(wavnum*wdepth)**2)
	ust=A*np.cosh(2*wavnum*(zst+wdepth))
		
	return (ust,zst)
