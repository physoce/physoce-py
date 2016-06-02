""" 
Functions for analysis of surface gravity waves
Tom Connolly 2014 
"""

import numpy as np

def wavedisp(wavper,h):
	""" 
(omega,k,Cph,Cg) = wavedisp(wavper,h)
------------------
Returns [omega,k,Cph,Cg]
	
Inputs (can use arrays): 	
			wavper	- wave period [s]
			h	- water depth [m]
   
Outputs: 	
			omega	- angular wave frequency [radians/s]
			k	- angular wave number	 [radians/m]
			Cph	- phase speed [m/s]
			Cg	- group velocity [m/s]
    """
			
	""" T Connolly 2014
	based on Matlab function wavedisp.m from S Lentz """
	
    # make sure inputs are arrays
	wavper=np.array(wavper)
	h=np.array(h)
	
	omega = 2*np.pi/wavper
	g = 9.8
	c = omega**2*h/g
	
	x = np.sqrt(c)	
	
   	d = 100*np.ones(np.shape(wavper))
	tol = 5.*np.finfo(float).eps
	while (d>tol).any():
		f1=x*np.tanh(x)-c
		f2=x*(1/np.cosh(x))**2+np.tanh(x)
		x=x-f1/f2
		d=np.abs((c-x*np.tanh(x))/c)
	k=x/h
	Cph=omega/k
	Cg=(g*x*(1/np.cosh(x))**2+g*np.tanh(x))/(2*np.sqrt(g*k*np.tanh(x)))
	
	return (omega,k,Cph,Cg)

def ustokes(Hsig,wavper,h,zst=()):
	""" 
ust,zst = ustokes(Hsig,wavper,h,zst=()):
Stokes drift magnitude
---------------------
Inputs:     Hsig	- significant wave height [m]
            wavper	- wave period [s]
            h     - bottom depth [m]
            zst	- depths of calculation [m] 
            (optional, every 1m if not given)
            
Outputs:	ust	- Stokes drift [m/s]
			zst	- depths of calculation [m]
    """
	
	"""T Connolly 2014
	based on Matlab function ustokes.m from S Lentz"""
	
	if not any(zst):
		zst = np.arange(0,-h-1,-1)
		
	(omega,k,Cph,Cg)=wavedisp(wavper,h)
	A=Hsig**2*omega*k/(16*np.sinh(k*h)**2)
	ust=A*np.cosh(2*k*(zst+h))
		
	return (ust,zst)
 
if __name__ == '__main__':
    
    ### Test wavedisp ###

    # check just one value instead of list
    mat_omega = 0.8976
    mat_k = 0.0990
    mat_Cph = 9.0631
    mat_Cg = 6.5488
    omega,k,Cph,Cg = wavedisp(7,12)
    test = np.isclose(np.array([omega,k,Cph,Cg]),
                      np.array([mat_omega,mat_k,mat_Cph,mat_Cg]),
                          atol = 1e-4)             
    if test.all():
        print('wavedisp test #1: passed')
    else:
        raise ValueError('wavedisp test #1: failed')    
    
    # values from original Matlab function
    mat_omega = np.array([0.8976,0.6283])
    mat_k = np.array([0.0990,0.0630])
    mat_Cph = np.array([9.0631,9.9667])
    mat_Cg = np.array([6.5488,8.4740])
    omega,k,Cph,Cg = wavedisp([7,10],12)
    test = np.isclose(np.array([omega,k,Cph,Cg]),
                      np.array([mat_omega,mat_k,mat_Cph,mat_Cg]),
                          atol = 1e-4)             
    if test.all():
        print('wavedisp test #2: passed')
    else:
        raise ValueError('wavedisp test #2: failed')
        
    ### Test ustokes ###
    # values from original Matlab function
    ust0_mat = 0.0545
    ust,zst = ustokes(2,7,12)
    ust0 = ust[0]
    test = np.isclose(ust0_mat,ust0,atol=1e-4)
    if test:
        print('ustokes test: passed')
    else:
        raise ValueError('ustokes test: failed')      