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

def ubwave(Hsig,wavper,h):
    """
Ub = ubwave(Hsig,wavper,h)

Calculate bottom wave orbital velocity 

Inputs: 
Hsig: significant wave height [m]
per: wave period [s]
h: bottom depth [m]

Output:
Ub: wave orbital velocity near the bottom [m/s]
    """
    
    a = Hsig/np.sqrt(8)    # wave amplitude
    (omega,k,Cph,Cg) = wavedisp(wavper,h)
    Ub = a*2*np.pi*wavper**-1.*(np.sinh(k*h))**-1
    return Ub

def ubstream(Hsig,wavper,h,formula='LH',kN=0.03):
    """
Uo = ubstream(Hsig,wavper,h,kN,formula='LH')
    
Calculate bottom streaming velocity (wave-averaged Eulerian velocity at the top
of the wave boundary layer).
    
Inputs: 
Hsig: significant wave height [m]
per: wave period [s]
h: bottom depth [m]
formula: 'LH' - Longuet-Higgins (1953) [default]
         'K' - Kranenburg (2012)
kN: Nikuradse roughness length (kN = 30*zo)
    [default = 0.03m, not needed for Longuet-Higgins formula]

Output:
Uo = velocity at top of wave boundary layer [m/s]

References:
Longuet-Higgins, M. S. (1953) Mass transport in water waves, Philos. Trans. R.
    Soc. London, Ser. A, 245(903), 535-581, doi:10.1098/rsta.1953.0006.
Kranenburg, W. M. et al. (2012) Net currents in the wave bottom boundary layer:
    On waveshape streaming and progressive wave streaming, 117, F03005, 
    doi:10.1029/2011JF002070.
    """

    (omega,k,Cph,Cg) = wavedisp(wavper,h)
    Ub = ubwave(Hsig,wavper,h) 
    Ab = Ub*(2*np.pi*wavper**-1)**-1 # near bottom wave excursion
    if formula == 'LH':
        fac = 0.75
    elif formula == 'K':
        fac = 0.345 + 0.7*(Ab/kN)**-0.9 - 0.25*(np.sinh(k*h))**-2
    else:
        raise ValueError('specified formula not understood')
        
    Uo = fac*Ub**2*(Cph**-1)
    return Uo

if __name__ == '__main__':
    
    ### Test wavedisp ###
    # check just one value
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
    
    # test multiple values
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
        
    ### Test bottom streaming ###
    # values from original Matlab function
    mat_UoK = np.array([0.0058,3.3409e-04])
    UoK = ubstream(np.array([2.,1.]),np.array([7.,10.]),12,'K')
    test = np.isclose(mat_UoK,
                      UoK, atol = 1e-4)             
    if test.all():
        print('ubstream test: passed')
    else:
        raise ValueError('ubstream test: failed')