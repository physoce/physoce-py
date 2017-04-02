import unittest
import numpy as np
from physoce import tseries
'''
Author: Patrick Daniel
Date: 2016-07-14
First attempt to write and use a unit test. 
Specifically for the principal axis code, tseries.princax and therotation code, 
tseries.rot
Using https://jeffknupp.com/blog/2013/12/09/improve-your-python-understanding-unit-testing/ as a guideline

Updated: Tom Connolly
Date: 2017-04-02
Added tests for new depthavg function
'''

class TseriesTestCase(unittest.TestCase):
    '''
    Tests for 'tseries.princax() and tseries.rot()
    '''
    def test_for_princax_horizontal_line(self):
        '''
        Identify the the principal axis along a horizonal line, y = 2, which
        should be 0 degrees
        '''
        u = np.linspace(-3,3,50)
        v = np.ones(50)*2
        theta,major,minor = tseries.princax(u,v)
        self.assertEqual(theta, 0,msg='That is not the major axis')
    def test_for_princax_vertical_line(self):
        '''
        Identify the the principal axis along a vertical line, x = 2, which
        should be 90 degrees
        '''
        u = np.ones(50) * 2
        v = np.linspace(-3,3,50)
        theta,major,minor = tseries.princax(u,v)
        self.assertEqual(theta,90,msg='That is not the major axis: ' + str(theta))
    
    def test_for_rot_ninety_degrees(self):
        '''
        Test if a horizontal line is correctly rotated 90 degrees using tseries.rot        
        Modified version of in situ test from tseries.py
        '''        
        u = np.linspace(-3,3,4)
        v = np.zeros(4)
        ur,vr = tseries.rot(u,v,90)
        test1 = np.isclose(ur,v)
        test2 = np.isclose(vr,u)
        test_all = np.array([test1.all(),test2.all()])
        self.assertTrue(test_all.all(),msg='They are not equal!')

    def test_for_depthavg_1d(self):
        '''
        Test for correct answers with 1D arrays input as lists
        '''    
        
        x = [1,2,2,4]
        z = [-4,-3,-2,-1]
        h = 6

        # Test #1
        xda = tseries.depthavg(x,z,h,bottom='extrap',surface='mixed')
        ans = np.trapz([-1,1,2,2,4,4],[-6,-4,-3,-2,-1,0])/6
        test1 = (xda == ans)
        
        # Test #2
        xda = tseries.depthavg(x,z,h,bottom='zero',surface='extrap')
        ans = np.trapz([0,1,2,2,4,6],[-6,-4,-3,-2,-1,0])/6
        test2 = (xda == ans)
        
        # Test #3: optional ssh input
        xda = tseries.depthavg(x,z,h,ssh=1,bottom='mixed',surface='extrap')
        ans = np.trapz([1,1,2,2,4,8],[-6,-4,-3,-2,-1,1])/7
        test3 = (xda == ans)
        
        # Test #4: NaNs near surface
        x = [1,2,2,4,np.nan]
        z = [-4,-3,-2,-1,0]
        xda = tseries.depthavg(x,z,h,ssh=1,bottom='mixed',surface='extrap')
        ans = np.trapz([1,1,2,2,4,8],[-6,-4,-3,-2,-1,1])/7
        test4 = (xda == ans)
        
        test_all = np.array([test1,test2,test3,test4]).all()
        self.assertTrue(test_all,msg='depthavg 1D tests failed')
        
    def test_for_depthavg_1d(self):
        '''
        Test for correct answers with 2D array
        '''    
        
        x = np.array([[1,2,2,4],
                      [1,2,2,np.nan],
                      [1,2,np.nan,np.nan],
                      [np.nan,np.nan,np.nan,np.nan],
                      [1,2,2,2]])
        z = [-4,-3,-2,-1]
        ssh = np.array([0.,0.,0.,0.,1.])
        h = 6

        
        xda = tseries.depthavg(x,z,h,ssh,bottom='extrap',surface='extrap')
        
        # Test #1
        ans = np.trapz([-1,1,2,2,4,6],[-6,-4,-3,-2,-1,0])/6
        test1 = (xda[0] == ans)
        
        # Test #2
        ans = np.trapz([-1,1,2,2,2],[-6,-4,-3,-2,0])/6
        test2 = (xda[1] == ans)
        
        # Test #3
        ans = np.trapz([-1,1,2,5],[-6,-4,-3,0])/6
        test3 = (xda[2] == ans)
        
        # Test #4
        ans = np.nan
        test4 = np.isnan(xda[3])
        
        # Test #5
        ans = np.trapz([-1,1,2,2,2,2],[-6,-4,-3,-2,-1,1])/7
        test5 = np.isclose(xda[4],ans)

        test_all = np.array([test1,test2,test3,test4,test5]).all()
        self.assertTrue(test_all,msg='depthavg 1D tests failed')        
        
if __name__ == '__main__':
    unittest.main()
    