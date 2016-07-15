import unittest
import numpy as np
import tseries
'''
Author: Patrick Daniel
Date: 2016-07-14
First attempt to write and use a unit test. 
Specifically for the principal axis code, tseries.princax and therotation code, 
tseries.rot
Using https://jeffknupp.com/blog/2013/12/09/improve-your-python-understanding-unit-testing/ as a guideline
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
    
    def test_for_rot_ninty_degrees(self):
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
        
if __name__ == '__main__':
    unittest.main()
    