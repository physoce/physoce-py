import unittest
import numpy as np
from physoce import util

class UtilTestCase(unittest.TestCase):
    '''
    Tests for functions in physoce.util
    '''
    def test_compasstransform(self):
        theta_input = [0,30,-30,280,115,180]
        theta_check = [90,60,120,170,335,270]
        theta_out = util.compasstransform(theta_input)
        test = np.isclose(theta_out,theta_check)
        self.assertTrue(test.all(),msg='compasstransform failed')
        
if __name__ == '__main__':
    unittest.main()

