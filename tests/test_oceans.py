import unittest
import numpy as np
from physoce import oceans

class OceansTestCase(unittest.TestCase):
    '''
    Tests for functions in physoce.oceans
    '''
    def test_wavedisp1(self):
        ### Test wavedisp ###
        # check just one value
        # values from original Matlab function
        mat_omega = 0.8976
        mat_k = 0.0990
        mat_Cph = 9.0631
        mat_Cg = 6.5488
        omega,k,Cph,Cg = oceans.wavedisp(7,12)
        test = np.isclose(np.array([omega,k,Cph,Cg]),
                                   np.array([mat_omega,mat_k,mat_Cph,mat_Cg]),
                                            atol = 1e-4)     
        self.assertTrue(test.all(),msg='wavedisp test #1 (single value): failed')         

    def test_wavedisp2(self):
        # test multiple values
        # values from original Matlab function
        mat_omega = np.array([0.8976,0.6283])
        mat_k = np.array([0.0990,0.0630])
        mat_Cph = np.array([9.0631,9.9667])
        mat_Cg = np.array([6.5488,8.4740])
        omega,k,Cph,Cg = oceans.wavedisp([7,10],12)
        test = np.isclose(np.array([omega,k,Cph,Cg]),
                          np.array([mat_omega,mat_k,mat_Cph,mat_Cg]),
                              atol = 1e-4)           
        self.assertTrue(test.all(),msg='wavedisp test #2 (multiple values): failed')
        
    def test_ustokes(self):
        ### Test ustokes ###
        # values from original Matlab function
        ust0_mat = 0.0545
        ust,zst = oceans.ustokes(2,7,12)
        ust0 = ust[0]
        test = np.isclose(ust0_mat,ust0,atol=1e-4)
        self.assertTrue(test,msg='ustokes test: failed')
        
    def test_bstream(self):
        ### Test bottom streaming ###
        # values from original Matlab function
        mat_UoK = np.array([0.0058,3.3409e-04])
        UoK = oceans.ubstream(np.array([2.,1.]),np.array([7.,10.]),12,'K')
        test = np.isclose(mat_UoK,
                      UoK, atol = 1e-4)             
        self.assertTrue(test.all(),msg='bstream test: failed')
        
if __name__ == '__main__':
    unittest.main()

