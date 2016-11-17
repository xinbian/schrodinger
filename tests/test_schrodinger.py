#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_schrodinger
----------------------------------

Tests for `schrodinger` module.
"""


import sys
import unittest
import numpy as np


from schrodinger.schrodinger import *


class TestSchrodinger(unittest.TestCase):
    period=2
    def setUp(self):    
        #generate wave function
        xx=np.linspace(-self.period/2,self.period/2,10)
        yy=np.sin(2*np.pi*xx/self.period)
        return xx, yy    
    def tearDown(self):
        pass
    #this test function isn't testing the code, but testing the setUpfunction in the unittest
    def test_unit(self):
        xx, yy=self.setUp()
        #the area of the period function should be 0
        self.assertTrue(abs(yy.sum())<0.00001)
    #test wave_cf function for Fourier series
    def test_wave_for_coeffcient(self):
        xx, yy=self.setUp()
        #if given a sin(x), the coeffcient for e^0 should be 0
        self.assertTrue(abs(wave_cf(xx,yy,1,10,self.period)[0]),0.00001)
    #for Fourier hamilton matrix, the non-diagonal elements are 0, test sum of them should be 0
    def test_hamilton_for(self):
        hij=hamilton_matrix(1, 10)
        self.assertEqual(hij.sum()-np.matrix.trace(hij),0j)
    #test Laplacian on coeffcients
    def test_laplacian(self):
        xx=np.linspace(-1,1,20)
        yy=xx**2
        #the first coeffcient should be -2 due to given function
        self.assertTrue(abs(after_cf(xx,yy,2,10,2)[0]+2)<0.00001)
    #test norm matrix for Fourier
    def test_norm_for(self):
        self.assertEqual(np.real(norm_matrix(1,10,2)[5][5]),2)
        self.assertEqual(np.real(norm_matrix(1,10,2)[1][2]),0)
    
    

if __name__ == '__main__':
    sys.exit(unittest.main())