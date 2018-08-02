import unittest

import numpy as np

from lib.fftso3 import vee, hat

class TestMiscMatrixFunc(unittest.TestCase):
    def test_hat(self):
        vec = [1, 2, 3]
        exp_mat = np.array([[0, -3, 2],[3, 0, -1],[-2, 1, 0]])

        np.testing.assert_array_equal(hat(vec), exp_mat)

    def test_vee(self):
        mat = np.array([[0, -3, 2],[3, 0, -1],[-2, 1, 0]])
        exp_vec = [1, 2, 3]

        np.testing.assert_array_equal(vee(mat), exp_vec)
