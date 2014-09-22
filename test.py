import unittest
import fft

class TestFFT(unittest.TestCase):
    def _check_result(self, expected, actual):
        self.assertEqual(len(expected), len(actual), 'Result length is not correct')
        for exp_pt, act_pt in zip(expected, actual):
            if abs(exp_pt[0] - act_pt[0]) > 1e-6:
                return False
            if abs(exp_pt[1] - act_pt[1]) > 1e-6:
                return False
        return True

    def test_shuffle(self):
        fin = [(0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (0, 0)]
        inverse = 0
        fout_exp = [
            (4.0, 0.0), (-1.5, -0.8660254037844388),
            (-0.5, -0.8660254037844386), (0.0, 0.0),
            (-0.5, 0.8660254037844388), (-1.5, 0.8660254037844386)
        ]
        
        fout_act = fft.fft(fin, inverse)
        good_result = self._check_result(fout_exp, fout_act)
        self.assertTrue(good_result, '%s vs %s' % (fout_exp, fout_act))


if __name__ == '__main__':
    unittest.main()
