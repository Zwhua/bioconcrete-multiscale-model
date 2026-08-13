import unittest
from bioconcrete.validation_3d import dimensional_reduction_checks

class DimensionalReduction3DTests(unittest.TestCase):
    def test_dual_adapters(self):
        checks=dimensional_reduction_checks()
        for value in checks.values():
            if isinstance(value,dict): self.assertTrue(value['passed'],value)

if __name__=='__main__': unittest.main()
