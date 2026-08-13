import unittest
from bioconcrete.boundary_3d import BoundaryCondition3D, validate_boundaries

class Boundary3DTests(unittest.TestCase):
    def test_blind_rejects_outlet(self):
        with self.assertRaises(ValueError): validate_boundaries({'x_max':BoundaryCondition3D('outlet')},'blind_crack')
    def test_through_requires_explicit_faces(self):
        with self.assertRaises(ValueError): validate_boundaries({},'through_crack')
        validate_boundaries({'x_min':BoundaryCondition3D('inlet'),'x_max':BoundaryCondition3D('outlet')},'through_crack')
    def test_negative_robin_coefficient_rejected(self):
        with self.assertRaises(ValueError): BoundaryCondition3D('robin',1,-1)

if __name__ == '__main__': unittest.main()
