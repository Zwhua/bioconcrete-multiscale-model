import unittest
from bioconcrete.cli import build_parser
class Cli3DTests(unittest.TestCase):
 def test_3d_options_parse(self):
  a=build_parser().parse_args(['simulate','--level','3d','--resume','--checkpoint','c.npz','--memory-limit-gb','2','--linear-solver','cg'])
  self.assertEqual(a.level,'3d');self.assertTrue(a.resume);self.assertEqual(a.linear_solver,'cg')
 def test_validate_3d_options_parse(self):
  a=build_parser().parse_args(['validate-3d','--output','out','--full'])
  self.assertEqual(a.command,'validate-3d');self.assertTrue(a.full)
 def test_render_3d_options_parse(self):
  a=build_parser().parse_args(['render-3d','--run','run'])
  self.assertEqual(a.command,'render-3d')
if __name__=='__main__':unittest.main()
