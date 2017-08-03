import unittest
import pymcell as m

class SimpleSpeciesTestCase(unittest.TestCase):
    def test_mol_is_3d(self):
        self.vm = m.Species("vm", 1e-6)
        assert self.vm.surface == False, "Volume molecule has 'surface' True"

    def test_mol_is_2d(self):
        self.sm = m.Species("sm", 1e-6, surface=True)
        assert self.sm.surface == True, "Surface molecule has 'surface' False"

if __name__ == "__main__":
    unittest.main()
                            
