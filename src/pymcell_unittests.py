import unittest
import pymcell as m

class SimpleSpeciesTestCase(unittest.TestCase):
    def setUp(self):
        self.vm1 = m.Species("vm1", 1e-6)
        self.vm1_down = m.OrientedSpecies(self.vm1, m.Orient.down)
        self.sm1 = m.Species("sm1", 1e-6, surface=True)
        self.sm1_up = m.OrientedSpecies(self.sm1, m.Orient.up)
        self.vm2 = m.Species("vm2", 1e-6)
        self.vm2_down = m.OrientedSpecies(self.vm2, m.Orient.down)

    def test_mol_is_3d(self):
        assert self.vm1.surface == False, "Volume molecule has 'surface' True"

    def test_mol_is_2d(self):
        assert self.sm1.surface == True, "Surface molecule has 'surface' False"


class SimpleReactionTestCase(SimpleSpeciesTestCase):
    def test_vm1_to_null(self):
        rxn = m.Reaction(self.vm1.down(), None, 1e8)

    def test_vm1_to_vm2(self):
        rxn = m.Reaction(self.vm1, self.vm2, 1e8)

    def test_vm1sm1_to_vm2(self):
        rxn = m.Reaction((self.vm1.down(), self.sm1.up()), self.vm2.down(), 1e8)


if __name__ == "__main__":
    unittest.main()
                            
