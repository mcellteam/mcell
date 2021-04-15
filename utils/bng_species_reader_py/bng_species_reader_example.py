import sys
import os
import pandas as pd

MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
    sys.path.append(os.path.join(MCELL_PATH, 'lib'))
else:
    print("Error: variable MCELL_PATH that is used to find the mcell library was not set.")
    sys.exit(1)

import mcell as m

# utility module to load ASCII viz_output .dat file from 
import species_reader


def get_bound_em(complex, src_em, comp_name):
    # assuming that there is a single component in src_em 
    # (source elementary molecule) with comp_name that is bound
    
    bond_index = -1
    for c in src_em.components:
        if c.component_type.name == comp_name:
            # not allowed cases, we need the number
            assert c.bond != m.BOND_UNBOUND  # no bond
            assert c.bond != m.BOND_ANY      # !?
            assert c.bond != m.BOND_BOUND    # !+
            bond_index = c.bond 
    
    assert bond_index != -1, "Did not find " + comp_name + " in " + src_em.to_bngl_str()
    
    # find this bond in the complex (might be slow)
    for em in complex.elementary_molecules:
        if em is src_em:
            continue
        
        for c in em.components:
            if c.bond == bond_index:
                return em
    
    assert False, "Did not find paired bond " + str(bond_index) + " in " + complex.to_bngl_str()
    

def convert_species_file(file_name):
    # read the .species file and parse complexes to the internal MCell
    # representation 
    complex_counts = species_reader.read_species_file(file_name)
    
    # prepare the component type that we will append to the CaMKII molecules
    # the new component means 'upper' and has a boolean state
    ct_u = m.ComponentType('u', ['0','1'])
    
    # process the data by adding component 'u' to the CaMKII elementary molecules 
    for (complex, count) in complex_counts:
        
        # if this a CaMKII dodecamer? 
        camkiis = []
        for em in complex.elementary_molecules:
            em_name = em.elementary_molecule_type.name
            if em_name == 'CaMKII':
                camkiis.append(em)
                
        if len(camkiis) != 12:
            # output it directly
            print(complex.to_bngl_str() + " " + str(count))
            continue
        
        # ok, we have the holoenzyme, how can we figure out which 
        # of the CaMKIIs belong to the upper and the lower 
        # ring?
        
        # let's say that the first one is in the upper ring,
        # go along the 'r' 
        upper_first = camkiis[0]
        upper_ring = [ upper_first ]
        curr = upper_first
        for i in range(1, 6):
            curr = get_bound_em(complex, curr, 'r')
            upper_ring.append(curr)
        
        assert upper_first is get_bound_em(complex, curr, 'r'), "A ring must be formed"
        
        # then go down along 'c' and go again along 'r' to get molecules of the lower ring
        lower_first = get_bound_em(complex, camkiis[0], 'c')
        lower_ring = [ lower_first ]
        curr = lower_first
        for i in range(1, 6):
            curr = get_bound_em(complex, curr, 'r')
            lower_ring.append(curr)
            
        assert lower_first is get_bound_em(complex, curr, 'r'), "A ring must be formed"

        # now the modifications - add components by instatiating the component type 'u'
        # the way how the complexes were is parsed was that each complex has its own instance 
        # of the elementary molecule type for CaMKII so lets change one of them
        upper_first.elementary_molecule_type.components.append(ct_u)
        
        for em in upper_ring:
            em.components.append(ct_u.inst('1'))
        
        for em in lower_ring:
            em.components.append(ct_u.inst('0'))
        
        print(complex.to_bngl_str() + " " + str(count))
        
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Expected .species file name as argument")
        sys.exit(1)
        
    convert_species_file(sys.argv[1])
    