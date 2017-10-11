""" Forward-mapping an mbuild compound

1) Specify a mapping procedure
    a) Provide mapping files in .json format
    b) Provide heuristic like 3:1
2) Provide a topology of molecules (so the system has knowledge of
    which atoms group together to form a molecule)
    a) Gromacs already has molecule information
    b) Hoomd/lammps do not have molecule information
    c) If given heuristic information then start grouping together atoms 
        independent of molecular topology 
3) Do converting
    a) Molecule by molecule, convert the corresponding atoms to beads
    b) If given heuristic, each 3:1 grouping becomes its own molecule
4) Create bonds
    a) If given mapping information, then bonding patterns will already be 
        specified
    b) If given heuristics, then assume linear bonding pattern?

"""
import numpy as np
import mbuild as mb
import mdtraj

def coarse_grain(fine_grained, heuristic=None, mapping_file=None):
    # Molecular conversions should be a dictionary whose keys are cg beads
    # and indices correspond to the global atom indices
    if not heuristic:
        molecular_conversions = _apply_heuristic(fine_grain, heuristic)
    elif not mapping_file:
        molecule_info = _extract_molecules(fine_grained)
        molecular_conversions = _apply_mapping(fine_grain, mapping_file)
    else:
        return fine_grained
    coarse_grained = _convert(fine_grained, molecular_conversions)


def _apply_heuristic(fine_grained, heuristic):
    """ Apply a heuristic to outline the molecular conversions
    i.e. 3:1 means this method will group 3 heavy atoms to 1 CG bead"""
    return None

def _apply_mapping(fine_grained, mapping_file):
    """ Read a mapping.json file to outline molecular conversions


    """
    return None

def _extract_molecules(fine_grained):
    """ From the structure file,
    be able to determine what molecule or residue each atom belongs to
    This should be easy with gro files,
    With hoomd/lammps files, there are no molecules/residues so this 
    has to be provided or somehow inferred
    """

    return None


def _convert(fine_grained, molecular_conversions):
    """ Use molecular_converions to fwd map the fine_grained

    """
    for child in fine_grained.children:
        # check if it can map
        if child in molecular_conversions:
            # newbead = center of mass of these corresponding indices
            return 0
        # If not, convert the children
        else: 
            _convert(child, molecular_conversions)

