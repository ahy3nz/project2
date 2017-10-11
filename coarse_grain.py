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
    d) Should return some sort of list that relates coarse grain beads 
        to global atomic indices (regardless of heuristic or mapping file)
3) Do converting
    a) Molecule by molecule, convert the corresponding atoms to beads
4) Create bonds
    a) If given mapping information, then bonding patterns will already be 
        specified
    b) If given heuristics, then assume linear bonding pattern?

"""
import pdb
import numpy as np
import mbuild as mb
import mdtraj
import json
from parmed.periodic_table import Mass
from collections import OrderedDict


def coarse_grain(fine_grained, heuristic=None, mapping_files=[], 
        table_of_contents=None):
    """
     Molecular conversions should be a dictionary whose keys are cg beads
     and indices correspond to the global atom indices
     Table of contents matches molecules to global atom indices
    """
    # Heuristic denotes something like general 3:1 heavy atom mapping
    if heuristic:
        molecular_conversions, cg_bonds = _apply_heuristic(fine_grained, heuristic)

    # Mapping files denote a list of json files to create a large dictionary
    elif mapping_files:

        # Load in the mapping files into mapping_info
        mapping_info = {}
        for json_file in mapping_files:
            print("Loading json mapping <{}>".format(json_file))
            mapping_info.update(json.load(open(json_file,'r'), object_pairs_hook=OrderedDict))

        # Create a table of contents to match molecuels to global atom indices
        if not table_of_contents:
            print("Generating table of contents")
            table_of_contents = _extract_molecules(fine_grained, mapping_info)
        else:
            print("Loading table of contents <{table_of_contents}>".format(**locals()))
            table_of_contents = json.load(open(table_of_contents,'r'))


        # Combine table of contents and mapping info
        # Generate dictionary whose keys are cg beads and values are global 
        # atom indices
        # And generate a list of cg bonds
        molecular_conversions, cg_bonds = _apply_mapping( 
                table_of_contents, mapping_info)
    else:
        print("More coarse-graining parameters necessary, returning original compound")
        return fine_grained

    coarse_grained = _convert_xyz(fine_grained, molecular_conversions)
    # Apply bonding
    coarse_grained_bonded = _apply_bonding(coarse_grained, cg_bonds)

    # Save
    coarse_grained.save('cg.mol2',overwrite=True)

def _apply_bonding(coarse_grained, cg_bonds):
    """ Apply bonds between CG beads

    Parameters
    ---------
    coarse_grained : mb.Compound
    cg_bonds : list of tuples (2, n_beads)

    Returns
    -------
    cg_bonded : mb.Compound
    """
    for bead_i, bead_j in cg_bonds:
        pdb.set_trace()
        coarse_grained.add_bond([coarse_grained.children[bead_i], 
                                coarse_grained.children[bead_j]])
    return coarse_grained

def _convert_xyz(fine_grained, molecular_conversions):
    """ Use molecular_converions to fwd map the fine_grained

    Parameters
    ----------
    molecular_conversions : list of lists
        [CG_bead, [global atom indices]]
        index of list corresponds to bead index

    """
    cg_system = mb.Compound()
    for bead, atom_indices in molecular_conversions:
        new_bead = mb.Compound(name=bead.split("-")[0], pos=_compute_center_of_mass(fine_grained, atom_indices))
        cg_system.add(new_bead)

    return cg_system


    

def _apply_mapping(table_of_contents, mapping_info):
    """ Read mapping.json files to outline molecular conversions

    table_of_contents : dict
        {resname : [global_atom_indices]}
    mapping_info : dict
        loaded-in json file

    Returns
    -------
    molecular_conversions : list of lists
        [CG_bead, [global atom indices]]
        index of list corresponds to bead index
    """
    molecular_conversions = []
    cg_bonds = []
    # Using mapping_info's local atom indices,
    # Shift them by the global atom indices to 
    # Add the CG bead to a list of global atom indices
    for molecule in table_of_contents.keys():
        molecule_name = molecule.split("-")[0]
        global_atom_indices = table_of_contents[molecule]
        first_bead_index = len(molecular_conversions)

        for bead in mapping_info[molecule_name]['map'].keys():
            updated_atom_indices = [local_index + global_atom_indices[0] for local_index in mapping_info[molecule_name]['map'][bead]['children']]
            molecular_conversions.append([bead,updated_atom_indices])

        # Update bonding for global BEAD indices
        for bead_i, bead_j in mapping_info[molecule_name]['bond']:
            cg_bonds.append([bead_i + first_bead_index, 
                        bead_j + first_bead_index])



    return molecular_conversions, cg_bonds


def _extract_molecules(fine_grained, mapping_info):
    """ From the structure file,
    be able to determine what molecule or residue each atom belongs to
    This should be easy with gro files,
    With hoomd/lammps files, there are no molecules/residues so this 
    has to be provided or somehow inferred

    Currently this is being implemented as a dict
    """
    has_res_information = True
    residues = [item for item in mapping_info.keys()]
    table_of_contents = OrderedDict()
    if has_res_information:
        traj = mdtraj.load('lots.gro')
        for residue in traj.topology.residues:
            indices = [at.index for at in residue.atoms]
            table_of_contents.update({residue.name + "-" + str(residue.index): indices})

    return table_of_contents


def _apply_heuristic(fine_grained, heuristic):
    """ Apply a heuristic to outline the molecular conversions
    i.e. 3:1 means this method will group 3 heavy atoms to 1 CG bead"""
    return None

def _compute_center_of_mass(fine_grained, atom_indices):
    """ Compute center of mass"""
    masses = [Mass[fine_grained.children[index].name] for index in atom_indices]
    total_mass = sum(masses)
    com = np.ndarray(3)
    for i in range(3):
        com[i] = sum([fine_grained.children[index].pos[i]*mass/total_mass for index,mass in zip(atom_indices,masses)])
    return com 



if __name__ == "__main__":
    #fine_grained = mb.load('propane.mol2')
    fine_grained = mb.load('two_propane.gro')
    fine_grained.name="PR3"
    #toc = None
    toc = "two_propane_toc.json"
    coarse_grained = coarse_grain(fine_grained, mapping_files=["propane.json"], 
            table_of_contents=toc)

