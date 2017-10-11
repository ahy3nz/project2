import pdb
import numpy as np
import mbuild as mb
import mdtraj
import json
from parmed.periodic_table import Mass
from collections import OrderedDict

""" Forward-mapping an mbuild compound


1) Specify a mapping procedure
    a) Provide mapping files in .json format
    b) Provide heuristic like 3:1
2) Provide a "table of contetns" of molecules (so the system has knowledge of
    which molecules contain which atoms)
    a) Generate table of contents from structure file
        i) Gromacs is easy, already has molecule/residue info
        ii) Hoomd/lammps do not have molecule/residue info
    b) Supply a table of contents file
        i) This is json file mapping- 
            str(residuename-residueindex) : [global_atom_indices]
    c) If given heuristic information then start grouping together atoms 
        solely based on the fine-grained bonding patterns
        i) Utilize bondgraphs and finding neighbors
    d) All methods should create a cg hierarchy 
        cg-system (mb.Compound)-> molecules (mb.Compound)-> 
        cg-beads(str) -> atom_indices([int])
3) Do converting
    a) Iterate through the cg hierarchy, filling in the cg-system with particles
4) Create bonds
    a) If given mapping information, then bonding patterns will already be 
        specified
    b) If given heuristics, then assume linear bonding pattern?

cg_hierarchy
------------
cg_system (mb.Compound)
    cg_molecule (mb.Compound)
        cg_bead (str(bead_type-local_bead_index))
            global_atom_indices ([int])
        cg_bead (str)
        cg_bead (str)
    cg_molecule (mb.Compound)
    cg_molecule (mb.Compound)



"""

def forward_map(fine_grained, heuristic=None, mapping_files=[], 
        table_of_contents=None, dump_toc=True):
    """ Coarse grain an mb.Compound

    Parameters
    ---------
    fine_grained : mb.Compound
        original structure
    heuristic : str
        General mapping scheme i.e. '3:1'
    mapping_files : list of .json files
        json file specify residue names, mappings, and bonds
    dump_toc : bool
        Dump out table of contents to table_of_contents.json

     Molecular conversions should be a dictionary whose keys are cg beads
     and indices correspond to the global atom indices
     Table of contents matches molecules to global atom indices
    """
    # Heuristic denotes something like general 3:1 heavy atom mapping
    if heuristic:
        cg_system, cg_bonds, cg_hierarchy = _apply_heuristic(fine_grained, heuristic)

    # Mapping files denote a list of json files to create a large dictionary
    elif mapping_files:

        # Load in the mapping files into mapping_template
        mapping_template = {}
        for json_file in mapping_files:
            print("Loading json mapping <{}>".format(json_file))
            mapping_template.update(json.load(open(json_file,'r'), 
                object_pairs_hook=OrderedDict))

        # Create a table of contents to match molecules to global atom indices
        if not table_of_contents:
            print("Generating table of contents")
            table_of_contents = _extract_molecules(fine_grained, mapping_template)
            json.dump(table_of_contents, open('table_of_contents.json','w'), 
                    indent=4,separators=(',', ': '),ensure_ascii=False)
        else:
            print("Loading table of contents <{table_of_contents}>".format(**locals()))
            table_of_contents = json.load(open(table_of_contents,'r'))


        # Combine table of contents and mapping template
        # Create cg_system container for molecules and related beads
        # Identify bonds in the cg_system
        # Create the cg_hierarchy
        print("Generating hierarchy")
        cg_system, cg_bonds, cg_hierarchy = _create_cg_outline(table_of_contents, 
                mapping_template)
    else:
        print("More coarse-graining parameters needed, returning original compound")
        return fine_grained

    # Calculate centers of masses of groups of atoms to add particles
    # to cg_system
    print("Filling compound with particles")
    cg_system = _fill_cg_outline(fine_grained, cg_system, cg_hierarchy)

    # Apply bonding
    print("Applying bonds")
    cg_system = _apply_bonding(cg_system, cg_bonds)
    return cg_system

    
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
    particles = [p for p in coarse_grained.particles()]
    for bead_i, bead_j in cg_bonds:
        coarse_grained.add_bond([particles[bead_i], 
                                particles[bead_j]])
    return coarse_grained

def _fill_cg_outline(fine_grained, cg_system, cg_hierarchy):
    """ Add particles to cg system

    Parameters
    ----------
    fine_grained : mb.Compound
        original structure
    cg_system : mb.compound
        At this stage, no particles/beads added, but hierarachy established
    cg_hierarchy : OrderedDict()
        {molecule : {bead : [global atom indices]}}
    """
    
    for cg_molecule in cg_hierarchy.keys():
        for bead in cg_hierarchy[cg_molecule]:
            new_bead = mb.Compound(name=bead.split("-")[0], 
                    pos=_compute_center_of_mass(fine_grained, cg_hierarchy[cg_molecule][bead]))
            cg_molecule.add(new_bead)
        #cg_system.add(cg_molecule)

    return cg_system


    

def _create_cg_outline(table_of_contents, mapping_template):
    """ Read mapping.json files to outline molecular conversions

    table_of_contents : dict
        {resname-resindex : [global_atom_indices]}
    mapping_template : dict
        loaded-in json file

    Returns
    -------
    cg_system : mb.compound
        At this stage, no particles/beads added, but hierarachy established
    cg_bonds : list of tuples (2, n_beads)
        A list of 2-tuples that specifies bonded beads
    cg_hierarchy : OrderedDict()
        {mb.Compound molecule : {str bead : [global atom indices]}}

    """
    cg_bonds = []
    cg_system= mb.Compound()
    cg_hierarchy = OrderedDict()

    # Iterate molecule by molecule
    bead_counter = 0
    for molecule in table_of_contents.keys():
        molecule_name = molecule.split("-")[0]

        cg_molecule = mb.Compound(name=molecule_name)
        cg_system.add(cg_molecule)
        cg_hierarchy.update(OrderedDict({cg_molecule:OrderedDict()}))

        global_atom_indices = table_of_contents[molecule]


        # Update bonding for global BEAD indices
        # Take the mapping's local bonding
        # And shift appropriately by the global bead indices
        # Which is just the number of beads we've added 
        # before anything in this new molecule
        for bead_i, bead_j in mapping_template[molecule_name]['bond']:
            cg_bonds.append([bead_i + bead_counter, 
                        bead_j + bead_counter])

        # Take the mapping's local atom indices
        # And shift them appropriately by the global atom indices
        # Which is just adding the global atom index of the first
        # atom in the new molecule
        for bead in mapping_template[molecule_name]['map'].keys():
            updated_atom_indices = [local_index + global_atom_indices[0] for local_index in mapping_template[molecule_name]['map'][bead]]
            cg_hierarchy[cg_molecule].update({bead: updated_atom_indices})
            bead_counter+=1

    return cg_system, cg_bonds, cg_hierarchy


def _extract_molecules(fine_grained, mapping_template):
    """ From the structure file,
    be able to determine what molecule or residue each atom belongs to
    This should be easy with gro files,
    With hoomd/lammps files, there are no molecules/residues so this 
    has to be provided or somehow inferred

    Currently this is being implemented as a dict
    """
    has_res_information = True
    residues = [item for item in mapping_template.keys()]
    table_of_contents = OrderedDict()
    if has_res_information:
        traj = mdtraj.load('testing/bulk_DSPC_900K.gro')
        for residue in traj.topology.residues:
            indices = [at.index for at in residue.atoms]
            table_of_contents.update({residue.name + "-" + str(residue.index): indices})

    return table_of_contents


def _apply_heuristic(fine_grained, heuristic):
    """ Apply a heuristic to outline the molecular conversions
    i.e. 3:1 means this method will group 3 heavy atoms to 1 CG bead"""
    # Identify groups of atoms according to heuristic
    return None

def _compute_center_of_mass(fine_grained, atom_indices):
    """ Compute center of mass"""
    masses = [Mass[fine_grained.children[index].name[0:1]] for index in atom_indices]
    total_mass = sum(masses)
    com = np.ndarray(3)
    for i in range(3):
        com[i] = sum([fine_grained.children[index].pos[i]*mass/total_mass 
            for index,mass in zip(atom_indices,masses)])
    return com 



if __name__ == "__main__":
    fine_grained = mb.load('testing/two_propane.gro')
    mapping_files =['testing/propane.json']
    fine_grained.name = ""
    toc = 'testing/two_propane_toc.json'

    coarse_grained = forward_map(fine_grained, mapping_files=mapping_files, table_of_contents=toc)

    # Save
    coarse_grained.save('fwd_map.top',overwrite=True)
    coarse_grained.save('fwd_map.gro',overwrite=True, residues=['PR3'])
    coarse_grained.save('fwd_map.mol2',overwrite=True, residues=['PR3'])


    ## here comes the reversemapping
