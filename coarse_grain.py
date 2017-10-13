import pdb
import numpy as np
import mbuild as mb
import mdtraj
import json
from mapping_moieties.ch2_aa import CH2_aa
from mapping_moieties.ch3_aa import CH3_aa

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
        ?Parmed?
        *Get molecules based on bonding networks
        *Bead mbuild classes with atom children
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
        cg_bead (str(bead_type + "-" + local_bead_index))
            global_atom_indices ([int])
        cg_bead (str)
        cg_bead (str)
    cg_molecule (mb.Compound)
    cg_molecule (mb.Compound)

"""
""" Reverse-mapping an mbuild compound

1) Specify a mapping procedure
    a) .json format
    b) Should relate a cg_local_index to atomtypes
2) Obtain "table of contents" to that relates molecule to respective beads
    a) str(residuename-residueindex): [global_bead_indices]
3) Do converting
    a) For each bead, replace it with atom-particles centered around that point
        i) This could be accomplished using mb.Spherepattern
4) Create bonds
"""

def reverse_map(coarse_grained, heuristic=None, 
        table_of_contents=None, dump_toc=True):
    """ Reverse map an mb.Compound

    Parameters
    ---------
    coarse_grained : mb.Compound
        original structure
    table_of_contents : dictionary
        map molecules to global atom indices

    """
    # Get molecular information through bonding 
    table_of_contents = _detect_molecules(coarse_grained)
    # For now, just load it in 
    #table_of_contents = json.load(open(table_of_contents,'r'),
        #object_pairs_hook=OrderedDict)

    # Establish the correct mappings
    mapping_moieties = {'CH3':CH3_aa, 
                        'CH2':CH2_aa}

    aa_system = mb.Compound()
    # CG to AA relates the CG bead to its AA bead
    cg_to_aa = OrderedDict()

    # For each bead, replace it with the appropraite mb compound
    particles = [p for p in coarse_grained.particles()]
    #for molecule in table_of_contents.keys():
    for molecule in table_of_contents:
        new_molecule =  mb.Compound()
        #for bead_index in table_of_contents[molecule]:
        for bead in molecule:
            #new_atom = mapping_moieties[particles[bead_index].name]()
            new_atom = mapping_moieties[bead.name]()
            cg_to_aa[bead] = new_atom
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)

    # Go back and include bonds
    for p_i, p_j in coarse_grained.bonds():
        mb.force_overlap(cg_to_aa[p_i],
                from_positions=cg_to_aa[p_i].available_ports()[0],
                to_positions=cg_to_aa[p_j].available_ports()[0])

    # Translate atoms after they've been bonded
    for cg,aa in cg_to_aa.items():
        aa.translate_to(cg.pos)
    return aa_system
    

def _detect_molecules(coarse_grained):
    """ Get each molecule using bond graph

    Returns
    ------
    molecules: list of lists
        Each list contains the connected components
        However it isn't ordered by connectivity
    """
    molecules = coarse_grained.bond_graph.connected_components()

    return molecules



















    #if mapping_files:

    #    # Load in the mapping files into mapping_template
    #    mapping_template = OrderedDict()
    #    for json_file in mapping_files:
    #        print("Loading json mapping <{}>".format(json_file))
    #        mapping_template.update(json.load(open(json_file,'r'), 
    #            object_pairs_hook=OrderedDict))

    #    # Create a table of contents to match molecules to global atom indices
    #    if not table_of_contents:
    #        print("Generating table of contents")
    #        table_of_contents = _extract_molecules(fine_grained, mapping_template)
    #        json.dump(table_of_contents, open('table_of_contents.json','w'), 
    #                indent=4,separators=(',', ': '),ensure_ascii=False)
    #    else:
    #        print("Loading table of contents <{table_of_contents}>".format(**locals()))
    #        table_of_contents = json.load(open(table_of_contents,'r'),
    #                object_pairs_hook=OrderedDict)

    #    print("Generating hierarchy")
    #    aa_system, aa_bonds, aa_hierarchy = _create_aa_outline(coarse_grained, 
    #            table_of_contents, mapping_template)

    #else: 
    #    print("Not enough reverse-mapping parameters specified, returning original compound")
    #    return coarse_grained

    ## AA system and AA hierarchy have atoms in it
    ## But need to place them around the cg bead center
    #print("Inflating system with atoms")
    #aa_system = _restore_atoms(coarse_grained, aa_system, aa_hierarchy)

    ## Apply bonding
    #print("Applying bonds")
    #aa_system = _apply_bonding(aa_system, aa_bonds)
    #return aa_system

def _extract_molecules(fine_grained, mapping_template):
    """ Return a dictionary mapping resname-resindex to global indices
    
    Parameters
    ---------
    mapping_template : OrderedDict()
        mapping template from json
    From the structure file,
    be able to determine what molecule or residue each atom belongs to
    This should be easy with gro files,
    With hoomd/lammps files, there are no molecules/residues so this 
    has to be provided or somehow inferred

    Currently this is being implemented as a dict
    Notes
    -----
    This can generate a table of contents for both cg and aa compounds,
    depending on what is the structure provided
    """
    has_res_information = True
    residues = [item for item in mapping_template.keys()]
    table_of_contents = OrderedDict()
    if has_res_information:
        traj = mdtraj.load('testing/two_propane.gro')
        for residue in traj.topology.residues:
            indices = [at.index for at in residue.atoms]
            table_of_contents.update({residue.name + "-" + str(residue.index): indices})

    return table_of_contents


def _create_aa_outline(coarse_grained, table_of_contents, mapping_template):
    """ Read mapping.json files to outline conversions

    table_of_contents : dict
        {resname-resindex : [global_bead_indices]}
    mapping_template : dict
        loaded-in json file

    Returns
    -------
    aa_system : mb.compound
        Hierarchy established, positions still unset
        mb.Compound molecule > mb.Compound atom

    aa_bonds : list of tuples (2, n_bonds)
        A list of 2-tuples that specifies bonded atoms
    aa_hierarchy : OrderedDict()
        {mb.Compound molecule : {mb.Compound bead : [mb.Compound atom ]}}

    """
    aa_bonds = []
    aa_system = mb.Compound()
    aa_hierarchy = OrderedDict()
    atom_counter = 0
    bead_counter = 0
    particles = [p for p in coarse_grained.particles()]
    for molecule in table_of_contents.keys():
        molecule_name = molecule.split("-")[0]

        aa_molecule = mb.Compound(name=molecule_name)
        aa_system.add(aa_molecule)
        #aa_hierarchy.update(OrderedDict({aa_molecule:OrderedDict()}))
        aa_hierarchy.update([ (aa_molecule, OrderedDict()) ])

        global_bead_indices = table_of_contents[molecule]
        local_atom_list = [""]*int(mapping_template[molecule_name]['n_atoms'])

        # Update bonding for global atom indices
        # Take the mapping's local aa bonding
        # And shift appropriately by the global atom indices
        # Which is just the number of atoms we've added
        # before anything in thie new molecule
        for atom_i, atom_j in mapping_template[molecule_name]['aa_bond']:
            aa_bonds.append([atom_i + atom_counter, atom_j + atom_counter])


        # For each bead in the molecule's mapping
        # Create a local atom list corresponding to all the atoms
        # in that molecule
        # Be careful about bookkeeping because the the atom order
        # doesn't follow the bead order perfectly
        # Add the local atom list to the aa system
        # aa hierarchy stil keeps molecules -> beads -> atoms organized
        # aa system is molecules -> atoms
        for bead in mapping_template[molecule_name]['map'].keys():
            updated_bead_index = int(bead.split("-")[1]) + int(global_bead_indices[0])
            cg_bead = mb.Compound(name=bead, 
                    pos=particles[updated_bead_index].pos)
            #aa_hierarchy[aa_molecule].update(OrderedDict({cg_bead:[]}))
            aa_hierarchy[aa_molecule].update([ (cg_bead, []) ])
            for atom in mapping_template[molecule_name]['map'][bead]:
                atom_name, local_index = atom.split("-")
                local_atom = mb.Compound(name=atom_name, 
                        pos=particles[updated_bead_index].pos)
                local_atom_list[int(local_index)] = local_atom

                aa_hierarchy[aa_molecule][cg_bead].append(local_atom)
                atom_counter+=1
            bead_counter+=1


        for atom in local_atom_list:
            aa_molecule.add(atom)
    return aa_system, aa_bonds, aa_hierarchy


def _restore_atoms(coarse_grained, aa_system, aa_hierarchy):
    """ Optimize aa coordinates

    Parameters
    ----------
    coarse_grained : mb.Compound
        CG structure
    aa_system : mb.Compound
        AA structure with unoptimized coordinates
    aa_hierarchy : OrderedDict()
        {mb.Compound molecule : {mb.Compound bead : [mb.Compound atom ]}}


    Returns
    -------
    aa_system : mb.Compound
        AA structure with optimized coordinates

    Notes
    -----
    Generates a sphere of points around each bead

        """
    for molecule in aa_hierarchy.keys():
        for bead in aa_hierarchy[molecule]:
            spherical_pattern = mb.SpherePattern(
                    n=len(aa_hierarchy[molecule][bead])+1,scale=0.1)
            sphere = spherical_pattern.apply(bead) 
            for index, atom in enumerate(aa_hierarchy[molecule][bead]):
                atom.translate(sphere[index].pos)



    return aa_system




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


def _create_cg_outline(table_of_contents, mapping_template):
    """ Read mapping.json files to outline conversions

    table_of_contents : dict
        {resname-resindex : [global_atom_indices]}
    mapping_template : dict
        loaded-in json file

    Returns
    -------
    cg_system : mb.compound
        At this stage, no particles/beads added, but hierarachy established
    cg_bonds : list of tuples (2, n_bonds)
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
        #cg_hierarchy.update(OrderedDict({cg_molecule:OrderedDict()}))
        cg_hierarchy.update([ (cg_molecule, OrderedDict()) ])

        global_atom_indices = table_of_contents[molecule]


        # Update bonding for global BEAD indices
        # Take the mapping's local bonding
        # And shift appropriately by the global bead indices
        # Which is just the number of beads we've added 
        # before anything in this new molecule
        for bead_i, bead_j in mapping_template[molecule_name]['cg_bond']:
            cg_bonds.append([bead_i + bead_counter, 
                        bead_j + bead_counter])

        # Take the mapping's local atom indices
        # And shift them appropriately by the global atom indices
        # Which is just adding the global atom index of the first
        # atom in the new molecule
        for bead in mapping_template[molecule_name]['map'].keys():
            updated_atom_indices = [int(local_index.split('-')[1]) + global_atom_indices[0] for local_index in mapping_template[molecule_name]['map'][bead]]
            #cg_hierarchy[cg_molecule].update({bead: updated_atom_indices})
            cg_hierarchy[cg_molecule].update([ (bead, updated_atom_indices) ])
            bead_counter+=1

    return cg_system, cg_bonds, cg_hierarchy

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

    return cg_system

if __name__ == "__main__":
    #fine_grained = mb.load('testing/two_propane.gro')
    mapping_files =['testing/propane.json']
    #fine_grained.name = ""
    residues=['PR3']
    #toc = 'testing/two_propane_aa_toc.json'
    ##toc = None

    #coarse_grained = forward_map(fine_grained, mapping_files=mapping_files, table_of_contents=toc)

    ## Save
    #coarse_grained.save('fwd_map.top',overwrite=True)
    #coarse_grained.save('fwd_map.gro',overwrite=True, residues=residues)
    #coarse_grained.save('fwd_map.mol2',overwrite=True, residues=residues)


    ## here comes the reversemapping
    #toc = None
    #toc = 'cg_toc.json'
    toc = 'testing/two_propane_cg_toc.json'
    coarse_grained = mb.load('sample_cg_two_propane.mol2')
    recovered = reverse_map(coarse_grained, table_of_contents=toc)
    recovered.save('rev_map.gro',overwrite=True, residues=residues)
    recovered.save('rev_map.mol2',overwrite=True, residues=residues)



