#!/usr/bin/python3

import mdtraj
from openmmtools import forces
from openmm import unit


def set_restraints(topology,
                   system,
                   restraint_dict=None):
    """Generate different restraint types for receptor-ligand complexes."""

    # Print message for setting restraints.
    print('Setting restraints ...')
    top = mdtraj.Topology.from_openmm(topology)
    # Check for restraints added to the system.
    restraint_force = forces.find_forces(system, r'\bFlatBottom|\bHarmonicRestraint', only_one=False)
    # If restraint dict is defined and there is a restraint definition in system.
    if restraint_dict is not None and len(restraint_force) == 0:
        # If harmonic restraints are defined in the dictionary.
        if 'harmonic' in restraint_dict.keys():
            # Create a Harmonic restraint object using selected receptor and ligand atoms.
            restraint = forces.HarmonicRestraintForce(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                                      restrained_atom_indices1=top.select(restraint_dict['harmonic']['receptor']),
                                                      restrained_atom_indices2=top.select(restraint_dict['harmonic']['ligand']))
        # If flat-bottom restraints are defined in the dictionary.
        if 'flatbottom' in restraint_dict.keys():
            # Create a FlatBottom restraint object using selected receptor and ligand atoms.
            restraint = forces.FlatBottomRestraintForce(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                                        well_radius=10.0 / unit.angstrom,
                                                        restrained_atom_indices1=top.select(restraint_dict['flatbottom']['receptor']),
                                                        restrained_atom_indices2=top.select(restraint_dict['flatbottom']['ligand']))
        system.addForce(restraint)
        return system
    else:
        return system
