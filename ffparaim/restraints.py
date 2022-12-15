#!/usr/bin/python3

from yank.restraints import Harmonic, FlatBottom, Boresch, RestraintParameterError
from openmmtools import states
from yank.yank import Topography
from simtk import unit


def set_restraints(top,
                   system,
                   positions,
                   restraint_dict=None,
                   ligand_atom_list=None):
    """Generate different restraint types for receptor-ligand complexes."""

    # Print message for setting restraints.
    print('Setting restraints ...')
    # If restraint dict is defined.
    if restraint_dict is not None:
        # Create a state with the system at a reference temperature.
        thermodynamic_state = states.ThermodynamicState(system, 298.15 * unit.kelvin)
        # Extract coordinates from this state.
        sampler_state = states.SamplerState(positions)
        # If ligand atom index list is defined.
        if ligand_atom_list is not None:
            # Create a Topography class from system's topology.
            topography = Topography(top, ligand_atoms=ligand_atom_list)
        # If harmonic restraints are defined in the dictionary.
        if 'harmonic' in restraint_dict.keys():
            # Create a Harmonic restraint object using selected receptor and ligand atoms.
            restraint = Harmonic(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                 restrained_receptor_atoms=topography.select(restraint_dict['harmonic']['receptor']),
                                 restrained_ligand_atoms=topography.select(restraint_dict['harmonic']['ligand']))
        # If flat-bottom restraints are defined in the dictionary.
        if 'flatbottom' in restraint_dict.keys():
            # Create a FlatBottom restraint object using selected receptor and ligand atoms.
            restraint = FlatBottom(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                   well_radius=10.0 / unit.angstrom,
                                   restrained_receptor_atoms=topography.select(restraint_dict['flatbottom']['receptor']),
                                   restrained_ligand_atoms=topography.select(restraint_dict['flatbottom']['ligand']))
        # If Boresch-Karplis restraints are defined in the dictionary.
        if 'boresch' in restraint_dict.keys():
            # Create a Boresch restraint object using selected receptor and ligand atoms.
            restraint = Boresch(restrained_receptor_atoms=topography.select(restraint_dict['boresch']['receptor']),
                                restrained_ligand_atoms=topography.select(restraint_dict['boresch']['ligand']),
                                K_r=4184.0 * unit.kilojoule_per_mole / unit.nanometer ** 2,
                                K_thetaA=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_thetaB=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiA=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiB=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiC=41.84 * unit.kilojoule_per_mole / unit.radians ** 2)
        # If they aren't missing parameters.
        try:
            # Create the restrainted state.
            restraint.restrain_state(thermodynamic_state)
        # If missing parameters are present.
        except RestraintParameterError:
            # Print message for choosing parameters.
            print('Choosing equilibrium distance, angles and torsions automatically.')
            # Generate missing parameters for the restraint.
            restraint.determine_missing_parameters(thermodynamic_state, sampler_state, topography)
            # Create the restrained state.
            restraint.restrain_state(thermodynamic_state)

    return thermodynamic_state.get_system()
