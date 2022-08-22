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
    print('Setting restraints ...')
    if restraint_dict is not None:
        thermodynamic_state = states.ThermodynamicState(system, 298.15 * unit.kelvin)
        sampler_state = states.SamplerState(positions)
        if ligand_atom_list is not None:
            topography = Topography(top, ligand_atoms=ligand_atom_list)
        if 'harmonic' in restraint_dict.keys():
            restraint = Harmonic(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                 restrained_receptor_atoms=topography.select(restraint_dict['harmonic']['receptor']),
                                 restrained_ligand_atoms=topography.select(restraint_dict['harmonic']['ligand']))
        if 'flatbottom' in restraint_dict.keys():
            restraint = FlatBottom(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                   well_radius=10.0 / unit.angstrom,
                                   restrained_receptor_atoms=topography.select(restraint_dict['flatbottom']['receptor']),
                                   restrained_ligand_atoms=topography.select(restraint_dict['flatbottom']['ligand']))

        if 'boresch' in restraint_dict.keys():
            restraint = Boresch(restrained_receptor_atoms=topography.select(restraint_dict['boresch']['receptor']),
                                restrained_ligand_atoms=topography.select(restraint_dict['boresch']['ligand']),
                                K_r=4184.0 * unit.kilojoule_per_mole / unit.nanometer ** 2,
                                K_thetaA=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_thetaB=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiA=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiB=41.84 * unit.kilojoule_per_mole / unit.radians ** 2,
                                K_phiC=41.84 * unit.kilojoule_per_mole / unit.radians ** 2)
        try:
            restraint.restrain_state(thermodynamic_state)
        except RestraintParameterError:
            print('Choosing restraint parameters automatically.')
            restraint.determine_missing_parameters(thermodynamic_state, sampler_state, topography)
            restraint.restrain_state(thermodynamic_state)

    return
