"""
Unit and regression test for the ffparaim package. Testing restraints.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
import ffparaim
import ffparaim.mdtools as mdt
import ffparaim.restraints as restr

from openmm import System
from openmmtools import forces
from collections import OrderedDict
from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose


def test_set_restraints_boresch(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    positions = system_structure.positions
    ligand_atom_list = mdt.get_atoms_idx(system_structure, ':MOL')
    restraint_dict = {'boresch': {'receptor': [1627, 1640, 1642], 'ligand': [0, 11, 9]}}
    system = restr.set_restraints(system_structure.topology,
                                  system,
                                  positions,
                                  restraint_dict,
                                  ligand_atom_list)
    assert isinstance(system, System) is True
    restraint_force = forces.find_forces(system, r'\bCustom', only_one=False)
    assert isinstance(restraint_force, OrderedDict)
    assert list(restraint_force.values())[0].getBondParameters(0)[0] == (1627, 1640, 1642, 0, 11, 9)
    assert_equal(list(restraint_force.values())[0].getForceGroup(), 3)


def test_set_restraints_harmonic(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    positions = system_structure.positions
    ligand_atom_list = mdt.get_atoms_idx(system_structure, ':MOL')
    restraint_dict = {'harmonic': {'receptor': 'resid 99', 'ligand': 'resname MOL'}}
    system = restr.set_restraints(system_structure.topology,
                                  system,
                                  positions,
                                  restraint_dict,
                                  ligand_atom_list)
    assert isinstance(system, System) is True
    restraint_force = forces.find_forces(system, r'\bHarmonicRestraint', only_one=False)
    assert isinstance(restraint_force, OrderedDict)
    assert_allclose(list(restraint_force.values())[0].getBondParameters(0)[1][0], 83.68)
    assert_equal(list(restraint_force.values())[0].getForceGroup() == 3)


def test_set_restraints_flatbottom(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    positions = system_structure.positions
    ligand_atom_list = mdt.get_atoms_idx(system_structure, ':MOL')
    restraint_dict = {'flatbottom': {'receptor': 'resid 99', 'ligand': 'resname MOL'}}
    system = restr.set_restraints(system_structure.topology,
                                  system,
                                  positions,
                                  restraint_dict,
                                  ligand_atom_list)
    assert isinstance(system, System) is True
    restraint_force = forces.find_forces(system, r'\bFlatBottom', only_one=False)
    assert isinstance(restraint_force, OrderedDict)
    assert_allclose(list(restraint_force.values())[0].getBondParameters(0)[1][1], 100.0)
    assert_equal(list(restraint_force.values())[0].getForceGroup() == 3)


def test_set_restraints_invalid(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    pytest.raises(TypeError, restr.set_restraints)
    pytest.raises(TypeError, restr.set_restraints, system_structure.topology)
    pytest.raises(TypeError, restr.set_restraints, system_structure.topology, system)
    pytest.raises(TypeError,
                  restr.set_restraints,
                  system_structure.topology,
                  system)
