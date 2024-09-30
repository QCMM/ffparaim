"""
Unit and regression test for the ffparaim package. Testing restraints.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
import ffparaim
import ffparaim.restraints as restr

from openmm import System
from openmmtools import forces
from collections import OrderedDict
from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose


def test_set_restraints_harmonic(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    restraint_dict = {'harmonic': {'receptor': 'resid 99', 'ligand': 'resname MOL'}}
    system = restr.set_restraints(system_structure.topology,
                                  system,
                                  restraint_dict)
    assert isinstance(system, System) is True
    restraint_force = forces.find_forces(system, r'\bHarmonicRestraint', only_one=False)
    assert isinstance(restraint_force, OrderedDict)
    assert_allclose(list(restraint_force.values())[0].getBondParameters(0)[1][0], 83.68)
    assert_equal(list(restraint_force.values())[0].getForceGroup(), 0)


def test_set_restraints_flatbottom(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    restraint_dict = {'flatbottom': {'receptor': 'resid 99', 'ligand': 'resname MOL'}}
    system = restr.set_restraints(system_structure.topology,
                                  system,
                                  restraint_dict)
    assert isinstance(system, System) is True
    restraint_force = forces.find_forces(system, r'\bFlatBottom', only_one=False)
    assert isinstance(restraint_force, OrderedDict)
    assert_allclose(list(restraint_force.values())[0].getBondParameters(0)[1][1], 100.0)
    assert_equal(list(restraint_force.values())[0].getForceGroup(), 0)


def test_set_restraints_invalid(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('complex.pdb')) as infile:
        molecule, lig_strucure, system_structure, system = ffp.prepare('CC(=O)OCc1ccccc1', str(infile))
    pytest.raises(TypeError, restr.set_restraints)
    pytest.raises(TypeError, restr.set_restraints, system_structure.topology)
