"""
Unit and regression test for the ffparaim package. Testing mdtools.py.
"""

# Import package, test suite, and other packages as needed

import os
import pytest
import filecmp
import ffparaim.mdtools as mdt

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from parmed.structure import Structure
from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose


def test_separate_components(tmpdir):
    os.chdir(tmpdir)
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        mdt.separate_components(infile, ':MOL')
    assert filecmp.cmp(os.path.join(tmpdir, 'env.pdb'), files('ffparaim.data').joinpath('env.pdb')) == 0
    assert filecmp.cmp(os.path.join(tmpdir, 'lig.pdb'), files('ffparaim.data').joinpath('lig.pdb')) == 0


def test_separate_components_invalid():
    pytest.raises(TypeError, mdt.separate_components)
    pytest.raises(TypeError, mdt.separate_components, ligand_selection=':MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        pytest.raises(TypeError, mdt.separate_components, pdb_file=infile)
        pytest.raises(AttributeError, mdt.separate_components, infile, 0)


def test_define_molecule(tmpdir):
    os.chdir(tmpdir)
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        mdt.separate_components(infile, ':MOL')
    mol = mdt.define_molecule('c1ccc(cc1)O')
    assert isinstance(Molecule, mol)


def test_define_molecule_invalid():
    pytest.raises(TypeError, mdt.define_molecule)
    pytest.raises(OSError, mdt.define_molecule, 'CO', 'lig.pdb')
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as infile:
        pytest.raises(ValueError, mdt.define_molecule, 'CO', infile)


def test_define_forcefield():
    ff = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    assert isinstance(ForceField, ff)


def test_define_forcefield_invalid():
    pytest.raises(TypeError, mdt.define_forcefield)
    pytest.raises(OSError, mdt.define_forcefield, 'test')


def test_prepare_ligand(tmpdir):
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    assert isinstance(Structure, lig_structure)


def test_prepare_ligand_invalid(tmpdir):
    pytest.raises(TypeError, mdt.prepare_ligand)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    pytest.raises(TypeError, mdt.prepare_ligand, molecule)
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    pytest.raises(FileNotFoundError, mdt.prepare_ligand, molecule, forcefield)
