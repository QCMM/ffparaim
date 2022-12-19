"""
Unit and regression test for the ffparaim package. Testing mdtools.py.
"""

# Import package, test suite, and other packages as needed

import os
import pytest
import filecmp
import ffparaim.mdtools as mdt

from openff.toolkit.topology import Molecule
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
