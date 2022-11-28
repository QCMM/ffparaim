"""
Unit and regression test for the ffparaim package. Testing atomdb.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
# import filecmp
import ffparaim


from openff.toolkit.topology import Molecule

from importlib_resources import files, as_file
# from numpy.testing import assert_equal, assert_allclose


def test_atomdb(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        molecule, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(pdb))
    atomdb = ffparaim.AtomDB(molecule, ffp.method, ffp.basis)
    assert isinstance(atomdb.molecule, Molecule) is True
    assert atomdb.method == 'B3LYP'
    assert atomdb.basis == 'def2-TZVP'


def test_atomdb_invalid(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    pytest.raises(TypeError, ffparaim.AtomDB)
    pytest.raises(ValueError, ffparaim.AtomDB, 'CO', ffp.method, ffp.basis)


'''def test_create_table(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        molecule, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(pdb))
    atomdb = ffparaim.AtomDB(molecule, ffp.method, ffp.basis)
    rcubed_table = atomdb.create_table()
    rcubed_table_ref = {"8": 3.2555507175615275, "1": 1.1713382554315122, "6": 5.179509879206654}
    assert filecmp.cmp(os.path.join(tmpdir, 'rcubed_table.json'), files('ffparaim.data').joinpath('rcubed.json')) == 0
    assert rcubed_table == rcubed_table_ref
'''
