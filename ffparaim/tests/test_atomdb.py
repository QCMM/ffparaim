"""
Unit and regression test for the ffparaim package. Testing atomdb.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
import ffparaim


from openff.toolkit.topology import Molecule

from importlib_resources import files, as_file


def test_atomdb(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(pdb))
    atomdb = ffparaim.atomdb.AtomDB(molecule, ffp.method, ffp.basis)
    assert isinstance(atomdb.molecule, Molecule) is True
    assert atomdb.method == 'B3LYP'
    assert atomdb.basis == 'def2-TZVP'


def test_atomdb_invalid(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    pytest.raises(TypeError, ffparaim.atomdb.AtomDB)
    pytest.raises(ValueError, ffparaim.atomdb.AtomDB, 'CO', ffp.method, ffp.basis)
