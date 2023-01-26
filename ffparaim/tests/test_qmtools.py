"""
Unit and regression test for the ffparaim package. Testing output.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
import filecmp
import ffparaim.qmtools as qmt

from ffparaim import utils
from importlib_resources import files, as_file


def test_set_qm_atoms():
    with as_file(files('ffparaim.data').joinpath('output.pdb')) as pdb_file:
        lig_atoms_idx = qmt.set_qm_atoms(':MOL', pdb_file)
        assert lig_atoms_idx[:4] == [0, 1, 2, 3]


def test_set_qm_atoms_invalid(tmpdir):
    os.chdir(tmpdir)
    pytest.raises(TypeError, qmt.set_qm_atoms)
    pytest.raises(FileNotFoundError, qmt.set_qm_atoms, ':MOL', 'test.pdb')


def test_write_qmmmm_pdb(tmpdir):
    os.chdir(tmpdir)
    with as_file(files('ffparaim.data').joinpath('output.pdb')) as pdb_file:
        lig_atoms_idx = qmt.set_qm_atoms(':MOL', pdb_file)
    qmt.write_qmmm_pdb(lig_atoms_idx)
    assert filecmp.cmp(os.path.join(tmpdir, 'output_qmmm.pdb'),
                       files('ffparaim.data').joinpath('output_qmmm.pdb')) == 0


def test_write_qmmm_pdb_invalid(tmpdir):
    os.chdir(tmpdir)
    pytest.raises(TypeError, qmt.write_qmmm_pdb)
    lig_atoms_idx = list(range(0, 14))
    pytest.raises(FileNotFoundError, qmt.set_qm_atoms, lig_atoms_idx, 'test.pdb')


def test_write_orca_input_qmmm(tmpdir):
    os.chdir(tmpdir)
    qmt.write_orca_input('qmmm',
                         atom=None,
                         ligand_selection=':MOL',
                         pdb_file="output_recenter.pdb",
                         method='B3LYP',
                         basis='def2-TZVP',
                         qm_charge=0,
                         qm_mult=1)
    assert filecmp.cmp(os.path.join(tmpdir, 'orca_qmmm.inp'),
                       files('ffparaim.data').joinpath('orca_qmmm.inp')) == 0


def test_write_orca_input_pol_corr(tmpdir):
    os.chdir(tmpdir)
    qmt.write_orca_input('pol_corr',
                         atom=None,
                         ligand_selection=':MOL',
                         pdb_file="output_recenter.pdb",
                         method='B3LYP',
                         basis='def2-TZVP',
                         qm_charge=0,
                         qm_mult=1)
    assert filecmp.cmp(os.path.join(tmpdir, 'orca_pol_corr.inp'),
                       files('ffparaim.data').joinpath('orca_pol_corr.inp')) == 0


def test_write_orca_input_uks(tmpdir):
    os.chdir(tmpdir)
    qmt.write_orca_input('uks',
                         atom=utils.elements[6],
                         qm_mult=utils.mult_table[6])
    assert filecmp.cmp(os.path.join(tmpdir, 'orca_uks.inp'),
                       files('ffparaim.data').joinpath('orca_uks.inp')) == 0


def test_write_orca_input_invalid():
    pytest.raises(TypeError, qmt.write_orca_input)
    pytest.raises(ValueError, qmt.write_orca_input, 'md')
