"""
Unit and regression test for the ffparaim package. Testing output.py.
"""

# Import package, test suite, and other packages as needed
import os
import shutil
import pytest
import filecmp
import ffparaim.qmtools as qmt
import ffparaim.mdtools as mdt

from ffparaim import ffparaim, utils
from importlib_resources import files, as_file


def test_write_qmmmm_pdb(tmpdir):
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    ligand_atom_list = mdt.get_atoms_idx(system_structure, 'resname MOL')
    with as_file(files('ffparaim.data').joinpath('output_recenter.pdb')) as pdb_file:
        shutil.copy(pdb_file, tmpdir)
    os.chdir(tmpdir)
    qmt.write_qmmm_pdb(ligand_atom_list)
    assert filecmp.cmp(os.path.join(tmpdir, 'output_qmmm.pdb'),
                       files('ffparaim.data').joinpath('output_qmmm.pdb')) == 0


def test_write_qmmm_pdb_invalid(tmpdir):
    os.chdir(tmpdir)
    pytest.raises(TypeError, qmt.write_qmmm_pdb)
    ligand_atom_list = list(range(0, 14))
    pytest.raises(RuntimeError, qmt.write_qmmm_pdb, ligand_atom_list, 'test.pdb')


def test_write_orca_input_qmmm(tmpdir):
    with as_file(files('ffparaim.data').joinpath('output_qmmm.pdb')) as pdb_file:
        shutil.copy(pdb_file, tmpdir)
    os.chdir(tmpdir)
    qmt.write_orca_input('qmmm',
                         atom=None,
                         ligand_selection='resname MOL',
                         pdb_file="output_qmmm.pdb",
                         method='B3LYP',
                         basis='def2-TZVP',
                         qm_charge=0,
                         qm_mult=1)
    assert filecmp.cmp(os.path.join(tmpdir, 'orca_qmmm.inp'),
                       files('ffparaim.data').joinpath('orca_qmmm.inp')) == 0


def test_write_orca_input_pol_corr(tmpdir):
    with as_file(files('ffparaim.data').joinpath('output_qmmm.pdb')) as pdb_file:
        shutil.copy(pdb_file, tmpdir)
    os.chdir(tmpdir)
    qmt.write_orca_input('pol_corr',
                         atom=None,
                         ligand_selection='resname MOL',
                         pdb_file="output_qmmm.pdb",
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
