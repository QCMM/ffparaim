"""
Unit and regression test for the ffparaim package. Testing ffparaim.py.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest
import filecmp
import ffparaim


from parmed import Structure
from openmm import System
from openff.toolkit.topology import Molecule

from importlib_resources import files, as_file
from numpy.testing import assert_equal


def test_ffparaim_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ffparaim" in sys.modules


def test_ffparaim():
    ffp = ffparaim.FFparAIM()
    assert ffp.qm_charge == 0
    assert ffp.ligand_selection == ':1'
    assert ffp.ligand_atom_list is None
    assert ffp.n_updates == 3
    assert ffp.sampling_time == 25
    assert ffp.total_qm_calculations == 100
    assert ffp.method == 'B3LYP'
    assert ffp.basis == 'def2-TZVP'
    assert any(ffp.data) is False


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_ffparaim_prepare(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        molecule, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(pdb))
    assert ffp.smiles == 'c1ccc(cc1)O'
    assert 'solvent.pdb' in ffp.pdb_file
    assert ffp.forcefield == 'openff_unconstrained-2.0.0.offxml'
    assert filecmp.cmp(os.path.join(tmpdir, 'env.pdb'), files('ffparaim.data').joinpath('env.pdb')) == 0
    assert filecmp.cmp(os.path.join(tmpdir, 'lig.pdb'), files('ffparaim.data').joinpath('lig.pdb')) == 0
    assert isinstance(molecule, Molecule) is True
    assert isinstance(system_structure, Structure) is True
    assert isinstance(system, System) is True


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_ffparaim_prepare_invalid(tmpdir):
    ffp = ffparaim.FFparAIM()
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        pytest.raises(TypeError, ffp.prepare)
        pytest.raises(TypeError, ffp.prepare, 1, str(pdb))
        pytest.raises(ValueError, ffp.prepare, 'CO', str(pdb))


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_ffparaim_run(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(n_updates=1, sampling_time=0.001, total_qm_calculations=1)
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as pdb:
        molecule, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(pdb))
    ffp.run(molecule, system_structure, system, off=True)
    assert_equal(ffp.ligand_atom_list, range(0, 13))
