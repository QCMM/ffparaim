"""
Unit and regression test for the ffparaim package. Testing output.py.
"""

# Import package, test suite, and other packages as needed
import os
import pytest
import filecmp
import pandas as pd
import ffparaim.mdtools as mdt

from ffparaim import output
from importlib_resources import files, as_file


def test_to_pickle(tmpdir):
    os.chdir(tmpdir)
    data = dict()
    output.to_pickle(data)
    assert os.path.exists(os.path.join(tmpdir, 'ffparaim.pickle')) is True


def test_to_pickle_invalid():
    pytest.raises(TypeError, output.to_pickle)


def test_to_csv(tmpdir):
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        lig_structure = mdt.prepare_ligand(molecule, forcefield, lig_pdb_file)
    with as_file(files('ffparaim.data').joinpath('ffparaim.dat')) as data:
        df = pd.read_csv(data)
    charges = df.atomic_charges.to_list()
    sig = df.sigma_nanometer.to_list()
    eps = df.epsilon_kjmol.to_list()
    output.to_csv(lig_structure, charges, sig, eps)
    assert filecmp.cmp(os.path.join(tmpdir, 'ffparaim.dat'),
                       files('ffparaim.data').joinpath('ffparaim.dat')) == 0


def test_to_csv_invalid():
    pytest.raises(TypeError, output.to_csv)


def test_to_dat(tmpdir):
    os.chdir(tmpdir)
    output.to_dat(0.5, 0.25, 'test.out')
    assert os.path.exists(os.path.join(tmpdir, 'test.out')) is True
    with open(os.path.join(tmpdir, 'test.out'), 'r') as f:
        line = f.readline()
    assert line.split()[0] == '0.500000'
    assert line.split()[1] == '0.250000'


def test_to_dat_invalid():
    pytest.raises(TypeError, output.to_dat)
