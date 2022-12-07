"""
Unit and regression test for the ffparaim package. Testing ffderiv.py.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest
import filecmp
import ffparaim
import json

import numpy as np

from iodata import IOData
from grid.molgrid import MolGrid
from desnpart.vh import ProModel
from grid.basegrid import Grid
from openff.toolkit.topology import Molecule

from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose
from numpy.core._exceptions import UFuncTypeError


def test_ffderiv():
    ffd = ffparaim.ForceFieldDerivation()
    assert ffd.data is None
    assert ffd.grid is None
    assert ffd.rho is None
    assert ffd.pro_model is None
    assert ffd.localgrids is None
    assert ffd.results is None


def test_ffderiv_invalid():
    pytest.raises(TypeError, ffparaim.ForceFieldDerivation, 'test')


def test_load_data():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    assert isinstance(data, IOData) is True


def test_load_data_invalid():
    ffd = ffparaim.ForceFieldDerivation()
    pytest.raises(ValueError, ffd.load_data, 'data')


def test_set_molgrid():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    assert isinstance(ffd.grid, MolGrid) is True
    assert isinstance(ffd.data, dict) is True


def test_set_molgrid_invalid():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    pytest.raises(TypeError, ffd.set_molgrid)
    pytest.raises(AttributeError, ffd.set_molgrid, 'data')
    pytest.raises(TypeError, ffd.set_molgrid, data, nrad='150')
    pytest.raises(TypeError, ffd.set_molgrid, data, nang='194')
    pytest.raises(TypeError, ffd.set_molgrid, data, chunk_size='10000')
    pytest.raises(TypeError, ffd.set_molgrid, data, gradient='False')
    pytest.raises(TypeError, ffd.set_molgrid, data, orbitals='False')
    pytest.raises(TypeError, ffd.set_molgrid, data, store_atgrids='False')


def test_do_partitioning():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    results = ffd.do_partitioning(data)
    assert isinstance(results, dict)
    assert isinstance(ffd.pro_model, ProModel)
    assert isinstance(ffd.localgrids, Grid)
    assert_equal(ffd.data['density'][:3], np.array([123.69134467, 123.69134469, 123.69134469]))
    assert_equal(results['charges'], np.array([-1.23994956e-05]))
    assert_equal(results['radial_moments'][:3], np.array([5.99998444, 7.15516811, 13.94956547]))
    assert_equal(results['multipole_moments'][0][:3], np.array([-2.63543186e-06, 9.01856960e-06, 4.42787335e-06]))
    assert_equal(results['gtol'], 1e-8)
    assert_equal(results['maxiter'], 1000)
    assert_equal(results['density_cutoff'], 1e-10)
    assert_allclose(sum(ffd.pro_model.charges), 0, atol=1.e-4)


def test_do_partitioning_invalid():
    ffd = ffparaim.ForceFieldDerivation()
    pytest.raises(TypeError, ffd.do_partitioning)
    pytest.raises(TypeError, ffd.do_partitioning, 'data')
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    pytest.raises(AttributeError, ffd.do_partitioning, 'data')
    pytest.raises(ValueError, ffd.do_partitioning, data, 'hi')
    pytest.raises(TypeError, ffd.do_partitioning, data, gtol='1e-8')
    pytest.raises(TypeError, ffd.do_partitioning, data, maxiter='1000')
    pytest.raises(TypeError, ffd.do_partitioning, data, density_cutoff='1e-10')


def test_get_charges():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    ffd.do_partitioning(data)
    assert_equal(ffd.get_charges(), np.array([-1.23994956e-05]))


def test_get_epol():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    ffd.do_partitioning(data)
    assert_equal(ffd.get_epol(), 5.438830183804583)


def test_get_rcubed():
    ffd = ffparaim.ForceFieldDerivation()
    with as_file(files('ffparaim.data').joinpath('orca_uks.molden.input')) as infile:
        data = ffd.load_data(infile)
    ffd.set_molgrid(data)
    ffd.do_partitioning(data)
    assert_equal(ffd.get_charges(), np.array([5.17950988]))


def test_normalize_atomic_charges():
    mol = Molecule('C=O')
    norm_atcharges = ffparaim.ffderiv.normalize_atomic_charges(mol, 0, np.array([-0.3, 1, -0.6, -0.5]))
    assert_allclose(norm_atcharges, np.array([-0.2, 1.1, -0.5, -0.4]))


def test_normalize_atomic_charges_invalid():
    mol = Molecule('C=O')
    pytest.raises(TypeError, ffparaim.ffderiv.normalize_atomic_charges)
    pytest.raises(TypeError, ffparaim.ffderiv.normalize_atomic_charges, mol)
    pytest.raises(TypeError, ffparaim.ffderiv.normalize_atomic_charges, mol, 0)
    pytest.raises(AttributeError, ffparaim.ffderiv.normalize_atomic_charges, 'mol', 0, np.array([-0.3, 1, -0.6, -0.5]))
    pytest.raises(UFuncTypeError, ffparaim.ffderiv.normalize_atomic_charges, mol, '0', np.array([-0.3, 1, -0.6, -0.5]))
    pytest.raises(TypeError, ffparaim.ffderiv.normalize_atomic_charges, mol, 0, [-0.3, 1, -0.6, -0.5])


def test_symmetrize():
    mol = Molecule.from_smiles('C=O')
    symm = ffparaim.ffderiv.symmetrize(mol, [0.5627, -0.5151, -0.0240, -0.0236])
    assert isinstance(symm, list)
    assert_allclose(symm, [0.5627, -0.5151, -0.0238, -0.0238])


def test_symmetrize_invalid():
    mol = Molecule.from_smiles('C=O')
    pytest.raises(TypeError, ffparaim.ffderiv.symmetrize)
    pytest.raises(TypeError, ffparaim.ffderiv.symmetrize, mol)
    pytest.raises(AttributeError, ffparaim.ffderiv.symmetrize, 'mol', [0.5627, -0.5151, -0.0240, -0.0236])
    pytest.raises(TypeError, ffparaim.ffderiv.symmetrize, mol, '[0.5627, -0.5151, -0.0240, -0.0236]')
    pytest.raises(IndexError, ffparaim.ffderiv.symmetrize, mol, [0.5627, -0.5151, -0.0240])


def test_get_lj_params():
    mol = Molecule.from_smiles('C=O')
    with as_file(files('ffparaim.data').joinpath('rcubed_table.json')) as infile:
        rcubed_table = {int(k): v for k, v in json.load(infile).items()}
    sig, eps = ffparaim.ffderiv.get_lj_params(mol, np.array([1, 2, 3, 4]), rcubed_table)
    assert_allclose(sig, [0.27003539347805483, 0.2842451596185054, 0.3395921721548409, 0.3538392956260175])
    assert_allclose(eps, [0.06457184713591461, 0.16089255961492702, 0.4006858034641198, 0.5566613685093755])


def test_get_lj_params_invalid():
    mol = Molecule.from_smiles('C=O')
    pytest.raises(TypeError, ffparaim.ffderiv.get_lj_params)
    pytest.raises(TypeError, ffparaim.ffderiv.get_lj_params, mol)
    pytest.raises(TypeError, ffparaim.ffderiv.get_lj_params, mol, np.array([1, 2, 3, 4]))
    with as_file(files('ffparaim.data').joinpath('rcubed_table.json')) as infile:
        rcubed_table = json.load(infile)
    pytest.raises(KeyError, ffparaim.ffderiv.get_lj_params, mol, np.array([1, 2, 3, 4]), rcubed_table)
    rcubed_table = {int(k): v for k, v in rcubed_table.items()}
    pytest.raises(AttributeError, ffparaim.ffderiv.get_lj_params, 'mol', np.array([1, 2, 3, 4]), rcubed_table)
    pytest.raises(TypeError, ffparaim.ffderiv.get_lj_params, mol, 'np.array([1, 2, 3, 4])', rcubed_table)
    pytest.raises(TypeError, ffparaim.ffderiv.get_lj_params, mol, np.array([1, 2, 3, 4]), 'rcubed_table')
    pytest.raises(KeyError, ffparaim.ffderiv.get_lj_params, mol, rcubed_table, np.array([1, 2, 3, 4]))

