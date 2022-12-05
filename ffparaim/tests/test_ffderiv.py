"""
Unit and regression test for the ffparaim package. Testing ffderiv.py.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest
import filecmp
import ffparaim
import numpy as np

from iodata import IOData
from grid.molgrid import MolGrid
from desnpart.vh import ProModel
from grid.basegrid import Grid

from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose


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
