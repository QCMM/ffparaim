"""
Unit and regression test for the ffparaim package. Testing ffderiv.py.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest
import filecmp
import ffparaim

from iodata import IOData
from grid.molgrid import MolGrid

from importlib_resources import files, as_file
from numpy.testing import assert_equal


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
    pytest.raises(TypeError, ffd.set_molgrid, 'data')
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
