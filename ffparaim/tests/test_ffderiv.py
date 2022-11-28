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
