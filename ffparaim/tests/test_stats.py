"""
Unit and regression test for the ffparaim package. Testing stats.py.
"""

# Import package, test suite, and other packages as needed
import pytest

from ffparaim import stats
from iodata import IOData
from numpy.testing import assert_allclose


def test_nb_stats():
    data = dict()
    data[0] = list()
    iod = IOData(atffparams={'charges': list(range(0, 14)),
                             'rcubed': list(range(14, 24))})
    data[0].append(iod)
    charges_mean, charges_std = stats.nb_stats(data[0], 'charges')
    rcubed_mean, rcubed_std = stats.nb_stats(data[0], 'rcubed')
    assert_allclose(charges_mean[-1], 13.0)
    assert_allclose(charges_std[0], 0.0)
    assert_allclose(rcubed_mean[-1], 23.0)
    assert_allclose(rcubed_std[0], 0.0)


def test_nb_stats_invalid():
    pytest.raises(TypeError, stats.nb_stats)
    data = dict()
    data[0] = list()
    iod = IOData(atffparams={'charges': list(range(0, 14)),
                             'rcubed': list(range(14, 24))})
    data[0].append(iod)
    pytest.raises(TypeError, stats.nb_stats, data)
    pytest.raises(AttributeError, stats.nb_stats, data, 'charges')
    pytest.raises(TypeError, stats.nb_stats, data[0])
    pytest.raises(KeyError, stats.nb_stats, data[0], 'test')


def test_epol_stats():
    data = dict()
    data[0] = list()
    for i in range(1, 6):
        iod = IOData(extra={'epol': i})
        data[0].append(iod)
    epol_mean, epol_std = stats.epol_stats(data[0])
    assert_allclose(epol_mean, 3.0)
    assert_allclose(epol_std, 1.4142135623730951)


def test_epol_stats_invalid():
    pytest.raises(TypeError, stats.epol_stats)
