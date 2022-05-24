#!/usr/bin/python3
import numpy as np


def nb_stats(datas, charges=False, rcubed=False):

    param = 'rcubed' if rcubed else 'charges'
    params = np.array([data.atffparams[param] for data in datas])
    mean, std = params.mean(
        axis=0, dtype=np.float64), params.std(axis=0, dtype=np.float64)
    return mean, std


def epol_stats(datas):

    epol = np.array([data.extra['epol'] for data in datas])
    epol_mean, epol_std = epol.mean(), epol.std()
    return epol_mean, epol_std
