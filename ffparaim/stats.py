#!/usr/bin/python3
import numpy as np


def nb_stats(data, charges=False, rcubed=False):

    param = 'rcubed' if rcubed else 'charges'
    params = np.array([data_update.atffparams[param] for data_update in data])
    mean, std = params.mean(
        axis=0, dtype=np.float64), params.std(axis=0, dtype=np.float64)
    return mean, std


def epol_stats(data_update):

    epol = np.array([element.extra['epol'] for element in data_update])
    epol_mean, epol_std = epol.mean(), epol.std()
    return epol_mean, epol_std
