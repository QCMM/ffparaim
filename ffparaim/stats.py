#!/usr/bin/python3
import numpy as np


def nb_stats(data, charges=False, rcubed=False):
    """Obtain the arimetic mean and standard deviation for non-bonded parameters."""

    # Assign respective value if correspond to atomic charges or effective volume.
    param = 'rcubed' if rcubed else 'charges'
    # Create an array for non-bonded parameters data of all updates.
    params = np.array([data_update.atffparams[param] for data_update in data])
    # Calculate arimetic mean and standard deviation.
    mean, std = params.mean(
        axis=0, dtype=np.float64), params.std(axis=0, dtype=np.float64)
    return mean, std


def epol_stats(data_update):
    """Obtain the arimetic mean and standard deviation for polarization correction energy."""

    # Create an array for polarization correction energy for every element in the update.
    epol = np.array([element.extra['epol'] for element in data_update])
    # Calculate arimetic mean and standard deviation.
    epol_mean, epol_std = epol.mean(), epol.std()
    return epol_mean, epol_std
