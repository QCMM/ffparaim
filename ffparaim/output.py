#!/usr/bin/python3

import pickle

from string import Template
from ffparaim import utils
from openmm import unit


def to_pickle(data):
    """Write MBIS partitioning data to .pkl file."""

    # Open .pkl file.
    with open('ffparaim.pickle', 'wb') as outfile:
        # Write data to file.
        pickle.dump(data, outfile)


def to_csv(lig_structure,
           charges,
           sig,
           eps):
    # Get sigma and epsilon values from system if only atomic charges were derived.
    if sig is None and eps is None:
        # Create a list with sigma values in compatible units.
        sig = [atom.usigma / unit.nanometer for atom in lig_structure.atoms]
        # Create a list with epsilon values in compatible units.
        eps = [atom.uepsilon / unit.kilojoule_per_mole for atom in lig_structure.atoms]
    # Replace values in formated string for every atom.
    dat_lines = [utils.dat_block.format(
                 atom.residue.name,
                 atom.name,
                 charges[i],
                 sig[i],
                 eps[i]) for i, atom in enumerate(lig_structure.atoms)]
    # Create a data block with formated strings.
    dat = "\n".join(dat_lines)
    # Add data to fields.
    fields = {
        "dat_block": dat,
    }
    # Write .dat output file.
    with open("ffparaim.dat", "w") as f:
        # Populate and write file.
        print(Template(utils.dat_template).substitute(fields), file=f)


def to_dat(mean,
           std,
           outfile):
    # Write mean and standard deviation to file.
    with open(outfile, 'w') as f:
        f.write(f'{mean:6f} {std:6f}\n')
        f.close()
