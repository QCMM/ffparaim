#!/usr/bin/python3

import pickle

from string import Template
from ffparaim import utils
from openmm import unit


def to_pickle(data):
    with open('ffparaim.pickle', 'wb') as outfile:
        pickle.dump(data, outfile)


def to_csv(lig_structure, charges, sig, eps):
    if sig is None and eps is None:
        sig = [atom.usigma / unit.nanometer for atom in lig_structure.atoms]
        eps = [atom.uepsilon / unit.kilojoule_per_mole for atom in lig_structure.atoms]
    dat_lines = [utils.dat_block.format(
                 atom.residue.name,
                 atom.name,
                 charges[i],
                 sig[i],
                 eps[i]) for i, atom in enumerate(lig_structure.atoms)]
    dat = "\n".join(dat_lines)
    fields = {
        "dat_block": dat,
    }
    # Write .dat output file.
    with open("ffparaim.dat", "w") as f:
        # Populate and write file.
        print(Template(utils.dat_template).substitute(fields), file=f)
