#!/usr/bin/python3

import os
import shutil
import pickle

from string import Template
from ffparaim import utils


def create_parm_dir(pdb_file, parm_dir, overwrite):
    # Create target directory and all intermediate directories.
    try:
        os.makedirs(parm_dir)
        print(f'Directory {parm_dir} created')
    except FileExistsError:
        print(f'Directory {parm_dir} already exists')
        if not overwrite:
            print(f'Not overwrite {parm_dir} ...')
            pass
    shutil.copyfile(pdb_file, f'{parm_dir}/{pdb_file}')


def to_pickle(data):
    with open('ffparaim.pickle', 'wb') as outfile:
        pickle.dump(data, outfile)


def to_dat(lig_structure, charges, sig, eps):
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
