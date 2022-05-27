#!/usr/bin/python3

import os
import shutil
import pickle

from ffparaim import stats


def create_output_files(charges_file='charges.out', epol_file='epol.out'):
    # define output files
    charges_std_file = f'{charges_file.split(".")[0]}_std.out'
    charges_out = open(charges_file, 'w')
    charges_std_out = open(charges_std_file, 'w')
    epol_out = open(epol_file, 'w')
    charges_out.write('# Atomic charges \n')
    charges_std_out.write('# Atomic charges (Std) \n')
    return charges_out, charges_std_out, epol_out


def create_parm_dir(pdb_file, parm_dir, overwrite):
    # create target directory & all intermediate directories if don't exists
    try:
        os.makedirs(parm_dir)
        print(f'Directory {parm_dir} created')
    except FileExistsError:
        print(f'Directory {parm_dir} already exists')
        if not overwrite:
            print(f'Not overwrite {parm_dir} ...')
            pass
    shutil.copyfile(pdb_file, f'{parm_dir}/{pdb_file}')


def write_charges(charges_mean, charges_std, charges_out, charges_std_out):
    for i, charge in enumerate(charges_mean):
        charges_out.write(f'{charge:.6f}  ')
        charges_std_out.write(f'{charges_std[i]:.6f}  ')
    charges_out.write('\n')
    charges_std_out.write('\n')


def write_epol(epol_mean, epol_std, epol_out):
    epol_out.write(f'{epol_mean:.3f}  {epol_std:.3f}\n')


def write_output(data):
    charges_out, charges_std_out, epol_out = create_output_files()
    # get charges and polarization energy statistics
    for update in data.keys():
        charges_mean, charges_std = stats.nb_stats(data[update], charges=True)
        epol_mean, epol_std = stats.epol_stats(data[update])
        write_charges(charges_mean, charges_std, charges_out, charges_std_out)
        write_epol(epol_mean, epol_std, epol_out)
    charges_out.close
    charges_std_out.close
    epol_out.close


def to_pickle(data):
    with open('ffparaim.pickle', 'wb') as outfile:
        pickle.dump(data, outfile)
