#!/usr/bin/python3
import os
import shutil
import sys

import subprocess as sp

from parmed.formats.pdb import PDBFile
from string import Template
from ffparaim import utils


def set_qm_atoms(ligand_selection, pdb_file='output_recenter.pdb'):
    pdb = PDBFile.parse(pdb_file)
    atom_list = pdb.atoms
    lig_atoms = pdb[ligand_selection].atoms
    lig_list = [idx for at in lig_atoms for idx, atom in enumerate(atom_list)
                if (at.xx, at.xy, at.xz) == (atom.xx, atom.xy, atom.xz)]
    return lig_list


def write_qmmm_pdb(lig_list, pdb_file='output_recenter.pdb'):
    pdb = PDBFile.parse(pdb_file)
    for atom in pdb.atoms:
        atom.bfactor = 0
        atom.occupancy = 1.0 if atom.idx in lig_list else 0
    pdb.write_pdb('output_qmmm.pdb')


def write_orca_input(orca_inp, atom=None, ligand_selection=None,
                     pdb_file="output_recenter.pdb", method='B3LYP',
                     basis='def2-TZVP', qm_charge=0, qm_mult=1):
    fields = {
        'nproc': utils.get_nproc(),
        'method': method,
        'basis': basis,
        'qm_charge': qm_charge,
        'qm_mult': qm_mult
    }
    if orca_inp is 'qmmm':
        template = utils.orca_qmmm_template
    elif orca_inp is 'pol_corr' and ligand_selection is not None:
        geometry = []
        pdb = PDBFile.parse(pdb_file)
        coords = pdb[ligand_selection].coordinates
        atoms = [utils.elements[atom.atomic_number]
                 for atom in pdb[ligand_selection].atoms]
        for atom, coord in zip(atoms, coords):
            geometry.append(f"{atom:>9} {coord[0]:10.6f} {coord[1]:10.6f} {coord[2]:10.6f}")
        fields["geometry"] = "\n".join(geometry)
        template = utils.orca_pol_corr_template
    elif orca_inp is 'uks':
        fields["atom"] = atom
        template = utils.orca_uks_template
    else:
        raise ValueError("Orca input not implemented.")
    with open(f'orca_{orca_inp}.inp', 'w') as f:
        # Populate files & write input
        print(Template(template).substitute(fields), file=f)


def exec_orca(uks=False, epol=False):
    """Generate commandline for Orca QM calculations."""
    cmdline = "cd " + os.getcwd() + "; "
    # Check if orca is in the PATH
    orca_cmd = shutil.which('orca')
    if orca_cmd is not None:
        if uks:
            cmdline += orca_cmd + " orca_uks.inp > orca_uks.out; "
            cmdline += orca_cmd + "_2mkl" + " orca_uks -molden > /dev/null "
        else:
            cmdline += orca_cmd + " orca_qmmm.inp > orca_qmmm.out; cp orca_qmmm.gbw orca_pol_corr.gbw; "
            cmdline += orca_cmd + "_2mkl" + " orca_qmmm -molden > /dev/null; "
            if epol:
                cmdline += orca_cmd + " orca_pol_corr.inp > orca_pol_corr.out "
        proc = sp.Popen(args=cmdline, shell=True)
        proc.wait()
        if proc.returncode != 0:
            sys.exit(proc.returncode)
