#!/usr/bin/python3
import os
import shutil
import sys

import subprocess as sp

from parmed.formats.pdb import PDBFile
from string import Template
from ffparaim import utils


def set_qm_atoms(ligand_selection, pdb_file='output_recenter.pdb'):
    """Define atom indices that are part of the QM region based on the user's selection."""

    # Read PDB file.
    pdb = PDBFile.parse(pdb_file)
    # Create a list of all atoms in the PDB file.
    atom_list = pdb.atoms
    # Create a list for atom selection.
    lig_atoms = pdb[ligand_selection].atoms
    # Create a list of atom indices based on the full atom list order.
    lig_list = [idx for at in lig_atoms for idx, atom in enumerate(atom_list)
                if (at.xx, at.xy, at.xz) == (atom.xx, atom.xy, atom.xz)]
    return lig_list


def write_qmmm_pdb(lig_list, pdb_file='output_recenter.pdb'):
    """Write a PDB File with B-factor and occupancy columns compatible with
    ORCA QM/MM calculation."""

    # Read PDB file.
    pdb = PDBFile.parse(pdb_file)
    # Iterate over every atom in the system.
    for atom in pdb.atoms:
        # Replace B-factor values for all atoms.
        atom.bfactor = 0
        # Replace occupancy values for atom in QM region, based on ligand'sselection.
        atom.occupancy = 1.0 if atom.idx in lig_list else 0
    # Write modified PDB file.
    pdb.write_pdb('output_qmmm.pdb')


def write_orca_input(orca_inp,
                     atom=None,
                     ligand_selection=None,
                     pdb_file="output_recenter.pdb",
                     method='B3LYP',
                     basis='def2-TZVP',
                     qm_charge=0,
                     qm_mult=1):
    """Write an ORCA input files for obtaining the molecular (or atomic)
    electron density and polarization correction calculations."""

    # Define a field dictionary with ORCA parameters.
    fields = {
        'nproc': utils.get_nproc(),
        'method': method,
        'basis': basis,
        'qm_charge': qm_charge,
        'qm_mult': qm_mult
    }
    # If ORCA run is a Kohn-Sham DFT QM/MM calculation.
    if orca_inp == 'qmmm':
        # Use Kohn-Sham DFT QM/MM calculation template.
        template = utils.orca_qmmm_template
    # If ORCA run is for polarization correction calculation and a ligand selection is defined.
    elif orca_inp == 'pol_corr' and ligand_selection is not None:
        # Empty list for populate with geometry's field lines.
        geometry = []
        # Read PDB file.
        pdb = PDBFile.parse(pdb_file)
        # Extract coordinates for ligand.
        coords = pdb[ligand_selection].coordinates
        # Generate a list with every element in the ligand.
        atoms = [utils.elements[atom.atomic_number]
                 for atom in pdb[ligand_selection].atoms]
        # Iterate over every atom elements and coordinates.
        for atom, coord in zip(atoms, coords):
            # Add string with formated string for geometry field.
            geometry.append(f"{atom:>9} {coord[0]:10.6f} {coord[1]:10.6f} {coord[2]:10.6f}")
        # Assign the geometry field variable with string block.
        fields["geometry"] = "\n".join(geometry)
        # Use polarization correction template.
        template = utils.orca_pol_corr_template
    # If ORCA run is an unrestricted Kohn-Sham DFT QM/MM calculation.
    elif orca_inp is 'uks':
        # Define element for atom.
        fields["atom"] = atom
        # Use unrestricted Kohn-Sham DFT QM/MM calculation template.
        template = utils.orca_uks_template
    # If ORCA input are not in the templates (qmmm, pol_corr, uks).
    else:
        # Raise error.
        raise ValueError("ORCA input not implemented.")
    # Open ORCA input template file.
    with open(f'orca_{orca_inp}.inp', 'w') as f:
        # Populate file and write input.
        print(Template(template).substitute(fields), file=f)


def exec_orca(cmd='ks', epol=False):
    """Generate command line for ORCA calculations."""

    # Change to work directory command line.
    cmdline = "cd " + os.getcwd() + "; "
    # Check if ORCA is in the path.
    orca_cmd = shutil.which('orca')
    # If ORCA is in the path.
    if orca_cmd is not None:
        # If ORCA run is an unrestricted Kohn-Sham DFT calculation.
        if cmd == 'uks':
            # ORCA run for unrestricted Kohn-Sham DFT.
            cmdline += orca_cmd + " orca_uks.inp > orca_uks.out; "
            # Create a molden file from ORCA output.
            cmdline += orca_cmd + "_2mkl" + " orca_uks -molden > /dev/null "
        # If ORCA run is a pre-processing step for QM/MM calculation.
        elif cmd == 'mm':
            # Create an ORCA forcefield file.
            cmdline += orca_cmd + "_mm" + " -convff -OPENMM system.xml "
        # If ORCA run is a Kohn-Sham DFT calculation
        else:
            # ORCA run for Kohn-Sham DFT and copy Gerbview file for polarization correction calculation.
            cmdline += orca_cmd + " orca_qmmm.inp > orca_qmmm.out; cp orca_qmmm.gbw orca_pol_corr.gbw; "
            # Create a molden file from ORCA output.
            cmdline += orca_cmd + "_2mkl" + " orca_qmmm -molden > /dev/null; "
            # If polarization correction calculation is needed.
            if epol:
                # ORCA run for polarization correction.
                cmdline += orca_cmd + " orca_pol_corr.inp > orca_pol_corr.out "
        # Run command lines.
        proc = sp.Popen(args=cmdline, shell=True)
        # Wait for process to finish.
        proc.wait()
        # If something goes wrong.
        if proc.returncode != 0:
            # Return error code.
            sys.exit(proc.returncode)
