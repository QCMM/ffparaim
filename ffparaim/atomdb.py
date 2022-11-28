#!/usr/bin/python3

import json

import ffparaim.qmtools as qmt

from ffparaim import utils
from ffparaim.ffderiv import ForceFieldDerivation
from openff.toolkit.topology import Molecule


class AtomDB(ForceFieldDerivation):
    """AtomDB main class for creating a third radial moment database
    for isoled atoms."""

    def __init__(self, molecule, method, basis):
        super().__init__()
        if isinstance(molecule, Molecule) is True:
            self.molecule = molecule
        else:
            raise ValueError("Variable molecule must be an instance of openff.toolkit.topology.Molecule.")
        self.method = method
        self.basis = basis

    def create_table(self):
        '''Create a table of third radial moment values
        for unique isolated atoms present in the molecule.
        This is necessary to derive Lennard-Jones parameters from
        molecular electron density partitioning.'''

        # Create a list of unique atomic numbers in the molecule.
        atnums = list(set(
            [atom.atomic_number for atom in self.molecule.atoms]))
        # Dict of third radial moment values for isolated atoms.
        rcubed_table = dict()
        # Iterate for every atomic number in list.
        for atn in atnums:
            # Create an ORCA input for unrestricted Kohn-Sham calculation.
            qmt.write_orca_input('uks',
                                 atom=utils.elements[atn],
                                 method=self.method,
                                 basis=self.basis,
                                 qm_mult=utils.mult_table[atn])
            # Run ORCA calculation.
            qmt.exec_orca(cmd='uks')
            # Load data of electron density from molden file.
            iodata = self.load_data('orca_uks.molden.input')
            # Set grid with 150 radial shells and 194 angular points per shell.
            self.set_molgrid(iodata)
            # Apply MBIS electronic density partition method.
            self.do_partitioning(iodata, method='mbis')
            # Add third radial moment value to dict.
            rcubed_table[atn] = self.get_rcubed()[0]
        # Write dict to JSON file.
        with open('rcubed_table.json', 'w') as json_out:
            json.dump(rcubed_table, json_out)
        return rcubed_table
