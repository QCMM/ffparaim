#!/usr/bin/python3

import json

import ffparaim.qmtools as qmt

from ffparaim import utils
from ffparaim.ffderiv import ForceFieldDerivation


class AtomDB(ForceFieldDerivation):
    """docstring for AtomDB."""

    def __init__(self, molecule, method, basis):
        super().__init__()
        self.molecule = molecule
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
            qmt.exec_orca(uks=True)
            # Load data of electron density from molden file.
            self.load_data('orca_uks.molden.input')
            # Set grid with 75 radial shells and 110 angular points per shell.
            self.set_molgrid(75, 110)
            # Apply MBIS electronic density partition method.
            self.do_partitioning(method='mbis')
            # Add third radial moment value to dict.
            rcubed_table[atn] = self.get_rcubed()[0]
        # Write dict to JSON file.
        with open('rcubed_table.json', 'w') as json_out:
            json.dump(rcubed_table, json_out)
        return rcubed_table
