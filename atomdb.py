#!/usr/bin/python3

import json
import ffparaim.qmtools as qmt
from ffparaim import utils
from ffparaim.ffderiv import ForceFieldDerivation


class AtomDB(ForceFieldDerivation):
    """docstring for ForceFieldDerivation."""

    def __init__(self, molecule, method, basis):
        super().__init__()
        self.molecule = molecule
        self.method = method
        self.basis = basis

    def create_table(self):
        atnums = list(set([atom.atomic_number for atom in self.molecule.atoms]))
        rcubed_table = dict()
        for atn in atnums:
            qmt.write_orca_input('uks', atom=utils.elements[atn], method=self.method,
                                 basis=self.basis, qm_mult=utils.mult_table[atn])
            qmt.exec_orca(uks=True)
            self.load_data('orca_uks.molden.input')
            self.set_molgrid(75, 110)
            self.do_partitioning(method='mbis')
            rcubed_table[atn] = self.get_rcubed()[0]
        with open('rcubed_table.json', 'w') as json_out:
            json.dump(rcubed_table, json_out)
        return rcubed_table
