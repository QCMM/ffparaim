#!/usr/bin/python3
import time
import os
import itertools

import ffparaim.mdtools as mdt
import ffparaim.qmtools as qmt

from ffparaim.ffderiv import ForceFieldDerivation
from ffparaim.ffderiv import symmetrize
from ffparaim.ffderiv import get_lj_params
from ffparaim.orcaff import OrcaForceField
from ffparaim.atomdb import AtomDB
from ffparaim import stats
from ffparaim import io
from ffparaim import utils
from iodata import IOData


class FFparAIM(object):
    """docstring for ffparaim."""

    def __init__(self, pdb_file, smiles, qm_charge=0, ligand_selection=':1',
                 receptor_selection=None, n_updates=3, sampling_time=25,
                 total_qm_calculations=100, method='B3LYP', basis='def2-TZVP',
                 forcefield='openff_unconstrained-1.3.0.offxml'):

        self.pdb_file = pdb_file
        self.smiles = smiles
        self.qm_charge = qm_charge  # Total qm charge.
        # Residue index for molecule to calculate charges.
        self.ligand_selection = ligand_selection
        if receptor_selection is not None:
            self.receptor_selection = receptor_selection
        self.n_updates = n_updates
        self.sampling_time = sampling_time  # Sampling time in ns.
        self.total_qm_calculations = total_qm_calculations
        self.method = method
        self.basis = basis
        self.forcefield = forcefield
        self.data = dict()

    def run(self, compl=False, output=True, json=False, charges=True, lj=False, symm=True):

        begin_time = time.time()
        # Get forcefield parameters.
        template, molecule = mdt.get_params(self.smiles, self.forcefield)
        # Create table with r^3 values for isolated atoms.
        if lj:
            atomdb = AtomDB(molecule, self.method, self.basis)
            rcubed_table = atomdb.create_table()
        # Generate serialized OpenMM system.
        ff = mdt.create_forcefield(template)
        pdb = mdt.read_pdb(self.pdb_file)
        ligand_atoms_idx = mdt.get_atoms_idx(pdb, self.ligand_selection)
        polar_h_idx = mdt.get_polar_hydrogens(molecule)
        # Apply restraints for complex simulations.
        restraint = True if compl else False
        receptor_atoms_idx = mdt.get_atoms_idx(pdb, self.receptor_selection) if compl else None
        system = mdt.create_system(ff, pdb, restraint, ligand_atoms_idx, receptor_atoms_idx)
        for update in range(self.n_updates):
            # Store charges and polarization energies for each update.
            self.data[update] = list()
            if update == 0:
                positions = pdb.positions
            # Create an OpenMM simulation object.
            simulation = mdt.setup_simulation(pdb, system, positions, update)
            # Write Orca forcefield file.
            orcaff = OrcaForceField(ff, system, pdb)
            params = orcaff.parse_params()
            orcaff.write_paramsfile(params)
            qm_calculations = int(
                self.total_qm_calculations / self.n_updates) * (update + 1)
            # Start loop to calculate atomic charges from different conformations.
            for i in range(qm_calculations):
                step = int(self.sampling_time * 500000 / qm_calculations)
                simulation.step(step)
                # Calculate charges and polarization energy for current configuration.
                print('Calculating charges ...')
                positions = mdt.get_positions(simulation)
                mdt.image_molecule()
                qm_region = qmt.set_qm_atoms(self.ligand_selection)
                qmt.write_qmmm_pdb(qm_region)
                for inp in ('qmmm', 'pol_corr'):
                    lig = self.ligand_selection if inp is 'pol_corr' else None
                    qmt.write_orca_input(inp, ligand_selection=lig, method=self.method,
                                         basis=self.basis, qm_charge=self.qm_charge)
                qmt.exec_orca()
                # Parameter derivation.
                ffderiv = ForceFieldDerivation()
                ffderiv.load_data('orca_qmmm.molden.input')
                ffderiv.set_molgrid(75, 110)
                ffderiv.do_partitioning(method='mbis')
                # Get data.
                charges = ffderiv.get_charges()
                epol = ffderiv.get_epol()
                rcubed = ffderiv.get_rcubed()
                # Store data.
                self.data[update].append(IOData(atffparams={'charges': charges, 'rcubed': rcubed},
                                                extra={'epol': epol}))
            # Update parameters.
            print('Updating parameters in forcefield ...')
            if charges:
                sig, eps = None, None
                new_charges = stats.nb_stats(self.data[update], charges=True)[0]
                if symm:
                    new_charges = symmetrize(molecule, new_charges)
            if lj:
                new_rcubed = stats.nb_stats(self.data[update], rcubed=True)[0]
                if symm:
                    new_rcubed = symmetrize(molecule, new_rcubed)
                sig, eps = get_lj_params(molecule, new_rcubed, rcubed_table)
            system = mdt.update_params(system, ligand_atoms_idx, polar_h_idx,
                                       charge=new_charges, sigma=sig, epsilon=eps)
        # Save serialized system.
        xml_file = f'{self.pdb_file[:-4]}.xml'
        mdt.serialize_system(system, xml_file)
        if output:
            io.write_output(self.data)
        if json:
            io.write_json(self.data)
        end_time = time.time()
        total_time = utils.get_time(begin_time, end_time)
        print(f'Total time: {round(total_time, 2)} hours')
        return

    def validation(self, overwrite=False, compl=False, **kwargs):

        sampling_time = kwargs['sampling_time'] if 'sampling_time' in kwargs else [
            self.sampling_time]
        n_updates = kwargs['n_updates'] if 'n_updates' in kwargs else [self.n_updates]
        total_qm_calculations = kwargs['total_qm_calculations'] if 'total_qm_calculations' in kwargs else [
            self.total_qm_calculations]
        method = kwargs['method'] if 'method' in kwargs else [self.method]
        basis = kwargs['basis'] if 'basis' in kwargs else [self.basis]
        params_list = [sampling_time, n_updates,
                       total_qm_calculations, method, basis]
        params_comb = list(itertools.product(*params_list))
        for params in params_comb:
            parm_dir = f's_{params[0]}_n_{params[1]}_t_{params[2]}_m_{params[3]}_b_{params[4]}'
            self.data = dict()
            self.sampling_time = params[0]
            self.n_updates = params[1]
            self.total_qm_calculations = params[2]
            self.method = params[3]
            self.basis = params[4]
            io.create_parm_dir(self.pdb_file, parm_dir, overwrite)
            os.chdir(parm_dir)
            self.run(compl, json=True)
            os.chdir('..')
        return
