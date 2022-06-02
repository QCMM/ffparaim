#!/usr/bin/python3

import time
import os
import itertools

import ffparaim.mdtools as mdt
import ffparaim.qmtools as qmt

from ffparaim.ffderiv import ForceFieldDerivation
from ffparaim.ffderiv import symmetrize, normalize_atomic_charges
from ffparaim.ffderiv import get_lj_params
from ffparaim.orcaff import OrcaForceField
from ffparaim.atomdb import AtomDB
from ffparaim import stats
from ffparaim import output
from ffparaim import utils


class FFparAIM(object):
    """docstring for FFparAIM."""

    def __init__(self,
                 pdb_file,
                 smiles,
                 qm_charge=0,
                 ligand_selection=':1',
                 n_updates=3,
                 sampling_time=25,
                 total_qm_calculations=100,
                 method='B3LYP',
                 basis='def2-TZVP',
                 forcefield='openff_unconstrained-2.0.0.offxml'):

        # PDB file of the complete system.
        self.pdb_file = pdb_file
        # SMILES string of the molecule to derive non-bonded parameters.
        self.smiles = smiles
        # Total charge for the ligand in the QM subsystem.
        self.qm_charge = qm_charge
        # Ligand AMBER mask residue index.
        self.ligand_selection = ligand_selection
        # Ligand atoms index.
        self.ligand_atom_list = None
        # Number of updates for non-bonded parameters.
        self.n_updates = n_updates
        # Sampling time for trajectories in nanoseconds.
        self.sampling_time = sampling_time
        # Number of total QM calculations for polarized electron density.
        self.total_qm_calculations = total_qm_calculations
        # DFT functional to use in QM calculations.
        self.method = method
        # Basis set for QM calculations.
        self.basis = basis
        # Small molecule force field.
        self.forcefield = forcefield
        # Data dict to store non-bonded parameters per update.
        self.data = dict()

    def run(self,
            restraint_dict=None,
            dat=True,
            pickle=False,
            charges=True,
            lj=False,
            symm=True,
            exhaustive=False):
        '''Run documentation.'''

        # Start of protocol execution.
        begin_time = time.time()
        # Separate componentes of the system.
        mdt.separate_components(self.pdb_file, self.ligand_selection)
        # Get small molecule definition.
        molecule = mdt.define_molecule(self.smiles)
        # Derivate Lennard-Jones parameters.
        if lj:
            # Create an atom database object.
            atomdb = AtomDB(molecule, self.method, self.basis)
            # Create table with third radial moment values for isolated atoms.
            rcubed_table = atomdb.create_table()
        # Create an OpenMM ForceField object from small molecule template generator.
        # ff = mdt.create_forcefield(template)
        # Read system coordinates from PDB file.
        # pdb = mdt.read_pdb(self.pdb_file)
        # Polar hydrogens index for ligand.
        # polar_h_idx = mdt.get_polar_hydrogens(molecule)
        # Generate serialized OpenMM system.
        lig_structure = mdt.prepare_ligand(molecule, self.forcefield)
        env_structure = mdt.prepare_enviroment()
        system_structure, system = mdt.create_system(lig_structure, env_structure)
        self.ligand_atom_list = mdt.get_atoms_idx(system_structure, self.ligand_selection)
        for update in range(self.n_updates):
            # Store charges and polarization energies for each update.
            self.data[update] = list()
            if update == 0:
                positions = system_structure.positions
            # Create an OpenMM simulation object.
            simulation = mdt.setup_simulation(system_structure,
                                              system,
                                              positions,
                                              update,
                                              restraint_dict,
                                              self.ligand_atom_list)
            # Write ORCA forcefield file.
            orcaff = OrcaForceField(system_structure, system)
            params = orcaff.parse_params()
            orcaff.write_paramsfile(params)
            qm_calculations = int(
                self.total_qm_calculations / self.n_updates) * (update + 1)
            # Starting loop to calculate atomic charges and third radial moments.
            for i in range(qm_calculations):
                step = int(self.sampling_time * 500000 / qm_calculations)
                simulation.step(step)
                # Calculate non-bonded parameters and polarization energy for current configuration.
                print('Calculating non-bonded parameters ...')
                positions = mdt.get_positions(simulation)
                mdt.image_molecule()
                qm_region = qmt.set_qm_atoms(self.ligand_selection)
                qmt.write_qmmm_pdb(qm_region)
                # QM/MM calculation.
                print('QM/MM calculation in progress ...')
                if exhaustive:
                    print('Exhaustive Polarization Energy calculation ...')
                    for inp in ('qmmm', 'pol_corr'):
                        pol_corr = True if inp is 'pol_corr' else False
                        lig = self.ligand_selection if inp is 'pol_corr' else None
                        qmt.write_orca_input(inp,
                                             ligand_selection=lig,
                                             method=self.method,
                                             basis=self.basis,
                                             qm_charge=self.qm_charge)
                        qmt.exec_orca(epol=pol_corr)
                else:
                    qmt.write_orca_input('qmmm',
                                         ligand_selection=None,
                                         method=self.method,
                                         basis=self.basis,
                                         qm_charge=self.qm_charge)
                    qmt.exec_orca()
                    # Polarization energy calculation.
                    if update + 1 == self.n_updates:
                        print('Polarization Energy calculation ...')
                        qmt.write_orca_input('pol_corr',
                                             ligand_selection=self.ligand_selection,
                                             method=self.method,
                                             basis=self.basis,
                                             qm_charge=self.qm_charge)
                        qmt.exec_orca(epol=True)
                # Parameter derivation.
                ffderiv = ForceFieldDerivation()
                iodata = ffderiv.load_data('orca_qmmm.molden.input')
                ffderiv.set_molgrid(iodata)
                ffderiv.do_partitioning(iodata, method='mbis')
                # Get data.
                charge = ffderiv.get_charges()
                rcubed = ffderiv.get_rcubed()
                if exhaustive:
                    epol = ffderiv.get_epol()
                else:
                    epol = ffderiv.get_epol() if update + 1 == self.n_updates else None
                # Store data.
                iodata.atffparams = {'charges': charge.tolist(),
                                     'rcubed': rcubed.tolist()}
                iodata.extra = {'epol': epol}
                self.data[update].append(iodata)
            # Update parameters.
            print('Updating parameters in forcefield ...')
            if charges:
                sig, eps = None, None
                new_atcharges = stats.nb_stats(self.data[update], charges=True)[0]
                if symm:
                    atcharges = symmetrize(molecule, new_atcharges)
                    norm_atcharges = normalize_atomic_charges(molecule,
                                                              self.qm_charge,
                                                              atcharges)
                    # print(new_charges)
            if lj:
                new_rcubed = stats.nb_stats(self.data[update], rcubed=True)[0]
                if symm:
                    new_rcubed = symmetrize(molecule, new_rcubed)
                    # print(new_rcubed)
                sig, eps = get_lj_params(molecule,
                                         new_rcubed,
                                         rcubed_table)
            system = mdt.update_params(system,
                                       self.ligand_atom_list,
                                       charge=norm_atcharges,
                                       sigma=sig,
                                       epsilon=eps)
        # Save serialized system.
        xml_file = f'{self.pdb_file[:-4]}.xml'
        mdt.serialize_system(system, xml_file)
        if output:
            output.to_dat(lig_structure, norm_atcharges, sig, eps)
        if pickle:
            output.to_pickle(self.data)
        # Get Polarization Energy value.
        epol_mean, epol_std = stats.epol_stats(self.data[update])
        print(f'Averaged Polarization Energy (kcal/mol) = {epol_mean:6f} +/- {epol_std:6f}')
        end_time = time.time()
        total_time = utils.get_time(begin_time, end_time)
        print(f'Total time: {round(total_time, 2)} hours')
        return

    def validation(self, overwrite=False, restraint_dict=None, **kwargs):

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
            output.create_parm_dir(self.pdb_file, parm_dir, overwrite)
            os.chdir(parm_dir)
            self.run(restraint_dict, lj=True, pickle=True)
            os.chdir('..')
        return
