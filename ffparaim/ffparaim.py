#!/usr/bin/python3

import time

import ffparaim.mdtools as mdt
import ffparaim.qmtools as qmt

from ffparaim.ffderiv import ForceFieldDerivation
from ffparaim.ffderiv import symmetrize, normalize_atomic_charges
from ffparaim.ffderiv import get_lj_params
from ffparaim.atomdb import AtomDB
from ffparaim import stats
from ffparaim import output
from ffparaim import utils


class FFparAIM(object):
    """FFparAIM main class for parameter's derivation D-MBIS workflow."""

    def __init__(self,
                 qm_charge=0,
                 ligand_selection='resname LIG',
                 n_updates=3,
                 sampling_time=25,
                 total_qm_calculations=100,
                 method='B3LYP',
                 basis='def2-TZVP'):

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
        # Data dict to store non-bonded parameters per update.
        self.data = dict()
        # SMILES string of the molecule to derive non-bonded parameters.
        self.smiles = None
        # PDB file of the complete system.
        self.pdb_file = None
        # Small molecule force field.
        self.forcefield = None

    def prepare(self,
                smiles,
                pdb_file,
                ff_lig='openff_unconstrained-2.0.0.offxml',
                ff_env=['amber14-all.xml', 'amber14/tip3p.xml']):
        """Prepare method for creating main objects for parameter's
        derivation with D-MBIS workflow based on Open Force Field environment."""
        # SMILES string of the molecule to derive non-bonded parameters.
        self.smiles = smiles
        # PDB file of the complete system.
        self.pdb_file = pdb_file
        # Small molecule force field.
        self.forcefield = ff_lig
        # Separate componentes of the system.
        mdt.separate_components(self.pdb_file, self.ligand_selection)
        # Get small molecule definition.
        molecule = mdt.define_molecule(self.smiles)
        # Create system.
        forcefield = mdt.define_forcefield(self.forcefield)
        lig_structure = mdt.prepare_ligand(molecule, forcefield)
        env_structure = mdt.prepare_enviroment(ff_env)
        system_structure, system = mdt.create_system(lig_structure, env_structure)
        return molecule, lig_structure, system_structure, system

    def run(self,
            molecule,
            lig_structure,
            system_structure,
            system,
            off=False,
            restraint_dict=None,
            csv=True,
            pickle=False,
            charges=True,
            lj=False,
            symm=True,
            exhaustive=False):
        """Run method for executing D-MBIS worklfow."""
        # Start of protocol execution.
        begin_time = time.time()
        self.ligand_atom_list = mdt.get_atoms_idx(system_structure, self.ligand_selection)
        # Derivate Lennard-Jones parameters.
        if lj:
            # Create an atom database object.
            atomdb = AtomDB(molecule, self.method, self.basis)
            # Create table with third radial moment values for isolated atoms.
            rcubed_table = atomdb.create_table()
        for update in range(self.n_updates):
            # Store charges and polarization energies for each update.
            self.data[update] = list()
            if update == 0:
                positions = system_structure.positions
            # Save snapshots in trajectory.
            frames = int(self.sampling_time * 500000 / self.total_qm_calculations)
            # Create an OpenMM simulation object.
            system, simulation = mdt.setup_simulation(system_structure,
                                                      system,
                                                      positions,
                                                      update,
                                                      frames,
                                                      restraint_dict)
            # Save serialized system.
            mdt.save_serialized_system(system, 'system.xml')
            # Write ORCA forcefield file.
            qmt.exec_orca(cmd='mm')
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
                        pol_corr = True if inp == 'pol_corr' else False
                        lig = self.ligand_selection if inp == 'pol_corr' else None
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
                try:
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
                except ValueError:
                    print('Encountered non-finite gradient. Please report this issue on https://github.com/theochem/denspart/issues.')
                    print('Skipping configuration ...')
            # Update parameters.
            print('Updating parameters in forcefield ...')
            if charges:
                sig, eps = None, None
                new_atcharges = stats.nb_stats(self.data[update], 'charges')[0]
                if symm:
                    atcharges = symmetrize(molecule, new_atcharges)
                    norm_atcharges = normalize_atomic_charges(molecule,
                                                              self.qm_charge,
                                                              atcharges)
            if lj:
                new_rcubed = stats.nb_stats(self.data[update], 'rcubed')[0]
                if symm:
                    new_rcubed = symmetrize(molecule, new_rcubed)
                sig, eps = get_lj_params(molecule,
                                         new_rcubed,
                                         rcubed_table)
            system = mdt.update_params(system,
                                       self.ligand_atom_list,
                                       charge=norm_atcharges,
                                       sigma=sig,
                                       epsilon=eps)
            # Save serialized system.
            mdt.save_serialized_system(system, 'system.xml')
        if off:
            smirks_dict = None
            off_ff = mdt.prepare_off_charges(molecule, self.forcefield, norm_atcharges)
            if lj:
                smirks_dict = mdt.prepare_off_lj(molecule, off_ff, sig, eps)
            mdt.save_forcefield(off_ff, smirks_dict, outfile=f'd-mbis_{self.forcefield}')
        if csv:
            output.to_csv(lig_structure, norm_atcharges, sig, eps)
        if pickle:
            output.to_pickle(self.data)
        # Get Polarization Energy value.
        epol_mean, epol_std = stats.epol_stats(self.data[update])
        print(f'Averaged Polarization Energy (kcal/mol) = {epol_mean:6f} +/- {epol_std:6f}')
        output.to_dat(epol_mean, epol_std, 'epol.dat')
        end_time = time.time()
        total_time = utils.get_time(begin_time, end_time)
        print(f'Total time: {round(total_time, 2)} hours')
        return
