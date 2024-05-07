#!/usr/bin/python3
import warnings
import mdtraj
import sys

import parmed as pmd
import pytraj
import numpy as np

from ffparaim.restraints import set_restraints
from openmm import openmm
from openmm import app
from openmm import unit
from openff import units
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler
from rdkit.Chem import AllChem, AddHs

# Avoid warnings.
warnings.filterwarnings('ignore')


def separate_components(pdb_file,
                        ligand_selection):
    """Create PDB files of ligand and molecular environment from PDB file
    input."""

    # Load PDB input file.
    pdb = pmd.load_file(pdb_file)
    # Write PDB file for ligand only.
    pdb[ligand_selection].write_pdb('lig.pdb')
    # Write PDB file for molecular environment.
    pdb[f'!{ligand_selection}'].write_pdb('env.pdb')


def fix_conect(lig_pdb_file='lig.pdb'):
    lig = pytraj.load(lig_pdb_file)
    lig.save('lig.pdb', options='pdbv3')


def define_molecule(smiles,
                    lig_pdb_file='lig.pdb'):
    """Create an openff-toolkit Molecule object from a RDKit template using a
    SMILES string and a PDB file of the protonated ligand."""

    # Fix CONNECT records.
    fix_conect(lig_pdb_file)
    # Create RDKit template from SMILES.
    template = AllChem.MolFromSmiles(smiles)
    # Add hydrogens.
    template = AddHs(template)
    # Read ligand PDB file with hydrogens.
    pdb = AllChem.MolFromPDBFile(lig_pdb_file, removeHs=False)
    # Assign atom conectivity and bond orders.
    rdmol = AllChem.AssignBondOrdersFromTemplate(template, pdb)
    return Molecule.from_rdkit(rdmol)


def define_forcefield(ff):
    """Create an openff-toolkit ForceField object from offxml filename."""

    return ForceField(ff)


def prepare_ligand(molecule,
                   forcefield,
                   lig_pdb_file='lig.pdb'):
    """Create a parametrized ParmEd Structure object for the ligand using an
    opeff-toolkit Molecule and ForceField object with a PDB file of the
    protonated ligand.."""

    # Read ligand PDB file.
    lig_pdb = app.PDBFile(lig_pdb_file)
    # Create ligand topology.
    off_topology = Topology.from_openmm(lig_pdb.getTopology(),
                                        unique_molecules=[molecule])
    # Create OpenMM System object for the ligand.
    lig_system = forcefield.create_openmm_system(off_topology)
    # Create ligand parametrized ParmEd Structure object.
    lig_structure = pmd.openmm.load_topology(lig_pdb.getTopology(),
                                             lig_system,
                                             xyz=lig_pdb.getPositions())
    return lig_structure


def prepare_enviroment(ff_env=['amber14-all.xml', 'amber14/tip3p.xml'],
                       env_pdb_file='env.pdb'):
    """Create a ParmEd Structure object of the molecular environment using
    OpenMM system building."""

    # Create OpenMM ForceField object.
    omm_ff = app.ForceField()
    # Iterate over forcefield files and load them.
    for ff in ff_env:
        omm_ff.loadFile(ff)
    # Create OpenMM PDBFile object.
    env_pdb = app.PDBFile(env_pdb_file)
    # Get topology from PDB file.
    omm_topology = env_pdb.getTopology()
    # Create OpenMM System object for molecular environment compatible with ParmEd.
    env_system = omm_ff.createSystem(omm_topology,
                                     rigidWater=False)
    # Create environment parametrized ParmEd Structure object.
    env_structure = pmd.openmm.load_topology(omm_topology,
                                             env_system,
                                             xyz=env_pdb.getPositions())
    return env_structure


def create_system(lig_structure,
                  env_structure):
    """Create OpenMM system from ParmEd Structure objects of ligand
    and molecular environment."""

    # Create system ParmEd Structure Object.
    system_structure = lig_structure + env_structure
    # Create OpenMM System Object.
    system = system_structure.createSystem(nonbondedMethod=app.PME,
                                           nonbondedCutoff=1 * unit.nanometer,
                                           constraints=app.HBonds)
    # Add Barostat.
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 298.15 * unit.kelvin))
    return system_structure, system


def save_serialized_system(system,
                           xml_file):
    """Save an OpenMM System Object in XML format."""
    # Serialize system.

    system_serialized = openmm.XmlSerializer.serialize(system)
    # Save XML file.
    with open(xml_file, 'w') as f:
        f.write(system_serialized)


def prepare_off_charges(molecule,
                        ff,
                        charges):
    """Include molecule's partial charges with normalized D-MBIS atomic charges on
    OpenForceField force field object using LibraryCharges."""

    off_ff = ForceField(ff)
    molecule.partial_charges = charges * units.unit.elementary_charge
    library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(molecule)
    off_ff["LibraryCharges"].add_parameter(parameter=library_charge_type, allow_duplicate_smirks=True)
    return off_ff


def prepare_off_lj(molecule,
                   off_ff,
                   sig,
                   eps):
    """Generate a SMIRKS-based dictionary that contains sigma and epsilon
    Lennard-Jones parameters for every pattern in the molecule."""

    # Get applied SMIRKS patterns from force field.
    applied_parameters = off_ff.label_molecules(molecule.to_topology())
    # Create a list with the applied patterns.
    patterns = [value.smirks for value in applied_parameters[0]['vdW'].values()]
    # Create SMIRKS dict with Lennard-Jones parameters.
    smirks_dict = dict()
    # Iterate over every pattern.
    for index, pattern in enumerate(patterns):
        # If pattern is not in dictionary, create a key for it and add the value.
        if pattern not in smirks_dict.keys():
            smirks_dict[pattern] = {'sigma': [sig[index]],
                                    'epsilon': [eps[index]]}
        # If exists, add the value to the list in the respective key.
        else:
            smirks_dict[pattern]['sigma'].append(sig[index])
            smirks_dict[pattern]['epsilon'].append(eps[index])
    return smirks_dict


def save_forcefield(off_ff,
                    smirks_dict=None,
                    outfile='d-mbis.offxml'):
    """Save a modified version of OpenForceField force field object with derived
    parameters."""

    if smirks_dict is not None:
        # Replace Lennard-Jones parameters for every SMIRKS.
        for pattern in smirks_dict.keys():
            # Rmin / 2.
            rmin_h = (np.array(smirks_dict[pattern]['sigma']).mean()) / (2 ** (5.0 / 6.0))
            # Epsilon.
            eps = np.array(smirks_dict[pattern]['epsilon']).mean()
            # Define the new Lennard-Jones parameters for pattern.
            vdw_handler = off_ff["vdW"].parameters[pattern]
            # Replace Rmin / 2 and epsilon values.
            vdw_handler.rmin_half = round(rmin_h, 6) * unit.nanometer
            vdw_handler.epsilon = round(eps, 6) * unit.kilojoule_per_mole
    # Save offxml file.
    off_ff.get_parameter_handler('Electrostatics').periodic_potential = 'PME'
    off_ff.to_file(outfile)


def setup_simulation(system_structure,
                     system,
                     positions,
                     update,
                     frames,
                     restraint_dict=None,
                     ligand_atom_list=None):
    """Setup the OpenMM simulation with the current topology and the input coordinates
     or the current positions, in an update-depending form."""

    # Apply restraint force for complex simulations.
    if restraint_dict is not None:
        system = set_restraints(system_structure.topology,
                                system,
                                positions,
                                restraint_dict,
                                ligand_atom_list)
    # Create an integrator instance.
    integrator = openmm.LangevinIntegrator(298.15 * unit.kelvin,
                                           1 / unit.picosecond,
                                           0.002 * unit.picoseconds)
    # Generate an OpenMM Simulation object.
    simulation = app.Simulation(system_structure.topology,
                                system,
                                integrator)
    # Extract simulation info for saved frames.
    simulation.reporters.append(app.StateDataReporter(sys.stdout,
                                                      frames,
                                                      step=True,
                                                      potentialEnergy=True,
                                                      temperature=True,
                                                      density=True))
    # Save frame in DCD file every 'frames' step.
    simulation.reporters.append(app.DCDReporter(f'traj_{update}.dcd', frames))
    # Set particle positions.
    simulation.context.setPositions(positions)
    # Energy minimization.
    simulation.minimizeEnergy()
    return system, simulation


def get_positions(simulation):
    """Extract coordinates from a configuration in the trayectory, including the
    periodic box vector and save them in a PDB file. Also return the particle's
    positions to continue the simulation."""

    # Get positions from simulation.
    positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=True).getPositions()
    # Extract box vectors information.
    _box_vectors = simulation.context.getState().getPeriodicBoxVectors()
    # Set box vectors in topology just to save it in the PDB file.
    simulation.topology.setPeriodicBoxVectors(_box_vectors)
    # Write PDB file.
    app.PDBFile.writeFile(simulation.topology, positions,
                          open('output.pdb', 'w'))
    return positions


def image_molecule(pdbfile='output.pdb'):
    '''Move system to the center of the box, wrapping atoms to the periodic
    unit cell.'''

    # Create MDTraj Trayectory object.
    traj = mdtraj.load(pdbfile)
    # Apply system recenter.
    try:
        traj.image_molecules()
    # If is it not possible, use first molecule of the system as reference.
    except ValueError:
        traj.image_molecules(anchor_molecules=[traj.top.find_molecules()[0]],
                             other_molecules=traj.top.find_molecules()[1:])
    # Transform positions in compatible units for PDB writing.
    pos = mdtraj.utils.in_units_of(traj.xyz[0],
                                   traj._distance_unit,
                                   "angstroms")
    # Same for unit cell vectors.
    ucl = mdtraj.utils.in_units_of(traj.unitcell_lengths[0],
                                   traj._distance_unit,
                                   "angstroms")
    # Create PDBTrajectoryFile object.
    pdb_recentered = mdtraj.formats.PDBTrajectoryFile("output_recenter.pdb", "w")
    # Write coordinates and box vectors to PDB file.
    pdb_recentered.write(pos,
                         traj.top,
                         modelIndex=None,
                         unitcell_lengths=ucl,
                         unitcell_angles=traj.unitcell_angles[0])


def get_atoms_idx(system_structure,
                  residue_selection):
    """Generate an atom index lists for a residue selection."""

    # Generate an atom atom list for entire system.
    atom_list = system_structure.atoms
    # Generate a residue list for selection.
    res = system_structure[residue_selection].atoms
    # Map atom ids for selected residues.
    atoms_idx = [idx for at in res for idx, atom in enumerate(atom_list) if (
        at.residue.number, at.name) == (atom.residue.number, atom.name)]
    return atoms_idx


def update_params(system,
                  ligand_atoms_idx,
                  charge=None,
                  sigma=None,
                  epsilon=None):
    """Update non-bonded parameters for the system replacing the values in
    OpenMM System object."""

    # Iterate over every force in system.
    for force in system.getForces():
        # If force correspond to NonbondedForce class.
        if force.getName() == 'NonbondedForce':
            # Iterate for every atom index in the list.
            for i, idx in enumerate(ligand_atoms_idx):
                # Get atomic charge value from OpenMM System object or from MBIS partitioning.
                q = force.getParticleParameters(idx)[0] if charge is None else round(charge[i], 6)
                # Get sigma value from OpenMM System object or from MBIS partitioning.
                sig = force.getParticleParameters(idx)[1] if sigma is None else round(sigma[i], 6)
                # Get epsilon value from OpenMM System object or from MBIS partitioning.
                eps = force.getParticleParameters(idx)[2] if epsilon is None else round(epsilon[i], 6)
                # Replace parameters in system.
                force.setParticleParameters(idx, charge=q, sigma=sig, epsilon=eps)
    return system
