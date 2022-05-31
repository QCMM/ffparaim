#!/usr/bin/python3

from string import Template
from ffparaim import utils


class OrcaForceField(object):
    """docstring for OrcaForceField."""

    def __init__(self, lig_structure, env_structure):
        '''_system = ff.createSystem(pdb.topology,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=1 * unit.nanometer,
                                  rigidWater=False)
        '''
        self.system_toppar = lig_structure + env_structure
        self._template = utils.orcaff_template

    def parse_params(self):
        # Non bonded LJ 1-4 scaling factor
        SCNB = 2.0
        # Add atom numbers
        for i, atom in enumerate(self.system_toppar.atoms):
            atom.number = i + 1
        # Write atom numbers, type, charge and LJ parameters.
        # Note the scaling of 1-4 parameters are performed here and that the LJ parameteres
        # are chosen to fit the CHARMM energy expression (sigma to rmin and epsilon * -1)
        atoms_header = utils.orcaff_blocks["header"].format(len(self.system_toppar.atoms), 1, 4)
        atoms_lines = [utils.orcaff_blocks["atoms"].format(
            atom.number,
            utils.elements[atom.atomic_number],
            atom.charge,
            atom.epsilon * -1,
            atom.rmin * 2,
            (atom.epsilon * -1) / SCNB,
            atom.rmin * 2) for atom in self.system_toppar.atoms]
        # Amber, GROMACS and CHARMM format are identical for bonds
        bonds_header = utils.orcaff_blocks["header"].format(len(self.system_toppar.bonds), 2, 2)
        bonds_lines = [utils.orcaff_blocks["bonds"].format(
            bond.atom1.number,
            bond.atom2.number,
            bond.type.req,
            bond.type.k) for bond in self.system_toppar.bonds]
        # Amber, GROMACS and CHARMM format identical for angles
        angles_header = utils.orcaff_blocks["header"].format(len(self.system_toppar.angles), 3, 2)
        angles_lines = [utils.orcaff_blocks["angles"].format(
            angle.atom1.number,
            angle.atom2.number,
            angle.atom3.number,
            angle.type.theteq,
            angle.type.k) for angle in self.system_toppar.angles]
        # Amber and CHARMM format identical for dihedrals.
        # GROMACS could have more than one dihedral type description for the same atoms.
        # This is the reason for use a try except statement.
        dihedrals_header = utils.orcaff_blocks["header"].format(
            len(self.system_toppar.dihedrals), 4, 3)
        dihedrals_lines = [utils.orcaff_blocks["dihedrals"].format(
            dihedral.atom1.number,
            dihedral.atom2.number,
            dihedral.atom3.number,
            dihedral.atom4.number,
            dihedral.type.phase,
            dihedral.type.phi_k,
            dihedral.type.per) for dihedral in self.system_toppar.dihedrals]
        # Assign str objects for writing.
        atoms = atoms_header + "\n".join(atoms_lines)
        bonds = bonds_header + "\n".join(bonds_lines)
        angles = angles_header + "\n".join(angles_lines)
        dihedrals = dihedrals_header + "\n".join(dihedrals_lines)
        return (atoms, bonds, angles, dihedrals)

    def write_paramsfile(self, params, template=utils.orcaff_template):
        # Generate fields dict.
        atoms, bonds, angles, dihedrals = params
        fields = {
            "atoms_block": atoms,
            "bonds_block": bonds,
            "angles_block": angles,
            "dihedrals_block": dihedrals
        }

        # Write Orca forcefield file.
        with open("ORCAFF.prms", "w") as f:
            # Populate files & write input
            print(Template(template).substitute(fields), file=f)
