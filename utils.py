#!/usr/bin/python3
"""Helper variables used for Orca 4.2.1 enviroment execution."""
import os

# List of elements according to atomic number.
elements = [
    None,
    'H', 'He',
    'Li', 'Be',
    'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
    'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub'
]

# Multiplicity for free atoms.
mult_table = {
    1: 2,
    3: 2,
    6: 3,
    7: 4,
    8: 3,
    9: 2,
    10: 0,
    14: 3,
    15: 4,
    16: 3,
    17: 2,
    18: 0,
    35: 2,
    36: 0,
    53: 2
}

# Values in atomic units, taken from:
# Chu, X., & Dalgarno, A. (2004).
# Linear response time-dependent density functional theory for van der Waals coefficients.
# The Journal of Chemical Physics, 121(9), 4083--4088. http://doi.org/10.1063/1.1779576
c6_table = {
    1: 6.5,
    3: 1393.0,
    6: 46.6,
    7: 24.2,
    8: 15.6,
    9: 9.5,
    10: 6.38,
    14: 305.0,
    15: 185.0,
    16: 134.0,
    17: 94.6,
    18: 64.3,
    35: 162.0,
    36: 130.0,
    53: 385
}
alpha_table = {
    1: 4.5,
    3: 164.0,
    6: 12.0,
    7: 7.4,
    8: 5.4,
    9: 3.8,
    10: 2.67,
    14: 37.0,
    15: 25.0,
    16: 19.6,
    17: 15.0,
    18: 11.1,
    35: 20.0,
    36: 16.7,
    53: 35.0
}

# Default template for Orca forcefield.
orcaff_template = """\
$$atoms
${atoms_block}
$$bonds
${bonds_block}
$$angles
${angles_block}
$$dihedrals
${dihedrals_block}\
"""

# String format for Orca forcefield blocks.
orcaff_blocks = {
    'header': '{0:1d} {1:1d} {2:1d}\n',
    'atoms': '{0:6d}   {1:2s}   {2:10.6f}   {3:10.6f}   {4:10.6f}   {5:10.6f}   {6:10.6f}',
    'bonds': '{0:6d} {1:10d} {2:10.6f} {3:10.6f}',
    'angles': '{0:6d} {1:10d} {2:10.6f} {3:10.6f}',
    'dihedrals': '{0:6d}   {1:6d}   {2:6d}   {3:6d}   {4:10.6f}   {5:10.6f}   {6:6d}'
}

# Default template for Orca QM/MM input.
orca_qmmm_template = """\
! QMMM ${method} ${basis} Grid4 TightSCF NOFINALGRID KeepDens
%output PrintLevel Mini Print[ P_Mulliken ] 1 Print[P_AtCharges_M] 1 end
%pal nprocs ${nproc} end
%qmmm
    ORCAFFFilename "ORCAFF.prms"
    Use_QM_InfoFromPDB true     # get QM atoms from pdb file
end
*pdbfile ${qm_charge} ${qm_mult} output_qmmm.pdb
"""

# Default template for Orca polarization correction input.
orca_pol_corr_template = """\
! ${method} ${basis} Grid4 TightSCF NOFINALGRID KeepDens
%output PrintLevel Mini Print[ P_Mulliken ] 1 Print[P_AtCharges_M] 1 end
%pal nprocs ${nproc} end
%coords
    CTyp xyz
    Charge ${qm_charge}
    Mult ${qm_mult}
    Units Angs
    coords
${geometry}
    end
end"""

# Default template for Orca UKS calculation input.
orca_uks_template = """\
! UKS ${method} ${basis} Grid4 TightSCF NOFINALGRID KeepDens
%output PrintLevel Mini Print[ P_Mulliken ] 1 Print[P_AtCharges_M] 1 end
%pal nprocs ${nproc} end
%coords
    CTyp xyz
    Charge ${qm_charge}
    Mult ${qm_mult}
    Units Angs
    coords
        ${atom} 0.0 0.0 0.0
    end
end"""


def get_nproc():
    """Get the number of processes for QM calculation."""
    if 'OMP_NUM_THREADS' in os.environ:
        nproc = int(os.environ['OMP_NUM_THREADS'])
    elif 'SLURM_NTASKS' in os.environ:
        nproc = int(os.environ['SLURM_NTASKS']) - 4
    else:
        nproc = 1
    return nproc


def get_time(begin_time, end_time):
    return (end_time - begin_time) / 3600
