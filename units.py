"""Unit conversion and physical constants.

Atomic units:

    Energy: hartree (Eh)
    Length: bohr (a0)

"""
import scipy.constants as spc

# Lenght
meter: float = 1 / spc.value(u'Bohr radius')
nanometer: float = 1e-9 * meter
# Energy
calmol: float = spc.calorie / spc.value('Avogadro constant') / spc.value('Hartree energy')
kcalmol: float = 1e3 * calmol
kjmol: float = kcalmol / 4.184
