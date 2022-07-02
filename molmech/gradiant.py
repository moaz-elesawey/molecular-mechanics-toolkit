from math import cos, radians, sqrt, pi, sin

from molmech.consts import const
from .mymath import euc

def get_bond_gradiant(bond):
    return 2.0 * bond.kb * (bond.r() - bond.eq_r)

def get_angle_gradiant(angle):
    return 2.0 * angle.ka * radians((angle.theta() - angle.eq_theta))

def get_torsion_gradiant(torsion):
    return torsion.nfolds * 0.5 * torsion.vn * (
        sin(torsion.nfolds*radians(torsion.dihedral()) - radians(torsion.gamma))
    )

def get_oof_gradiant(oof):
    return 0.5 * oof.vn * sin(radians(oof.phi()) - pi)

class Gradiant:
    def __init__(self, atoms, bonds, angles, torsions, oofs) -> None:
        self.atoms          = atoms
        self.bonds          = bonds
        self.angles         = angles
        self.torsions       = torsions
        self.oofs            = oofs

        self.total_gradiant   = self.tot_gradiant()

    def bonds_gradiant(self):
        """compute the enery of bonds"""
        e = 0.0

        for b in self.bonds:
            e += get_bond_gradiant(b)

        return e

    def angles_gradiant(self):
        """compute the enery of angles"""
        e = 0.0

        for a in self.angles:
            e += get_angle_gradiant(a)

        return e

    def cross_terms(self):
        """compute the cross terms gradiant"""
        e = 0.0

        for b in self.bonds:
            for a in self.angles:
                e += a.ka * b.kb * (b.r() - b.eq_r) * (radians(a.theta() - a.eq_theta))

        return e

    def torsions_gradiant(self):
        """compute the enery of torsions"""
        e = 0.0

        for t in self.torsions:
            e += get_torsion_gradiant(t)

        return e

    def outofplanes_gradiant(self):
        """compute the enery of out of planes"""
        e = 0.0

        for o in self.oofs:
            e += get_oof_gradiant(o)

        return e

    def vdws_gradiant(self):
        """compute the enery of vdws"""
        e = 0.0

        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                atom_i = self.atoms[i]
                atom_j = self.atoms[j]

                eps_ij = sqrt(atom_i.eps * atom_j.eps)
                ro_ij = 0.5 * (atom_i.ro + atom_j.ro)
                r_ij = euc(atom_i.pos, atom_j.pos)

                e += -12.0 *eps_ij * ((ro_ij/r_ij)**13 - 2*(ro_ij/r_ij)**6)

        return e

    def elsts_gradiant(self):
        """compute the enery of electrostatic"""
        e = 0.0

        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                atom_i = self.atoms[i]
                atom_j = self.atoms[j]
                r_ij = euc(atom_i.pos, atom_j.pos)

                e += (atom_i.charge*atom_j.charge*const.Q)/(4*pi*(r_ij**2))

        return e

    def tot_gradiant(self):
        """compute the self.enery of total gradiant"""

        self.e_bonds = self.bonds_gradiant()
        self.e_angles = self.angles_gradiant()
        self.e_torsions = self.torsions_gradiant()
        self.e_oofs = self.outofplanes_gradiant()
        self.e_bonded = self.e_bonds + self.e_angles + self.e_torsions + self.e_oofs
        self.e_cross = self.cross_terms()

        self.e_vdws = self.vdws_gradiant()
        self.e_elsts = self.elsts_gradiant()
        self.e_nonbonded = self.e_vdws + self.e_elsts

        self.e_total = self.e_bonds + self.e_angles + self.e_torsions + self.e_oofs + self.e_vdws + self.e_elsts

        return self.e_total

    def __repr__(self) -> str:
        return f'Gradiant({self.e_bonds}, {self.e_angles}, {self.e_torsions}, {self.e_vdws}, {self.e_elsts}, {self.e_total})'

