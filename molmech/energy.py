from math import cos, radians, sqrt, pi

from molmech.consts import const
from .mymath import euc

def get_bond_energy(bond):
    return bond.kb * (bond.r() - bond.eq_r)**2

def get_angle_energy(angle):
    return angle.ka * radians((angle.theta() - angle.eq_theta))**2

def get_torsion_energy(torsion):
    return 0.5 * torsion.vn * ( 
        1 + cos(torsion.nfolds*radians(torsion.dihedral()) - radians(torsion.gamma))
    )

def get_oof_energy(oof):
    return 0.5 * oof.vn * ( 1 + cos(radians(oof.phi()) - pi) )

class Energy:
    def __init__(self, atoms, bonds, angles, torsions, oofs) -> None:
        self.atoms          = atoms
        self.bonds          = bonds
        self.angles         = angles
        self.torsions       = torsions
        self.oofs            = oofs

        self.total_energy   = self.tot_energy()

    def bonds_energy(self):
        """compute the enery of bonds"""
        e = 0.0

        for b in self.bonds:
            e += get_bond_energy(b)

        return e

    def angles_energy(self):
        """compute the enery of angles"""
        e = 0.0

        for a in self.angles:
            e += get_angle_energy(a)

        return e

    def cross_terms(self):
        """compute the cross terms energy"""
        e = 0.0

        for b in self.bonds:
            for a in self.angles:
                e += a.ka * b.kb * (b.r() - b.eq_r) * (radians(a.theta() - a.eq_theta))

        return e

    def torsions_energy(self):
        """compute the enery of torsions"""
        e = 0.0

        for t in self.torsions:
            e += get_torsion_energy(t)

        return e

    def outofplanes_energy(self):
        """compute the enery of out of planes"""
        e = 0.0

        for o in self.oofs:
            e += get_oof_energy(o)

        return e

    def vdws_energy(self):
        """compute the enery of vdws"""
        e = 0.0

        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                atom_i = self.atoms[i]
                atom_j = self.atoms[j]

                eps_ij = sqrt(atom_i.eps * atom_j.eps)
                ro_ij = 0.5 * (atom_i.ro + atom_j.ro)
                r_ij = euc(atom_i.pos, atom_j.pos)

                e += eps_ij * ((ro_ij/r_ij)**12 - 2*(ro_ij/r_ij)**6)

        return e

    def elsts_energy(self):
        """compute the enery of electrostatic"""
        e = 0.0

        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                atom_i = self.atoms[i]
                atom_j = self.atoms[j]
                r_ij = euc(atom_i.pos, atom_j.pos)

                e += (atom_i.charge*atom_j.charge*const.Q)/(4*pi*r_ij)

        return e

    def tot_energy(self):
        """compute the self.enery of total energy"""
        
        self.e_bonds = self.bonds_energy()
        self.e_angles = self.angles_energy()
        self.e_torsions = self.torsions_energy()
        self.e_oofs = self.outofplanes_energy()
        self.e_bonded = self.e_bonds + self.e_angles + self.e_torsions + self.e_oofs
        self.e_cross = self.cross_terms()

        self.e_vdws = self.vdws_energy()
        self.e_elsts = self.elsts_energy()
        self.e_nonbonded = self.e_vdws + self.e_elsts
        
        self.e_total = self.e_bonds + self.e_angles + self.e_torsions + self.e_oofs + self.e_vdws + self.e_elsts

        return self.e_total
