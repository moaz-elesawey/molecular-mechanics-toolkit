from molmech.energy import Energy
from molmech.gradiant import Gradiant
from molmech.report import generate_report, print_report, write_report
from molmech.structure import COM, MOI, PMOI, Angle, Atom, Bond, OutOfPlane, Position, Torsion
from molmech.parser import convert_to_json, parse_xyz, parse_mxyz
from molmech.utils import cov_rads, elements
from molmech.mymath import euc
from molmech import force_fields as params

from time import time


class Molecule:
    def __init__(self, atoms: list) -> None:
        self.atoms = atoms
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.out_of_planes = []

        self.bonds_tree = {}

        self.n_atoms = len(atoms)
        self.n_bonds = 0
        self.n_angles = 0
        self.n_torsions = 0
        self.n_out_of_planes = 0
        self.energy = None
        self.gradiant = None

    @classmethod
    def from_xyz(cls, filename: str):
        _atoms = parse_xyz(filename)
        created_atoms = cls.create_atoms(_atoms)

        return cls(created_atoms)

    @classmethod
    def from_mxyz(cls, filename):
        _atoms = parse_mxyz(filename)
        created_atoms = cls.create_atoms_from_mxyz(_atoms)

        return cls(created_atoms)

    @staticmethod
    def create_atoms(atoms: list) -> list:
        created_atoms = []

        for i, a in enumerate(atoms):
            _pos = Position(a[1], a[2], a[3])
            atom = Atom(i+1, a[0], _pos)

            created_atoms.append(atom)

        return created_atoms

    @staticmethod
    def create_atoms_from_mxyz(atoms: list):
        created_atoms = []

        for i, a in enumerate(atoms):
            _pos = Position(a[1], a[2], a[3])
            _type = a[0]
            vdw = params.get_vdw_param(_type)
            atom = Atom(index=i+1, symbol=a[-1], pos=_pos, _type=_type, charge=a[4], ro=vdw[0], eps=vdw[1])

            created_atoms.append(atom)

        return created_atoms

    # set bonds
    def set_bonds(self):
        t = time()

        for i in range(self.n_atoms):
            for j in range(i+1, self.n_atoms):

                atom_i: Atom = self.atoms[i]
                atom_j: Atom = self.atoms[j]

                r_ij = euc(atom_i.pos, atom_j.pos)
                cov_ij = 1.2 * (cov_rads.get(atom_i.symbol) +
                                cov_rads.get(atom_j.symbol))

                if (r_ij < cov_ij):
                    bond_param = params.get_bond_param(atom_i._type, atom_j._type)

                    bond = Bond(atom_i, atom_j, kb=bond_param[0], eq_r=bond_param[1])
                    self.bonds.append(bond)

                    if atom_i in self.bonds_tree:
                        self.bonds_tree[atom_i].append(atom_j)
                    else:
                        self.bonds_tree[atom_i] = [atom_j]

                    if atom_j in self.bonds_tree:
                        self.bonds_tree[atom_j].append(atom_i)
                    else:
                        self.bonds_tree[atom_j] = [atom_i]

        self.n_bonds = len(self.bonds)

    def set_angles(self):
        for center_atom in self.bonds_tree:
            center_tree = self.bonds_tree[center_atom]

            for i in range(len(center_tree)):
                for j in range(i+1, len(center_tree)):
                    atom_i = center_tree[i]
                    atom_j = center_tree[j]
                    angle_param = params.get_angle_param(atom_i._type, center_atom._type, atom_j._type)

                    angle = Angle(atom_i, center_atom, atom_j, ka=angle_param[0], eq_theta=angle_param[1])
                    self.angles.append(angle)
        self.n_angles = len(self.angles)

    def set_torsions(self):
        for origin_i in self.bonds_tree:
            bonded_i = self.bonds_tree[origin_i]

            for origin_j in bonded_i:
                bonded_j = self.bonds_tree[origin_j]
                if origin_i == origin_j: continue

                for origin_k in bonded_j:
                    bonded_k = self.bonds_tree[origin_k]
                    if origin_j == origin_k or origin_i == origin_k: continue

                    for origin_l in bonded_k:
                        if origin_i == origin_l or origin_j == origin_l or origin_k == origin_l: continue
                        torsion_param = params.get_torision_param(origin_i._type,
                                origin_j._type, origin_k._type, origin_l._type)

                        torsion = Torsion(origin_i, origin_j, origin_k, origin_l, vn=torsion_param[0],
                                gamma=torsion_param[1], nfolds=torsion_param[-1])
                        self.torsions.append(torsion)

        self.n_torsions = len(self.torsions)//2
        self.torsions = self.torsions[:self.n_torsions]

    def set_out_of_planes(self):
        for origin, tree in self.bonds_tree.items():
            if len(tree) == 3:
                a0, a1, a2 = tree
                oof_param1 = params.get_oof_param(a0._type, a1._type, a2._type, origin._type)
                oof_param2 = params.get_oof_param(a2._type, a0._type, a1._type, origin._type)
                oof_param3 = params.get_oof_param(a1._type, a0._type, a2._type, origin._type)

                out1 = OutOfPlane(a0, a1, a2, origin, vn=oof_param1)
                out2 = OutOfPlane(a2, a0, a1, origin, vn=oof_param2)
                out3 = OutOfPlane(a1, a0, a2, origin, vn=oof_param3)

                for out in (out1, out2, out3):
                    self.out_of_planes.append(out)

        self.n_out_of_planes = len(self.out_of_planes)

    def set_center_of_mass(self):
        x_cm, y_cm, z_cm = 0, 0, 0
        mass = 0

        for atom in self.atoms:
            atomic_mass = elements.get(atom.symbol)[1]

            x_cm += atomic_mass * atom.pos.x
            y_cm += atomic_mass * atom.pos.y
            z_cm += atomic_mass * atom.pos.z

            mass += atomic_mass

        x_cm /= mass
        y_cm /= mass
        z_cm /= mass

        self.com = COM(x_cm, y_cm, z_cm)

    def set_moment_of_inertia(self):
        xx_cm, yy_cm, zz_cm = 0, 0, 0
        xy_cm, yz_cm, xz_cm = 0, 0, 0

        for atom in self.atoms:
            atomic_mass = elements.get(atom.symbol)[1]

            xx_cm += atomic_mass * ((atom.pos.y - self.com.y)**2 + (atom.pos.z - self.com.z)**2)
            yy_cm += atomic_mass * ((atom.pos.x - self.com.x)**2 + (atom.pos.z - self.com.z)**2)
            zz_cm += atomic_mass * ((atom.pos.x - self.com.x)**2 + (atom.pos.y - self.com.y)**2)

            xy_cm += -atomic_mass * ((atom.pos.x - self.com.x) * (atom.pos.y - self.com.y))
            xz_cm += -atomic_mass * ((atom.pos.x - self.com.x) * (atom.pos.z - self.com.z))
            yz_cm += -atomic_mass * ((atom.pos.y - self.com.y) * (atom.pos.z - self.com.z))

        self.moi = MOI(xx_cm, yy_cm, zz_cm, xy_cm, xz_cm, yz_cm)

        self.pmoi = PMOI(self.moi.tensor)

    @property
    def json(self):
        return convert_to_json(self)

    def set_energy(self):
        self.energy = Energy(self.atoms, self.bonds, self.angles, self.torsions, self.out_of_planes)

    def get_energy(self):
        if self.energy is not None:
            return self.energy

    def set_gradiant(self):
        self.gradiant = Gradiant(self.atoms, self.bonds, self.angles, self.torsions, self.out_of_planes)

    def get_gradiant(self):
        if self.gradiant is not None:
            return self.gradiant

    def __call__(self):
        start = time()
        self.set_bonds()
        self.set_angles()
        self.set_torsions()
        self.set_out_of_planes()
        end = time()
        # print("Took {:8.5f} to set up the Molecule".format(end-start))
        self.set_center_of_mass()
        self.set_moment_of_inertia()
        self.set_energy()
        self.set_gradiant()

    def set_energy(self):
        self.energy = Energy(self.atoms, self.bonds, self.angles, self.torsions, self.out_of_planes)

    def get_energy(self):
        if self.energy is not None:
            return self.energy

    def write(self):
        print(self.gradiant)
        write_report(self)

    def print_report(self):
        print_report(self)

    def get_report(self):
        return generate_report(self)

    def print_bond_tree(self):
        for b in self.bonds_tree:
            print("{:<6} => {:<8}".format(str(b), str(self.bonds_tree[b])))
