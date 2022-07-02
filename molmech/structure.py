import numpy as np
from molmech.energy import (
    get_angle_energy,
    get_bond_energy,
    get_torsion_energy,
    get_oof_energy,
)

from molmech.mymath import (
    euc,
    compute_angle,
    compute_dihedral_angle,
    compute_out_of_plane_angle,
)

from .utils import elements
from .consts import const


class Position:

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def pos(self):
        return self.x, self.y, self.z

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"


class Atom:

    def __init__(self,
                 index: int,
                 symbol: int,
                 pos: Position,
                 _type=None,
                 charge=0.0,
                 ro=0.0,
                 eps=0.0) -> None:
        self.index = index
        self.symbol = symbol
        self.pos = pos

        self._type = _type
        self.charge = charge
        self.ro = ro
        self.eps = eps

    @property
    def atomic_number(self) -> int:
        return elements.get(self.symbol, None)

    def obj(self):
        return {
            "index": self.index,
            "symbol": self.symbol,
            "string": "{} {}".format(self.symbol, self.index),
            "pos": (self.pos.x, self.pos.y, self.pos.z),
            "type": self._type,
            "charge": self.charge,
            "ro": self.ro,
            "eps": self.eps
        }

    def __str__(self) -> str:
        return f"{self.symbol}{self.index}"

    def __repr__(self) -> str:
        return f"{self.symbol}"


class Bond:

    def __init__(self, atom_i: Atom, atom_j: Atom, kb=0.0, eq_r=0.0) -> None:
        self.i_atom = atom_i
        self.j_atom = atom_j
        self.kb = kb
        self.eq_r = eq_r

    def r(self) -> float:
        return round(euc(self.i_atom.pos, self.j_atom.pos), 4)

    @property
    def eq_r(self):
        return self._eq_r

    @eq_r.setter
    def eq_r(self, value):
        self._eq_r = value

    @property
    def kb(self):
        return self._kb

    @kb.setter
    def kb(self, value):
        self._kb = value

    @property
    def b_en(self):
        return get_bond_energy(self)

    @property
    def index(self) -> str:
        return f"({self.i_atom.index}-{self.j_atom.index})"

    @property
    def symbol(self) -> str:
        return f"({self.i_atom.symbol}-{self.j_atom.symbol})"

    def obj(self):
        return {
            "index": (self.i_atom.index, self.j_atom.index),
            "string": "{}-{}".format(str(self.i_atom), str(self.j_atom)),
            "kb": self.kb,
            "eq_r": self.eq_r,
            "en": self.b_en
        }

    def __str__(self) -> str:
        return (f"({self.i_atom.index}-{self.j_atom.index}) \t"
                f"({self.i_atom.symbol}-{self.j_atom.symbol})")


class Angle:

    def __init__(self,
                 atom_i: Atom,
                 atom_j: Atom,
                 atom_k: Atom,
                 ka=0.0,
                 eq_theta=0.0) -> None:
        self.i_atom = atom_i
        self.j_atom = atom_j
        self.k_atom = atom_k

        self.ka = ka
        self.eq_theta = eq_theta

    def theta(self) -> float:
        return compute_angle(self.i_atom.pos, self.j_atom.pos, self.k_atom.pos)

    @property
    def eq_theta(self):
        return self._eq_theta

    @eq_theta.setter
    def eq_theta(self, value):
        self._eq_theta = value

    @property
    def ka(self):
        return self._ka

    @ka.setter
    def ka(self, value):
        self._ka = value

    @property
    def a_en(self):
        return get_angle_energy(self)

    @property
    def index(self) -> str:
        return f"({self.i_atom.index}-{self.j_atom.index}-{self.k_atom.index})"

    @property
    def symbol(self) -> str:
        return f"({self.i_atom.symbol}-{self.j_atom.symbol}-{self.k_atom.symbol})"

    def obj(self):
        return {
            "index": (self.i_atom.index, self.j_atom.index, self.k_atom.index),
            "string": "{}-{}-{}".format(str(self.i_atom), str(self.j_atom), str(self.k_atom)),
            "ka": self.ka,
            "eq_theta": self.eq_theta,
            "en": self.a_en
        }

    def __str__(self) -> str:
        return (
            f"({self.i_atom.index}-{self.j_atom.index}-{self.k_atom.index}) \t"
            f"({self.i_atom.symbol}-{self.j_atom.symbol}-{self.k_atom.symbol})"
        )


class Torsion:

    def __init__(self,
                 atom_i: Atom,
                 atom_j: Atom,
                 atom_k: Atom,
                 atom_l: Atom,
                 vn=None,
                 gamma=None,
                 nfolds=None) -> None:
        self.i_atom = atom_i
        self.j_atom = atom_j
        self.k_atom = atom_k
        self.l_atom = atom_l

        self.vn = vn
        self.gamma = gamma
        self.nfolds = nfolds

    def dihedral(self) -> float:
        return compute_dihedral_angle(self.i_atom.pos, self.j_atom.pos,
                                      self.k_atom.pos, self.l_atom.pos)

    @property
    def vn(self):
        return self._vn

    @vn.setter
    def vn(self, value):
        self._vn = value

    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    @property
    def nfolds(self):
        return self._nfolds

    @nfolds.setter
    def nfolds(self, value):
        self._nfolds = value

    @property
    def t_en(self):
        return get_torsion_energy(self)

    @property
    def index(self) -> str:
        return f"({self.i_atom.index}-{self.j_atom.index}-{self.k_atom.index}-{self.l_atom.index})"

    @property
    def symbol(self) -> str:
        return f"({self.i_atom.symbol}-{self.j_atom.symbol}-{self.k_atom.symbol}-{self.l_atom.symbol})"

    def obj(self):
        return {
            "index": (self.i_atom.index, self.j_atom.index, self.k_atom.index, self.l_atom.index),
            "string": "{}-{}-{}-{}".format(str(self.i_atom), str(self.j_atom), str(self.k_atom), str(self.l_atom)),
            "vn": self.vn,
            "gamma": self.gamma,
            "n": self.nfolds,
            "en": self.t_en
        }

    def __str__(self) -> str:
        return (
            f"({self.i_atom.index}-{self.j_atom.index}-{self.k_atom.index}-{self.l_atom.index}) \t"
            f"({self.i_atom.symbol}-{self.j_atom.symbol}-{self.k_atom.symbol}-{self.l_atom.symbol})"
        )


class OutOfPlane:

    def __init__(self, atom_i, atom_j, atom_k, atom_l, vn=0.0) -> None:
        self.i_atom = atom_i
        self.j_atom = atom_j
        self.k_atom = atom_k
        self.l_atom = atom_l

        self.vn = vn

    def phi(self) -> float:
        return compute_out_of_plane_angle(self.i_atom.pos, self.j_atom.pos,
                                          self.k_atom.pos, self.l_atom.pos)

    @property
    def o_en(self):
        return get_oof_energy(self)

    @property
    def index(self) -> str:
        return f"({self.i_atom.index}-{self.j_atom.index}-{self.k_atom.index}-{self.l_atom.index})"

    @property
    def symbol(self) -> str:
        return f"({self.i_atom.symbol}-{self.j_atom.symbol}-{self.k_atom.symbol}-{self.l_atom.symbol})"


class COM:

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __str__(self) -> str:
        return f"[{self.x, self.y, self.z}]"


class MOI:

    def __init__(self, xx: float, yy: float, zz: float, xy, xz, yz) -> None:
        self.xx = xx
        self.yy = yy
        self.zz = zz

        self.xy, self.yx = xy, xy
        self.xz, self.zx = xz, xz
        self.yz, self.zy = yz, yz

    @property
    def tensor(self) -> list:
        _moi = [
            [self.xx, self.xy, self.xz],
            [self.xy, self.yy, self.xz],
            [self.xz, self.yz, self.zz],
        ]

        return _moi

    def __str__(self) -> str:
        return f"[{self.xx, self.yy, self.zz}]"


class PMOI:

    def __init__(self, moi_tensor: np.ndarray) -> None:
        self.moi_tensor = np.array(moi_tensor)

        self.solve()

    def solve(self):
        v = np.linalg.eigvalsh(self.moi_tensor)

        self.Ia, self.Ib, self.Ic = v

        self.Ia = np.real(self.Ia)
        self.Ib = np.real(self.Ib)
        self.Ic = np.real(self.Ic)

        self.A = self.compute_rot_const(self.Ia)
        self.B = self.compute_rot_const(self.Ib)
        self.C = self.compute_rot_const(self.Ic)

        self.Abar = self.compute_rot_const_number(self.Ia)
        self.Bbar = self.compute_rot_const_number(self.Ib)
        self.Cbar = self.compute_rot_const_number(self.Ic)

    def compute_rot_const(self, I):
        try:
            _s = (10**(-6) * const.h * const.na * 10**3 *
                  10**20) / (8 * np.pi**2 * I)
        except Exception as e:
            _s = 0.0
        return _s

    def compute_rot_const_number(self, I):
        return (const.h * const.na * 10**3 * 10**20) / (8 * np.pi**2 * I *
                                                        const.c)

    @property
    def mol_type(self):
        _type = ""

        if round(self.A, 0) > round(self.B, 0) > round(self.C, 0) > 0:
            _type = "[Asymmetric Top]"
        elif round(self.A, 0) == round(self.B, 0) == round(self.C, 0) == 0:
            _type = "[Mono Atomic]"
        elif round(self.A, 0) == round(self.B, 0) > round(self.C, 0) == 0:
            _type = "[Linear]"
        elif round(self.A, 0) == round(self.B, 0) > round(self.C, 0) > 0:
            _type = "[Oblate Symmetric Top]"
        elif round(self.A, 0) > round(self.B, 0) == round(self.C, 0) > 0:
            _type = "[Prolate Symmetric Top]"
        elif round(self.A, 0) == round(self.B, 0) == round(self.C, 0) > 0:
            _type = "[Spherical Top]"

        return _type
