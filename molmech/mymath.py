import importlib
from math import radians, sqrt, acos, degrees, sin, asin

def euc(pos_i, pos_j) -> float:
    """
        compute the bond length between two atoms using euc method
    """
    d = sqrt(((pos_j.x - pos_i.x)**2) + ((pos_j.y - pos_i.y)**2) + ((pos_j.z - pos_i.z)**2))

    return d


def convert_to_unit_vectors(pos_i, pos_j):
    """
        convert a vector i to unit vector with respect to vector j.
    """
    Position = importlib.import_module("molmech.structure").Position

    r_ji = euc(pos_i, pos_j)
    unit = Position((pos_j.x - pos_i.x) / r_ji, (pos_j.y - pos_i.y) / r_ji, (pos_j.z - pos_i.z) / r_ji)

    return unit


def compute_dot_product(pos_i, pos_j):
    """
        calculate the dot product of atom vectors
    """

    dot = (pos_i.x * pos_j.x) + (pos_i.y * pos_j.y) + (pos_i.z * pos_j.z)

    return dot

def compute_cross_product(pos_i, pos_j):
    """
        calculate the cross product for atom vectors
    """

    Position = importlib.import_module("molmech.structure").Position

    x = ((pos_i.y * pos_j.z) - (pos_i.z * pos_j.y))
    y = ((pos_i.x * pos_j.z) - (pos_i.z * pos_j.x))
    z = ((pos_i.y * pos_j.x) - (pos_i.x * pos_j.y))

    cross = Position(x, y, z)

    return cross

def compute_angle(pos_i, pos_j, pos_k) -> float:
    """
        calculate the angles between bonds
    """
    u_ji = convert_to_unit_vectors(pos_i, pos_j)
    u_jk = convert_to_unit_vectors(pos_k, pos_j)

    dot = compute_dot_product(u_ji, u_jk)

    angle = degrees(acos(dot))

    return angle


def compute_dihedral_angle(pos_i, pos_j, pos_k, pos_l):
    """
        calculate the the torsion (dihedral) angles in the molecule
    """

    dihedral = 0.0

    u_ji = convert_to_unit_vectors(pos_i, pos_j)
    u_jk = convert_to_unit_vectors(pos_k, pos_j)

    u_kj = convert_to_unit_vectors(pos_j, pos_k)
    u_kl = convert_to_unit_vectors(pos_l, pos_k)

    pos_ijk = compute_cross_product(u_ji, u_jk)
    pos_jkl = compute_cross_product(u_kj, u_kl)

    dot  = compute_dot_product(pos_ijk, pos_jkl)

    angle_ijk_rad = radians(compute_angle(pos_i, pos_j, pos_k))
    angle_jkl_rad = radians(compute_angle(pos_j, pos_k, pos_l))

    try:
        dihedral = degrees(
            acos(dot / ( sin(angle_ijk_rad) * sin(angle_jkl_rad)))
        )
    except Exception as e:
        pass

    return dihedral


def compute_out_of_plane_angle(pos_i, pos_j, pos_k, pos_l) -> float:
    """
        calculate the out of plane angles between atoms i, j, k and l
    """

    phi = 0.0

    u_lj = convert_to_unit_vectors(pos_j, pos_l)
    u_lk = convert_to_unit_vectors(pos_k, pos_l)
    u_li = convert_to_unit_vectors(pos_i, pos_l)

    c_ik = compute_cross_product(u_lj, u_lk)
    theta_ilk = sin(radians(compute_angle(pos_j, pos_l, pos_k)))
    
    c_ik.x /= theta_ilk
    c_ik.y /= theta_ilk
    c_ik.z /= theta_ilk

    dot = compute_dot_product(c_ik, u_li)

    try:
        phi = degrees(asin(dot))
    except Exception as e:
        print(e)

    return phi
