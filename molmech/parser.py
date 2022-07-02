import json
import re

from . import force_fields as params

def parse_xyz(filename) -> list:
    """
        parse xyz files and convert to list of atoms
    """
    xyz = []

    multi_space = re.compile(r"\s+")

    with open(filename, 'r') as f:
        xyz_data = f.readlines()

    for l in  xyz_data[2:]:
        row = multi_space.sub(",", l.strip()).split(',')
        if len(row) != 4: continue
        row[1] = round(float(row[1]), 4)
        row[2] = round(float(row[2]), 4)
        row[3] = round(float(row[3]), 4)

        xyz.append(row)

    return xyz


def parse_mxyz(filename) -> list:
    """
        parse xyz files and convert to list of atoms
    """
    xyz = []

    multi_space = re.compile(r"\s+")

    with open(filename, 'r') as f:
        xyz_data = f.readlines()

    for l in  xyz_data[2:]:
        row = multi_space.sub(",", l.strip()).split(',')
        if len(row) != 5: continue
        row[1] = float(row[1])
        row[2] = float(row[2])
        row[3] = float(row[3])
        row[4] = float(row[4])

        xyz.append(row)

        _type = row[0]
        
        if len(_type) == 1:
            symbol = _type.upper()
        elif len(_type) == 2:
            if _type.lower().capitalize()[:-1] in params.elements: 
                symbol = _type.lower().capitalize()[:-1]
        elif len(_type) == 3:
            if _type.lower().capitalize()[:-2] in params.elements: 
                symbol = _type.lower().capitalize()[:-2]
        else:
            raise ValueError("elemtn with {} not exist".format(_type))

        row.append(symbol)

    return xyz

def convert_to_json(mol):
    """
        create an output file of json type of molecule

        Args:
            mol (Molecule) : the molecule to export it's structure
            filename (str) : the json filename
        Returns:
            None
    """

    output = dict()

    output['atoms'] = [atom.obj() for atom in mol.atoms]
    output['bonds'] = [bond.obj() for bond in mol.bonds]
    output['angles'] = [angle.obj() for angle in mol.angles]
    output['torsions'] = [torsion.obj() for torsion in mol.torsions]

    return output


def write_json(mol, filename):
    output = convert_to_json(mol)

    with open(filename, 'w') as f:
        json.dump(output, f)

    
