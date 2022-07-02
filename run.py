from argparse import ArgumentParser

from molmech import NAME, USAGE, Molecule
from molmech.parser import write_json
from molmech.report import report_energy, write_report, report_gradiant

parser = ArgumentParser(NAME, usage=USAGE)
parser.add_argument('--mxyz', '--minput', type=str, help="path to the mxyz file")
parser.add_argument('--xyz', '--input', type=str, help="path to the xyz file")
parser.add_argument('--out', '--output', type=str, help="path to the output file")

def main():
    args = parser.parse_args()

    if args.mxyz:
        mol = Molecule.from_mxyz(args.mxyz)
    else:
        mol = Molecule.from_xyz(args.xyz)
    mol()

    mol.print_report()
    report_energy(mol)
    # report_gradiant(mol)

    if args.out:
        write_report(mol, args.out)
        write_json(mol, args.out.split(".")[0]+'.json')

if __name__ == "__main__":
    main()
