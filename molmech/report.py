
def print_report(mol):

    """
        print report of calcualtion done

        Args:
            mol (Molecule) : the molecule to print its calculations
        
        Retruns:
            None
    """

    # print initial paramters and positions
    print("Parsing Molecule Initial Parameters")
    print("=" * 79)
    print("{:<4}{:>10}{:>15}{:>10}{:>10}{:>10}{:>10}{:>10}".format("#", "Symbol", "x", "y", "z", "ro", 'eps', 'q'))
    print("=" * 79)
    for idx, a in enumerate(mol.atoms):
        print("{:<4}{:>10}{:>15}{:>10}{:>10}{:>10}{:>10}{:>10}".format(
            idx+1, a.symbol, round(a.pos.x, 4), round(a.pos.y, 4), round(a.pos.z, 4), 
            round(a.ro, 4), round(a.eps, 4), round(a.charge, 4)))
    print("=" * 79, end="\n\n")

    # print bonds
    print("Found {} Bond(s) in Molecule".format(mol.n_bonds))
    print("=" * 65)
    print("{:<10} {:>10} {:>10} {:>10} {:>10} {:>10}".format("Index", "Symbol", "r(Å)", "ro(Å)", "kb", 'energy'))
    print("=" * 65)
    for b in mol.bonds:
        print("{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} ".format(
            b.index, b.symbol, round(b.r(), 4), round(b.eq_r, 4), round(b.kb, 4), round(b.b_en, 4))
        )
    print("=" * 65, end="\n\n")

    # print angles
    print("Found {} Angles(s) in Molecule".format(mol.n_angles))
    print("=" * 65)
    print("{:<10} {:>10} {:>10} {:>10} {:>10} {:>10}".format("Index", "Symbol", "\N{greek capital letter Theta}°", "\N{greek capital letter theta}o°", "ka", "energy"))
    print("=" * 65)
    for a in mol.angles:
        print("{:<10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
            a.index, a.symbol, round(a.theta(), 0), round(a.eq_theta, 1), round(a.ka, 4), round(a.a_en, 4))
        )
    print("=" * 65, end="\n\n")

    # print torsiosn
    print("Found {} Torsion(s) in Molecule".format(mol.n_torsions))
    print("=" * 86)
    print("{:<15} {:>15} {:>10} {:>10} {:>10} {:>10} {:>10}".format("Index", "Symbol", "\N{greek capital letter Phi}°", "\N{greek capital letter Gamma}°", "Vn", 'n', "energy"))
    print("=" * 86)
    for t in mol.torsions:
        print("{:<15} {:>15} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
            t.index, t.symbol, round(t.dihedral(), 0), round(t.gamma, 1), round(t.vn, 4), round(t.nfolds, 4), round(t.t_en, 4) )
        )
    print("=" * 86, end="\n\n")

    # print out of planes
    print("Found {} Out of Plane Angle(s) in Molecule".format(mol.n_out_of_planes))
    print("=" * 64)
    print("{:<15} {:>15} {:>10} {:>10} {:>10}".format("Index", "Symbol", "\N{greek capital letter Phi}°", "Vn", "energy"))
    print("=" * 64)
    for o in mol.out_of_planes:
        print("{:<15} {:>15} {:>10} {:>10} {:>10}".format(o.index, o.symbol, round(o.phi(), 0), round(o.vn, 4), round(o.o_en, 4)))
    print("=" * 64, end="\n\n")

    # print center of mass
    print("Molecule Center of Mass")
    print("=" * 36)
    print("{:>8} {:>12} {:>12}".format("x", "y", "z"))
    print("=" * 36)
    print("{:<8} {:>12} {:>12}".format(round(mol.com.x, 4), round(mol.com.y, 4), round(mol.com.z, 4)))
    print("=" * 36, end="\n\n")

    # print moment of interia
    print("Molecule Moment of Inertia")
    print("=" * 43)
    print("{:>4} {:>8} {:>14} {:>14}".format("", "x", "y", "z"))
    print("=" * 43)
    labels = ["x", "y", "z"]
    for l, row in zip(labels, mol.moi.tensor):
        print("{:<4} {:>8} {:>14} {:>14}".format(l, round(row[0], 4), round(row[1], 4), round(row[2], 4)))
    print("=" * 43, end="\n\n")

    # print principle moment of inertia
    print("Molecule Prniciple Moment of Inertia")
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format("Ia", "Ib", "Ic"))
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format(round(mol.pmoi.Ia, 4), round(mol.pmoi.Ib, 4), round(mol.pmoi.Ic, 4)))
    print("=" * 38, end="\n\n")

    # print molecule type
    print("Molecule Type is : {}".format(mol.pmoi.mol_type))
    print()

    # print rotational consts in cm-1
    print("Molecule Rotational Consts (MHZ)")
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format("A", "B", "C"))
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format(round(mol.pmoi.A, 4), round(mol.pmoi.B, 4), round(mol.pmoi.C, 4)))
    print("=" * 38, end="\n\n")

    # print rotational consts in cm-1
    print("Molecule Rotational Consts (cm -1)")
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format("A", "B", "C"))
    print("=" * 38)
    print("{:<10} {:>12} {:>14}".format(round(mol.pmoi.Abar, 4), round(mol.pmoi.Bbar, 4), round(mol.pmoi.Cbar, 4)))
    print("=" * 38, end="\n\n")


def write_report(mol, filename):
    """
        write calcualtions done to a text file

        Args:
            mol (Molecule) : the molecule to print its calculations
            filename (str) : path to the output filename
        
        Retruns:
            None
    """

    with open(filename, "w") as f:
        f.write("Parsing Molecule Initial Parameters\n")
        f.write("=" * 53 + "\n")
        f.write("{:<10} {:>10} {:>15} {:>15}\n".format("Symbol", "x", "y", "z"))
        f.write("=" * 53 + "\n")
        for a in mol.atoms:
            f.write("{:<10} {:>10} {:>15} {:>15}\n".format(
                a.symbol, a.pos.x, a.pos.y, a.pos.z))
        f.write("=" * 53 + "\n\n")

        # f.write bonds
        f.write("Found {} Bond(s) in Molecule\n".format(mol.n_bonds))
        f.write("=" * 39 + "\n")
        f.write("{:<12} {:>8} {:>15} Å\n".format("Index", "Symbol", "R"))
        f.write("=" * 39 + "\n")
        for b in mol.bonds:
            f.write("{:<12} {:>8} {:>15} Å\n".format(b.index, b.symbol, b.r()))
        f.write("=" * 39 + "\n\n")

        # f.write angles
        f.write("Found {} Angles(s) in Molecule\n".format(mol.n_angles))
        f.write("=" * 39 + "\n")
        f.write("{:<12} {:>8} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter theta}"))
        f.write("=" * 39 + "\n")
        for a in mol.angles:
            f.write("{:<12} {:>8} {:>15} °\n".format(a.index, a.symbol, round(a.theta(), 4)))
        f.write("=" * 39 + "\n\n")

        # f.write torsiosn
        f.write("Found {} Torsion(s) in Molecule\n".format(mol.n_torsions))
        f.write("=" * 51 + "\n")
        f.write("{:<16} {:>16} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter Phi}"))
        f.write("=" * 51 + "\n")
        for t in mol.torsions:
            f.write("{:<16} {:>16} {:>15} °\n".format(t.index, t.symbol, round(t.dihedral(), 4)))
        f.write("=" * 51 + "\n\n")

        # print out of planes
        f.write("Found {} Out of Plane Angle(s) in Molecule\n".format(mol.n_out_of_planes))
        f.write("=" * 51 + "\n")
        f.write("{:<16} {:>16} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter Phi}"))
        f.write("=" * 51 + "\n")
        for o in mol.out_of_planes:
            f.write("{:<16} {:>16} {:>15} °\n".format(o.index, o.symbol, round(o.phi(), 0)))
        f.write("=" * 51 + "\n\n")

        # f.write center of mass
        f.write("Molecule Center of Mass\n")
        f.write("=" * 36 + "\n")
        f.write("{:>8} {:>12} {:>12}\n".format("x", "y", "z"))
        f.write("=" * 36 + "\n")
        f.write("{:<8} {:>12} {:>12}\n".format(round(mol.com.x, 4), round(mol.com.y, 4), round(mol.com.z, 4)))
        f.write("=" * 36 + "\n\n")

        # f.write moment of interia
        f.write("Molecule Moment of Inertia\n")
        f.write("=" * 43 + "\n")
        f.write("{:>4} {:>8} {:>14} {:>14}\n".format("", "x", "y", "z"))
        f.write("=" * 43 + "\n")
        labels = ["x", "y", "z"]
        for l, row in zip(labels, mol.moi.tensor):
            f.write("{:<4} {:>8} {:>14} {:>14}\n".format(l, round(row[0], 4), round(row[1], 4), round(row[2], 4)))
        f.write("=" * 43 + "\n\n")

        # f.write principle moment of inertia
        f.write("Molecule Prniciple Moment of Inertia\n")
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format("Ia", "Ib", "Ic"))
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.Ia, 4), round(mol.pmoi.Ib, 4), round(mol.pmoi.Ic, 4)))
        f.write("=" * 38 + "\n\n")

        # f.write molecule type
        f.write("Molecule Type is : {}\n".format(mol.pmoi.mol_type))
        f.write("\n")

        # f.write rotational consts in cm-1
        f.write("Molecule Rotational Consts (MHZ)\n")
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format("A", "B", "C"))
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.A, 4), round(mol.pmoi.B, 4), round(mol.pmoi.C, 4)))
        f.write("=" * 38 + "\n\n")

        # f.write rotational consts in cm-1
        f.write("Molecule Rotational Consts (cm -1)\n")
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format("A", "B", "C"))
        f.write("=" * 38 + "\n")
        f.write("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.Abar, 4), round(mol.pmoi.Bbar, 4), round(mol.pmoi.Cbar, 4)))
        f.write("=" * 38 + "\n\n")

def generate_report(mol) -> list:
    """
        generate text string of the report

        Args:
            mol (Molecule) : the molecule to generate the report for

        Returns:
            report (str) : the report in list format.
    """

    report = []

    report.append("Parsing Molecule Initial Parameters\n")
    report.append("=" * 53 + "\n")
    report.append("{:<10} {:>10} {:>15} {:>15}\n".format("Symbol", "x", "y", "z"))
    report.append("=" * 53 + "\n")
    for a in mol.atoms:
        report.append("{:<10} {:>10} {:>15} {:>15}\n".format(
            a.symbol, a.pos.x, a.pos.y, a.pos.z))
    report.append("=" * 53 + "\n\n")

    # f.write bonds
    report.append("Found {} Bond(s) in Molecule\n".format(mol.n_bonds))
    report.append("=" * 39 + "\n")
    report.append("{:<12} {:>8} {:>15} Å\n".format("Index", "Symbol", "R"))
    report.append("=" * 39 + "\n")
    for b in mol.bonds:
        report.append("{:<12} {:>8} {:>15} Å\n".format(b.index, b.symbol, b.r()))
    report.append("=" * 39 + "\n\n")

    # f.write angles
    report.append("Found {} Angles(s) in Molecule\n".format(mol.n_angles))
    report.append("=" * 39 + "\n")
    report.append("{:<12} {:>8} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter theta}"))
    report.append("=" * 39 + "\n")
    for a in mol.angles:
        report.append("{:<12} {:>8} {:>15} °\n".format(a.index, a.symbol, round(a.theta(), 4)))
    report.append("=" * 39 + "\n\n")

    # f.write torsiosn
    report.append("Found {} Torsion(s) in Molecule\n".format(mol.n_torsions))
    report.append("=" * 51 + "\n")
    report.append("{:<16} {:>16} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter Phi}"))
    report.append("=" * 51 + "\n")
    for t in mol.torsions:
        report.append("{:<16} {:>16} {:>15} °\n".format(t.index, t.symbol, round(t.dihedral(), 4)))
    report.append("=" * 51 + "\n\n")

    # print out of planes
    report.append("Found {} Out of Plane Angle(s) in Molecule\n".format(mol.n_out_of_planes))
    report.append("=" * 51 + "\n")
    report.append("{:<16} {:>16} {:>15} °\n".format("Index", "Symbol", "\N{greek capital letter Phi}"))
    report.append("=" * 51 + "\n")
    for o in mol.out_of_planes:
        report.append("{:<16} {:>16} {:>15} °\n".format(o.index, o.symbol, round(o.phi(), 0)))
    report.append("=" * 51 + "\n\n")

    # f.write center of mass
    report.append("Molecule Center of Mass\n")
    report.append("=" * 36 + "\n")
    report.append("{:>8} {:>12} {:>12}\n".format("x", "y", "z"))
    report.append("=" * 36 + "\n")
    report.append("{:<8} {:>12} {:>12}\n".format(round(mol.com.x, 4), round(mol.com.y, 4), round(mol.com.z, 4)))
    report.append("=" * 36 + "\n\n")

    # f.write moment of interia
    report.append("Molecule Moment of Inertia\n")
    report.append("=" * 43 + "\n")
    report.append("{:>4} {:>8} {:>14} {:>14}\n".format("", "x", "y", "z"))
    report.append("=" * 43 + "\n")
    labels = ["x", "y", "z"]
    for l, row in zip(labels, mol.moi.tensor):
        report.append("{:<4} {:>8} {:>14} {:>14}\n".format(l, round(row[0], 4), round(row[1], 4), round(row[2], 4)))
    report.append("=" * 43 + "\n\n")

    # f.write principle moment of inertia
    report.append("Molecule Prniciple Moment of Inertia\n")
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format("Ia", "Ib", "Ic"))
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.Ia, 4), round(mol.pmoi.Ib, 4), round(mol.pmoi.Ic, 4)))
    report.append("=" * 38 + "\n\n")

    # f.write molecule type
    report.append("Molecule Type is : {}\n".format(mol.pmoi.mol_type))
    report.append("\n")

    # f.write rotational consts in cm-1
    report.append("Molecule Rotational Consts (MHZ)\n")
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format("A", "B", "C"))
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.A, 4), round(mol.pmoi.B, 4), round(mol.pmoi.C, 4)))
    report.append("=" * 38 + "\n\n")

    # f.write rotational consts in cm-1
    report.append("Molecule Rotational Consts (cm -1)\n")
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format("A", "B", "C"))
    report.append("=" * 38 + "\n")
    report.append("{:<10} {:>12} {:>14}\n".format(round(mol.pmoi.Abar, 4), round(mol.pmoi.Bbar, 4), round(mol.pmoi.Cbar, 4)))
    report.append("=" * 38 + "\n\n")

    return report


def report_energy(mol):
    """
        print a table with different energy values
    """

    en = mol.get_energy()

    e_bs    = en.e_bonds
    e_as    = en.e_angles
    e_ts    = en.e_torsions
    e_os    = en.e_oofs
    e_bonded = en.e_bonded

    e_vdws  = en.e_vdws
    e_elst  = en.e_elsts
    e_nonbonded = en.e_nonbonded

    e_tot   = en.e_total


    au_e_bs    = en.e_bonds * 0.0016
    au_e_as    = en.e_angles * 0.0016
    au_e_ts    = en.e_torsions * 0.0016
    au_e_os    = en.e_oofs * 0.0016
    au_e_bonded = en.e_bonded * 0.0016

    au_e_vdws  = en.e_vdws * 0.0016
    au_e_elst  = en.e_elsts * 0.0016
    au_e_nonbonded = en.e_nonbonded * 0.0016

    au_e_tot   = en.e_total * 0.0016


    kj_e_bs    = en.e_bonds * 4.184
    kj_e_as    = en.e_angles * 4.184
    kj_e_ts    = en.e_torsions * 4.184
    kj_e_os    = en.e_oofs * 4.184
    kj_e_bonded = en.e_bonded * 4.184

    kj_e_vdws  = en.e_vdws * 4.184
    kj_e_elst  = en.e_elsts * 4.184
    kj_e_nonbonded = en.e_nonbonded * 4.184

    kj_e_tot   = en.e_total * 4.184

    print("Molecule Energies")
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Energy", "[Kcal/mol]", "[KJ/mol]", "a.u."))
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Bonds Energy", round(e_bs, 5), round(kj_e_bs, 5), round(au_e_bs, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Angles Energy", round(e_as, 5), round(kj_e_as, 5), round(au_e_as, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Torsions Energy", round(e_ts, 5), round(kj_e_ts, 5), round(au_e_ts, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("OutOfPlane Energy", round(e_os, 5), round(kj_e_os, 5), round(au_e_os, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Bonded Energy", round(e_bonded, 5), round(kj_e_bonded, 5), round(au_e_bonded, 5)))
    print("="*59)
    # print("{:025}{1>2{:>13}{:>13}0}".format("Cross Term Energy", round(e_cross, 5), round(kj_e_cross, 5), round(au_e_cross, 5)))
    # print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Vender Vals Energy", round(e_vdws, 5), round(kj_e_vdws, 5), round(au_e_vdws, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("ElectroStatic Energy", round(e_elst, 5), round(kj_e_elst, 5), round(au_e_elst, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Non Bonded Energy", round(e_nonbonded, 5), round(kj_e_nonbonded, 5), round(au_e_nonbonded, 5)))
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Total Energy", round(e_tot, 5), round(kj_e_tot, 5), round(au_e_tot, 5)))
    print("="*59, end="\n\n")


def report_gradiant(mol):
    """
        print a table with different gradiant values
    """

    en = mol.get_gradiant()

    e_bs    = en.e_bonds
    e_as    = en.e_angles
    e_ts    = en.e_torsions
    e_os    = en.e_oofs
    e_bonded = en.e_bonded

    e_vdws  = en.e_vdws
    e_elst  = en.e_elsts
    e_nonbonded = en.e_nonbonded

    e_tot   = en.e_total


    au_e_bs    = en.e_bonds * 0.0016
    au_e_as    = en.e_angles * 0.0016
    au_e_ts    = en.e_torsions * 0.0016
    au_e_os    = en.e_oofs * 0.0016
    au_e_bonded = en.e_bonded * 0.0016

    au_e_vdws  = en.e_vdws * 0.0016
    au_e_elst  = en.e_elsts * 0.0016
    au_e_nonbonded = en.e_nonbonded * 0.0016

    au_e_tot   = en.e_total * 0.0016


    kj_e_bs    = en.e_bonds * 4.184
    kj_e_as    = en.e_angles * 4.184
    kj_e_ts    = en.e_torsions * 4.184
    kj_e_os    = en.e_oofs * 4.184
    kj_e_bonded = en.e_bonded * 4.184

    kj_e_vdws  = en.e_vdws * 4.184
    kj_e_elst  = en.e_elsts * 4.184
    kj_e_nonbonded = en.e_nonbonded * 4.184

    kj_e_tot   = en.e_total * 4.184

    print("Molecule Gradiants")
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Gradiant", "[Kcal/mol]", "[KJ/mol]", "a.u."))
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Bonds Gradiant", round(e_bs, 5), round(kj_e_bs, 5), round(au_e_bs, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Angles Gradiant", round(e_as, 5), round(kj_e_as, 5), round(au_e_as, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Torsions Gradiant", round(e_ts, 5), round(kj_e_ts, 5), round(au_e_ts, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("OutOfPlane Gradiant", round(e_os, 5), round(kj_e_os, 5), round(au_e_os, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Bonded Gradiant", round(e_bonded, 5), round(kj_e_bonded, 5), round(au_e_bonded, 5)))
    print("="*59)
    # print("{:025}{1>2{:>13}{:>13}0}".format("Cross Term Energy", round(e_cross, 5), round(kj_e_cross, 5), round(au_e_cross, 5)))
    # print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Vender Vals Gradiant", round(e_vdws, 5), round(kj_e_vdws, 5), round(au_e_vdws, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("ElectStatic Gradiant", round(e_elst, 5), round(kj_e_elst, 5), round(au_e_elst, 5)))
    print("{:<20}{:>13}{:>13}{:>13}".format("Non Bonded Gradiant", round(e_nonbonded, 5), round(kj_e_nonbonded, 5), round(au_e_nonbonded, 5)))
    print("="*59)
    print("{:<20}{:>13}{:>13}{:>13}".format("Total Gradiant", round(e_tot, 5), round(kj_e_tot, 5), round(au_e_tot, 5)))
    print("="*59, end="\n\n")
