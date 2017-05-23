from mdpy.entropy2 import get_atoms_id
import chemfiles
import os

class Data:

    def __init__(self, fname):
        self.n_Atoms = 0
        self.n_Bonds = 0
        self.n_Angles = 0
        self.n_Dihedrals = 0
        self.n_Impropers = 0
        self.n_AtomTypes = 0
        self.n_BondTypes = 0
        self.n_AngleTypes = 0
        self.n_DihedralTypes = 0
        self.n_ImproperTypes = 0
        self.xlo = 0
        self.xhi = 0
        self.ylo = 0
        self.yhi = 0
        self.zlo = 0
        self.zhi = 0

        self.Masses = []
        self.Bond_Coeffs = []
        self.Angle_Coeffs = []
        self.Dihedral_Coeffs = []
        self.Improper_Coeffs = []
        self.Atoms = []
        self.Velocities = []
        self.Bonds = []
        self.Angles = []
        self.Dihedrals = []
        self.Impropers = []

        with open(fname) as f:
            i_line = 0
            for line in f:

                i_line = i_line + 1

                if not line == "\n":
                    if line.strip().split()[-1] == 'atoms':
                        self.n_Atoms = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'bonds':
                        self.n_Bonds = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'angles':
                        self.n_Angles = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'dihedrals':
                        self.n_Dihedrals = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'impropers':
                        self.n_Impropers = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'types':
                        if line.strip().split()[1] == 'atom':
                            self.n_AtomTypes = int(line.strip().split()[0])
                        if line.strip().split()[1] == 'bond':
                            self.n_BondTypes = int(line.strip().split()[0])
                        if line.strip().split()[1] == 'angle':
                            self.n_AngleTypes = int(line.strip().split()[0])
                        if line.strip().split()[1] == 'dihedral':
                            self.n_DihedralTypes = int(line.strip().split()[0])
                        if line.strip().split()[1] == 'improper':
                            self.n_ImproperTypes = int(line.strip().split()[0])
                    if line.strip().split()[-1] == 'xhi':
                        self.xlo = float(line.strip().split()[0])
                        self.xhi = float(line.strip().split()[1])
                    if line.strip().split()[-1] == 'yhi':
                        self.ylo = float(line.strip().split()[0])
                        self.yhi = float(line.strip().split()[1])
                    if line.strip().split()[-1] == 'zhi':
                        self.zlo = float(line.strip().split()[0])
                        self.zhi = float(line.strip().split()[1])

                if line.strip() == "Masses":
                    self.i_line_Masses = i_line
                if line.strip() == "Bond Coeffs":
                    self.i_line_Bond_Coeffs = i_line
                if line.strip() == "Angle Coeffs":
                    self.i_line_Angle_Coeffs = i_line
                if line.strip() == "Dihedral Coeffs":
                    self.i_line_Dihedral_Coeffs = i_line
                if line.strip() == "Improper Coeffs":
                    self.i_line_Improper_Coeffs = i_line

                if line.strip() == "Atoms":
                    self.i_line_Atoms = i_line
                if line.strip() == "Bonds":
                    self.i_line_Bonds = i_line
                if line.strip() == "Angles":
                    self.i_line_Angles = i_line
                if line.strip() == "Dihedrals":
                    self.i_line_Dihedrals = i_line
                if line.strip() == "Impropers":
                    self.i_line_Impropers = i_line

        i_line = 0
        with open(fname) as f:
            for line in f:

                i_line = i_line + 1

                if (self.i_line_Masses+2) <= i_line <= (self.i_line_Masses+self.n_AtomTypes+1):
                    self.Masses.append(line)
                if (self.i_line_Bond_Coeffs+2) <= i_line <= (self.i_line_Bond_Coeffs+self.n_BondTypes+1):
                    self.Bond_Coeffs.append(line)
                if (self.i_line_Angle_Coeffs+2) <= i_line <= (self.i_line_Angle_Coeffs+self.n_AngleTypes+1):
                    self.Angle_Coeffs.append(line)
                if (self.i_line_Dihedral_Coeffs+2) <= i_line <= (self.i_line_Dihedral_Coeffs+self.n_DihedralTypes+1):
                    self.Dihedral_Coeffs.append(line)
                if (self.i_line_Improper_Coeffs+2) <= i_line <= (self.i_line_Improper_Coeffs+self.n_ImproperTypes+1):
                    self.Improper_Coeffs.append(line)

                if (self.i_line_Atoms+2) <= i_line <= (self.i_line_Atoms+self.n_Atoms+1):
                    self.Atoms.append(line)
                if (self.i_line_Bonds+2) <= i_line <= (self.i_line_Bonds+self.n_Bonds+1):
                    self.Bonds.append(line)
                if (self.i_line_Angles+2) <= i_line <= (self.i_line_Angles+self.n_Angles+1):
                    self.Angles.append(line)
                if (self.i_line_Dihedrals+2) <= i_line <= (self.i_line_Dihedrals+self.n_Dihedrals+1):
                    self.Dihedrals.append(line)
                if (self.i_line_Impropers+2) <= i_line <= (self.i_line_Impropers+self.n_Impropers+1):
                    self.Impropers.append(line)

    def __del__(self):
        class_name = self.__class__.__name__
        print(str(class_name) + " destroyed")
    
    def write(self, outname):
        with open(outname, 'w') as out:
            out.write('LAMMPS Description\n')
            out.write('\n')
            out.write('     %d  atoms\n' % self.n_Atoms)
            out.write('     %d  bonds\n' % self.n_Bonds)
            out.write('     %d  angles\n' % self.n_Angles)
            out.write('     %d  dihedrals\n' % self.n_Dihedrals)
            out.write('     %d  impropers\n' % self.n_Impropers)
            out.write('\n')
            out.write('     %d  atom types\n' % self.n_AtomTypes)
            out.write('     %d  bond types\n' % self.n_BondTypes)
            out.write('     %d  angle types\n' % self.n_AngleTypes)
            out.write('     %d  dihedral types\n' % self.n_DihedralTypes)
            out.write('     %d  improper types\n' % self.n_ImproperTypes)
            out.write('\n')
            out.write('  %.4f %.4f xlo xhi\n' % (self.xlo, self.xhi))
            out.write('  %.4f %.4f ylo yhi\n' % (self.ylo, self.yhi))
            out.write('  %.4f %.4f zlo zhi\n' % (self.zlo, self.zhi))
            out.write('\n')
            out.write('Masses\n')
            out.write('\n')
            for line in self.Masses:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Bond Coeffs\n')
            out.write('\n')
            for line in self.Bond_Coeffs:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Angle Coeffs\n')
            out.write('\n')
            for line in self.Angle_Coeffs:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Dihedral Coeffs\n')
            out.write('\n')
            for line in self.Dihedral_Coeffs:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Improper Coeffs\n')
            out.write('\n')
            for line in self.Improper_Coeffs:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Atoms\n')
            out.write('\n')
            for line in self.Atoms:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Bonds\n')
            out.write('\n')
            for line in self.Bonds:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Angles\n')
            out.write('\n')
            for line in self.Angles:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Dihedrals\n')
            out.write('\n')
            for line in self.Dihedrals:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')
            out.write('Impropers\n')
            out.write('\n')
            for line in self.Impropers:
                out.write('%s \n' % prettify_line(line))
            out.write('\n')

def prettify_line(line):
    return ' '.join(line.strip().split())

def replace_pos(f_data, f_traj, mol_types, out_data, translate_z=0):
    data = Data(f_data)

    traj = chemfiles.Trajectory(f_traj)
    frame = traj.read()

    mol_atomid = []
    for mol_type in mol_types:
        mol_atomid.extend(get_atoms_id(mol_type, frame))
    mol_atomid.sort() # array is 0-index

    positions = frame.positions()
    mol_pos = positions[mol_atomid]

    for loop, item in enumerate(mol_atomid):
        atom = data.Atoms[item].strip().split()
        pos = mol_pos[loop]
        pos[2] += translate_z
        data.Atoms[item] = '%s %s %s %s %.4f %.4f %.4f %s %s' % (atom[0], atom[1], atom[2], atom[3], pos[0], pos[1], pos[2], atom[7], atom[8])

    data.write(out_data)

def sort_dump(f_traj, f_out='sorted.lammpstrj', f_pz='sort_dump.pz'):

    print("Don't forget to activate python2!")

    with open(f_pz, 'w') as out:
        out.write('dump_file = dump("%s")\n' % f_traj)
        out.write('select_timestep = input("Choose timestep: ")\n')
        out.write('dump_file.tselect.one(select_timestep)\n')
        out.write('dump_file.sort()\n')
        out.write('dump_file.write("%s")\n' % f_out)

    os.system('python /Users/tonnamb/Documents/Cluster_Research/Tools/pizza_py/pizza-9Oct15/src/pizza.py -f %s' % f_pz)
