from Atom import Atom
class UnitCell:
    def __init__(self, filename):
        self.ident = filename
        self.file2unit_cell()
        
    def __str__(self):
        return "<Unit cell: {}>".format(self.ident)

    def file2unit_cell(self):
        with open(self.ident) as f:
            lines = f.readlines()
        atom_read = False
        self.atom_array = []
        mass_read = False
        self.mass_dict = {}
        zs = []
        for line in lines:
            split = line.split()
            if len(split) >= 1:
                
                if split[0] == 'LAMMPS':
                    continue
                if split[-1] == 'atoms':
                    self.nAtoms = split[0]
                    continue
                if split[-1] == 'xhi':
                    self.dx = split[1]
                    continue
                if split[-1] == 'yhi':
                    self.dy = split[1]
                    continue
                if split[-1] == 'zhi':
                    self.dz = split[1]
                    continue
                if split[-1] == 'yz':
                    self.dxy = split[0]
                    continue
                if split[-1] == 'Masses':
                    mass_read = True
                    continue
                if split[-1] == "Atoms":
                    mass_read = False
                    atom_read = True
                    continue
                if mass_read == True:
                    self.mass_dict[split[0]] = split[1]
                if atom_read == True:
                    a = Atom(int(split[0]),int(split[1]),(float(split[2]),float(split[3]),float(split[4])))
                    zs.append(float(split[4]))
                    self.atom_array += [a]
        self.z_max = max(zs)
        self.z_min = min(zs)
        return 'Unit Cell Created'