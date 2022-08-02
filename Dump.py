import numpy as np
from Atom import Atom

class Dump:
    def __init__(self,filename):
        self.filename = filename
        split = filename.split(".")
        self.descriptor = split[1]
#         self.timestep = int(split[2])
        
        self.basis = 5 # for Ti2
        
        self.parse()
        self.nAtoms = len(self.atom_array)
        
        [self.dx, self.dy, self.average_dz] = self.get_dimensions()
        self.area = self.get_area()
        self.simulation_box_volume = self.get_simulation_box_volume()
        self.approximate_volume = self.get_approximate_volume()
        
    def __str__(self):
        return "<Dump object from: {}>".format(self.filename)
    
    def parse(self):
        with open(self.filename) as f:
            lines = f.readlines()
        atom_read = False
        self.atom_array = []
        self.xs = []
        self.ys = []
        self.zs = []
        for line in lines:
            split = line.split()
            if len(split) > 1:
                if line.split()[1] == "ATOMS":
                    atom_read = True
                    continue
            if atom_read == True:
                a = Atom(int(split[0]),int(split[1]),(float(split[2]),float(split[3]),float(split[4])))
                self.atom_array += [a]
                x = float(split[2])
                y = float(split[3])
                z = float(split[4])
                self.xs.append(x)
                self.ys.append(y)
                self.zs.append(z)
                
        return self.atom_array
    
    def get_dimensions(self):
        self.zs.sort()
        top = np.average(self.zs[-1*self.nAtoms//self.basis:])
        bottom = np.average(self.zs[0:self.nAtoms//self.basis])
        self.average_dz = top - bottom
        
        self.dy = max(self.ys) - min(self.ys)
        
        self.dx = max(self.xs)
        
        return [self.dx, self.dy, self.average_dz]
    
    def get_area(self):
        return self.dx * self.dy
    
    
    def get_simulation_box_volume(self):
        with open(self.filename) as f:
            lines = f.readlines()

        a1 = float(lines[6-1].split()[1])
        b1 = float(lines[6-1].split()[0])
        b2 = float(lines[7-1].split()[1])
        c3 = float(lines[8-1].split()[1])
        a = [a1, 0, 0]
        b = [b1, b2, 0]
        c = [0, 0, c3]

        self.simulation_box_volume = np.linalg.det(np.dstack([a,b,c]))[0]
        return self.simulation_box_volume
    
    def get_approximate_volume(self):
        return self.dy*self.dx*self.average_dz
    
    
    
        