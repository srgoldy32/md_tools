import numpy as np
import itertools
from multiprocessing import Pool
from Atom import Atom

class Dump:
    def __init__(self,filename):
        self.filename = filename
        split = filename.split(".")
        self.descriptor = split[1]
#         self.timestep = int(split[2])
        
#         self.basis = 5 # for Ti2
        
        self.parse()
#         self.nAtoms = len(self.atom_array)
        
#         [self.dx, self.dy, self.average_dz] = self.get_dimensions()
#         self.area = self.get_area()
#         self.simulation_box_volume = self.get_simulation_box_volume()
#         self.approximate_volume = self.get_approximate_volume()
        
    def __str__(self):
        return "<Dump object from: {}>".format(self.filename)
    def parse(self):
        with open(self.filename) as f:
            lines = f.readlines()
        self.beg = lines[:8]
        self.atom_lines = lines[8:]
        self.get_box_dims()
        self.get_periodicity()
        self.atom_dict = {}
        
    
    def read_atoms(self):
        top = self.atom_lines.pop(0)
        split = top.split()
        if len(split) > 1:
            if split[1] == "ATOMS":
                extra_keys = []
                if len(split) > 7:
                    extra_keys = split[7:]
        for line in self.atom_lines:
            split = line.split()
            x = float(split[2])
            y = float(split[3])
            z = float(split[4])
            extra_values = []
            if len(split) > 5:
                extra_values = [float(i) for i in split[5:]]


            extra_dict = dict(zip(extra_keys, extra_values))
            a = Atom(int(split[0]),int(split[1]),(x,y,z),extra_dict)
            self.atom_dict[int(split[0])] = a
        
        
        
    def parse_old(self):
        with open(self.filename) as f:
            lines = f.readlines()
        atom_read = False
        self.atom_array = []
        self.periodicity = [True,True,True]
        self.get_box_dims(lines)
        
        self.xs = []
        self.ys = []
        self.zs = []
        for line in lines:
            split = line.split()
            if len(split) > 1:
                if split[1] == "ATOMS":
                    extra_keys = []
                    if len(split) > 7:
                        extra_keys = split[7:]
                        
                    atom_read = True
                    continue
                if split[1] == "BOX":
                    for i in [6,7,8]:
                        if split[i] != 'pp':
                            self.periodicity[i-6] = False
            if atom_read == True:
                x = float(split[2])
                y = float(split[3])
                z = float(split[4])
                self.xs.append(x)
                self.ys.append(y)
                self.zs.append(z)
                extra_values = []
                if len(split) > 5:
                    extra_values = [float(i) for i in split[5:]]
                
                
                extra_dict = dict(zip(extra_keys, extra_values))
                a = Atom(int(split[0]),int(split[1]),(x,y,z),extra_dict)
                self.atom_array += [a]
                
        return self.atom_array
        
    def get_box_dims(self):
        xlo = float(self.beg[5].split()[0])
        xhi = float(self.beg[5].split()[1])
        
        ylo = float(self.beg[6].split()[0])
        yhi = float(self.beg[6].split()[1])
        
        zlo = float(self.beg[7].split()[0])
        zhi = float(self.beg[7].split()[1])
        
        self.lx = xhi-xlo
        self.ly = yhi-ylo
        self.lz = zhi-zlo
        
    def get_periodicity(self):
        self.periodicity = [True,True,True]
        split = self.beg[4].split()
        for i in [6,7,8]:
            if split[i] != 'pp':
                self.periodicity[i-6] = False
    
    
    def get_atom_distance(self,a1,a2):
        rx = a1.x - a2.x
        ry = a1.y - a2.y
        rz = a1.z - a2.z
        lx = self.lx
        ly = self.ly
        lz = self.lz
        
        if self.periodicity[0]:
            if rx > lx/2:
                rx = rx-lx
            if rx <= -lx/2:
                rx = rx+lx
        
        if self.periodicity[1]:
            if ry > ly/2:
                ry = ry-ly
            if ry <= -ly/2:
                ry = ry+ly
                
        if self.periodicity[2]:
            if rz > lz/2:
                rz = rz-lz
            if rz <= -lz/2:
                rz = rz+lz
        
        distance = (rx**2+ry**2+rz**2)**0.5
        return distance
        
    
    def make_neighbor_dic(self,combos):
        values = list(itertools.repeat([], len(selt.atom_array)))
        self.neighbor_dic = dict(zip(self.atom_array,values))
        
        
        for c in combos:
            ai,aj = c
        
            d = self.get_atom_distance(ai,aj)
            if ai.atom_type == 1 and aj.atom_type == 1:
                if d < 4.07:
                    self.neighbor_dic[ai].append(aj)
                    self.neighbor_dic[aj].append(ai)
            if ai.atom_type == 1 and aj.atom_type == 2:
                if d < 3.52:
                    self.neighbor_dic[ai].append(aj)
                    self.neighbor_dic[aj].append(ai)
            if ai.atom_type == 2 and aj.atom_type == 3:
                if d < 2.86:
                    self.neighbor_dic[ai].append(aj)
                    self.neighbor_dic[aj].append(ai)
        return self.neighbor_dic
    
    # def get_dimensions(self):
    #     self.zs.sort()
    #     top = np.average(self.zs[-1*self.nAtoms//self.basis:])
    #     bottom = np.average(self.zs[0:self.nAtoms//self.basis])
    #     self.average_dz = top - bottom
        
    #     self.dy = max(self.ys) - min(self.ys)
        
    #     self.dx = max(self.xs)
        
    #     return [self.dx, self.dy, self.average_dz]
    
    # def get_area(self):
    #     return self.dx * self.dy
    
    
    # def get_simulation_box_volume(self):
    #     with open(self.filename) as f:
    #         lines = f.readlines()

    #     a1 = float(lines[6-1].split()[1])
    #     b1 = float(lines[6-1].split()[0])
    #     b2 = float(lines[7-1].split()[1])
    #     c3 = float(lines[8-1].split()[1])
    #     a = [a1, 0, 0]
    #     b = [b1, b2, 0]
    #     c = [0, 0, c3]

    #     self.simulation_box_volume = np.linalg.det(np.dstack([a,b,c]))[0]
    #     return self.simulation_box_volume
    
    # def get_approximate_volume(self):
    #     return self.dy*self.dx*self.average_dz

    # def in_box(self,xmin,xmax,ymin,ymax,zmin,zmax):
    #     atoms_in_box = []
    #     for atom in self.atom_array:
    #         x = atom.pos[0]
    #         y = atom.pos[1]
    #         z = atom.pos[2]
    #         if x < xmax and x > xmin:    
    #             if y < ymax and y > ymin:
    #                 if z < zmax and z > zmin:
    #                     atoms_in_box.append(atom)
    #     return atoms_in_box
    
    
    
        