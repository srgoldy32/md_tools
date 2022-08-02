import random
import numpy as np
from Atom import Atom

class Mesh:
    def __init__(self, ident, unit_cell, dupx, dupy, dupz):
        self.ident = ident
        self.unit_cell = unit_cell
        self.dupx = dupx
        self.dupy = dupy
        self.dupz = dupz
        self.create_mesh()
        self.nAtoms = len(self.atom_list)
        self.sort_atoms()
        
    def __str__(self):
        return "<Mesh: {}, nAtoms: {}>".format(self.ident,self.nAtoms)
    
    def generate_datafile(self,filename):
        self.nAtoms = len(self.atom_list)
        f = open(filename, "w")
        f.write("LAMMPS data file\n")
        f.write("\n")
        f.write("{} atoms\n".format(self.nAtoms))
        f.write("\n")
        f.write("{} atom types\n".format(len(self.atom_dict)))
        f.write("\n")
        f.write("0 {} xlo xhi\n".format(self.dupx*float(self.unit_cell.dx)))
        f.write("0 {} ylo yhi\n".format(self.dupy*float(self.unit_cell.dy)))
        f.write("0 {} zlo zhi\n".format(self.dupz*float(self.unit_cell.dz)))
        f.write("{} 0 0 xy xz yz\n".format(self.dupy*float(self.unit_cell.dxy))) # need to do a different calculation here - multiply by dupx (?)
        f.write("\n")
        f.write("Masses\n")
        f.write("\n")
        for mass in self.unit_cell.mass_dict:
            f.write("{} {}\n".format(mass,self.unit_cell.mass_dict[mass]))
        f.write("\n")
        f.write("Atoms\n")
        f.write("\n")
        for atom in self.atom_list:
            f.write(atom.line_string()+"\n")

        f.close()
    
    def sort_atoms(self):
        self.atom_dict = {}
        for atom in self.atom_list:
            if self.atom_dict.get(atom.atom_type) == None:
                self.atom_dict[atom.atom_type] = [atom]
            else:
                self.atom_dict[atom.atom_type].append(atom)
                
    def remove_by_ident(self,remove_ident):
        for atom in self.atom_list:
            if atom.ident == remove_ident:
                to_be = atom
                self.atom_list.remove(to_be)
                self.sort_atoms()
                return "Atom {} removed, nAtoms: {}".format(remove_ident,len(self.atom_list))
        return "Atom {} not found, nAtoms: {}".format(ident,len(self.atom_list))
    
    def remove_atom_type(self, remove_type):
        choice_list = []
        if self.atom_dict.get(remove_type) == None:
            return "No atoms of that type"
        for atom in self.atom_dict[remove_type]:
            choice_list.append(atom.ident)
        to_be_ident = random.choice(choice_list)
        return self.remove_by_ident(to_be_ident)

    def atom_distance(self,atom_a, atom_b):
        x = abs(atom_a.pos[0]-atom_b.pos[0])**2
        y = abs(atom_a.pos[1]-atom_b.pos[1])**2
        z = abs(atom_a.pos[2]-atom_b.pos[2])**2
        dist = (x+y+z)**.5
        return dist
        
    def check_overlap(self, tol):
        dists = []
        for a in self.atom_list:
            for b in self.atom_list:
                if a.ident != b.ident:
                    dist = self.atom_distance(a,b)
                    dists.append(dist)
                    if dist < tol:
                        print("{} between {} and {}".format(dist,a.ident,b.ident))
        return "Min distance: {}".format(min(dists))
    
    def create_mesh(self):
        unit_cell_atoms = self.unit_cell.atom_array
        dx = float(self.unit_cell.dx)
        dy = float(self.unit_cell.dy)
        dxy = float(self.unit_cell.dxy)
        count = 1
        self.atom_list = []

        for j in range(self.dupy):
            for i in range(self.dupx):
                for atom in self.unit_cell.atom_array:
                    new_atom = Atom(count, atom.atom_type, [atom.pos[0]+dx*i+j*dxy,atom.pos[1]+dy*j,atom.pos[2]])


                    self.atom_list.append(new_atom)
                    count+=1

        
                    
    def remove_percentage(self,remove_type, percent):
        if self.atom_dict.get(remove_type) == None:
            return "No atoms of that type"
        current = len(self.atom_dict.get(remove_type))
        cut = round((percent/100.0 * current))
        for i in range(cut):
            print(self.remove_atom_type(remove_type))
        
    def get_simulation_box_volume(self):
        a = [self.dupx*float(self.unit_cell.dx), 0, 0]
        b = [self.dupy*float(self.unit_cell.dxy), self.dupy*float(self.unit_cell.dy), 0]
        c = [0, 0, self.unit_cell.z_max]
        
        self.simulation_box_volume = np.linalg.det(np.dstack([a,b,c]))[0]
        return self.simulation_box_volume
    
    def get_dimensions(self):
        xs = []
        ys = []
        zs = []
        x_min = 0.0
        top_left_y = 0.0
        for atom in self.atom_list:
            x = atom.pos[0]
            y = atom.pos[1]
            z = atom.pos[2]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            if x < x_min:
                x_min = x
                top_left_y = y
            
        basis = 5 # for Ti2
        zs.sort()
        top = np.average(zs[-1*self.nAtoms//basis:])
        bottom = np.average(zs[0:self.nAtoms//basis])
        self.average_dz = top - bottom

        dy = max(ys) - min(ys)

        dx = max(xs) - min(xs)
        self.approximate_volume = dy*max(xs) * self.average_dz
        
        self.x_dim = max(xs)
        self.y_dim = (top_left_y**2 + x_min**2)**.5
        return [self.x_dim, self.y_dim]

        