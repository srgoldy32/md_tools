class Atom:
    def __init__(self, ident, atom_type, pos, extra_dict):
        self.ident = ident
        self.atom_type = atom_type
        self.pos = pos
        self.extra_dict = extra_dict
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
    def __str__(self):
        return "<Atom: {}, type: {}>".format(self.ident,self.atom_type)
    def line_string(self):
        return "{} {} {} {} {}".format(self.ident, self.atom_type, self.x, self.y, self.z)