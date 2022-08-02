class Atom:
    def __init__(self, ident, atom_type, pos):
        self.ident = ident
        self.atom_type = atom_type
        self.pos = pos
    def __str__(self):
        return "<Atom: {}, type: {}>".format(self.ident,self.atom_type)
    def line_string(self):
        return "{} {} {} {} {}".format(self.ident, self.atom_type, self.pos[0], self.pos[1], self.pos[2])