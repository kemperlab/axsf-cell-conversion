#!/usr/bin/env python3
import re
import numpy as np
import os.path
import sys

def error(s):
    print(s, file=sys.stderr)


def warn(s):
    print(s, file=sys.stderr)
################################################################################
class Atom:
    def __init__(self, name="Cu", x=0.0, y=0.0, z=0.0, Fx=0.0, Fy=0.0, Fz=0.0):
        "Atoms are defined to have a name, position vector, and force vector"
        self.name = name
        self.x,  self.y,  self.z  = x,  y,  z
        self.Fx, self.Fy, self.Fz = Fx, Fy, Fz
    def __repr__(self):
        return "%s: x=(%.6f,%.6f,%.6f) F=(%.6f,%.6f,%.6f)" % \
            (self.name, self.x, self.y, self.z, self.Fx, self.Fy, self.Fz)
    @classmethod
    def parse_pw_in(cls, line):
        "A method called when parsing pw.x input files, creates an atom \
            with 0 forces"
        vec = line.strip().split()
        if len(vec) != 4:
            error("error: couldn't parse atom position: '%s'"%line)
        try:
            name, x, y, z = vec[0], *[float(i) for i in vec[1:]]
        except Exception as e:
            error("error: couldn't parse float in atom position: '%s'"%line)
        return cls(name, x, y, z)
    @classmethod
    def copy(cls, atom):
        "Duplicates given atom"
        return cls(atom.name, atom.x, atom.y, atom.z, atom.Fx, atom.Fy, atom.Fz)
################################################################################
class Cell:
    def __init__(self, par="angstrom", pos="crystal", ux=1.0, uy=0.0, uz=0.0, \
            vx=0.0, vy=1.0, vz=0.0, wx=0.0, wy=0.0, wz=0.0, celldm=0.0):
        "Cells have basis vectors and contain some number of atoms"
        self.par = par #Sets the units of cell parameters (unit cell vectors)
        self.pos = pos #Sets the units of the atom positions
        #Sets values of unit cell vectors
        self.ux, self.uy, self.uz = ux, uy, uz
        self.vx, self.vy, self.vz = vx, vy, vz
        self.wx, self.wy, self.wz = wx, wy, wz
        #Sets alat parameter
        self.celldm = celldm
        if celldm == 0.0:
            #By default alat is set to the x magnitude of the first vector
            self.celldm = ux*self.get_par_convert_factor(par, "bohr")
        self.atoms = list()
    def __repr__(self):
        par = self.par
        par = "alat celldm(1) = %.6f"%self.celldm if par == "alat" else par
        #Sets how to print the cell in pw.x input format
        pars = """\
CELL_PARAMETERS %s
    %.6f  %.6f  %.6f
    %.6f  %.6f  %.6f
    %.6f  %.6f  %.6f

ATOMIC_POSITIONS %s
""" % (par, self.ux, self.uy, self.uz, self.vx, self.vy, self.vz, \
            self.wx, self.wy, self.wz, self.pos)
        atom_strs = list()
        #For nice formatting, adds a space if the number is not negative
        sp = lambda x: "" if x < 0 else " "
        for atom in self.atoms:
            atom_strs.append("    %-2s %s%.6f %s%.6f %s%.6f" % (atom.name, \
                sp(atom.x), atom.x, sp(atom.y), atom.y, sp(atom.z), atom.z))
        return pars + "\n".join(atom_strs)
    @classmethod
    def parse_pw_in(cls, txt):
        "A method which parses pw.x input files into a cell object"
        if os.path.isfile(txt):
            f = open(txt, "r")
            txt = f.read()
            f.close()
        txt = txt.split("\n\n")
        atom_strs = list()
        cell_pos = list()
        celldm = 0.0
        found_par = False
        found_pos = False
        #Parses text for celldm value and CELL_PAR, ATOMIC_POS sections
        for section in txt:
            #Search for celldm value
            r = re.search(r"celldm\(1\)\s*?\=\s*?([0-9\.]+)", section)
            if r:
                celldm = r.group(1)
                try:
                    celldm = float(celldm)
                except Exception as e:
                    error("error: couldn't parse celldm: '%s'"%celldm)
            #Search for CELL_PARAMETERS section and extract cell vectors
            r = re.search(r"^CELL_PARAMETERS\s+(\w+)\s*", section)
            if r:
                par = r.group(1)
                for line in section.split("\n")[1:]:
                    vec = line.split()
                    if len(vec) != 3:
                        error("error: couldn't parse cell parameter: '%s'"%line)
                    try:
                        vec = [float(i) for i in vec]
                    except Exception as e:
                        error("error: couldn't parse cell parameter: '%s'"%line)
                    cell_pos.append(vec)
                cell_pos_len = len(cell_pos)
                if cell_pos_len != 3:
                    error("error: found %d cell parameter vectors"%cell_pos_len)
                found_par = True
                #Search for ATOMIC_POSITIONS section and extract atom locations
            r = re.search(r"^ATOMIC_POSITIONS\s+(.*)", section)
            if r:
                pos = r.group(1)
                for line in section.split("\n")[1:]:
                    atom_strs.append(line)
                found_pos = True
        if found_par and found_pos:
            ux, uy, uz = cell_pos[0]
            vx, vy, vz = cell_pos[1]
            wx, wy, wz = cell_pos[2]
            if par == "alat" and celldm == 0.0:
                error("error: alat specified and celldm not found")
            cell = cls(par, pos, ux, uy, uz, vx, vy, vz, wx, wy, wz, celldm)
            for atom_str in atom_strs:
                cell.add_atom(Atom.parse_pw_in(atom_str))
            return cell
        else:
            error("error: couldn't parse text")
    @classmethod
    def copy(cls, cell):
        "Creates new, identical cell object"
        par, pos = cell.par, cell.pos
        ux, uy, uz = cell.ux, cell.uy, cell.uz
        vx, vy, vz = cell.vx, cell.vy, cell.vz
        wx, wy, wz = cell.wx, cell.wy, cell.wz
        celldm = cell.celldm
        copy = cls(par, pos, ux, uy, uz, vx, vy, vz, wx, wy, wz, celldm)
        for atom in cell.atoms:
            copy.add_atom(Atom.copy(atom))
        return copy
    @classmethod
    def equivalent(cls, cell1, cell2):
        "[DEPRECATED] For use with AXSF conversion, tests if two cells' atoms are in \
        symmetrically equivalent positions."
        if len(cell1.atoms) <= len(cell2.atoms):
            small_cell = Cell.copy(cell1)
            large_cell = Cell.copy(cell2)
        else:
            small_cell = Cell.copy(cell2)
            large_cell = Cell.copy(cell1)
        multiple = False
        for i in range(1,10):
            if not i != len(large_cell.atoms) / len(small_cell.atoms):
                multiple = True
                factor = i
        if not multiple:
            return False, -1
        small_cell.set_par('angstrom')
        large_cell.set_par('angstrom')
        small_cell.set_pos('angstrom')
        large_cell.set_pos('angstrom')
        for i in range(0,len(small_cell.atoms)):
            s_atom = small_cell.atoms[i]
            l_atom = large_cell.atoms[i]
            cell = small_cell
            if s_atom.name != l_atom.name:
                return False, i
            elif s_atom.z - cell.wz < -.001 and l_atom.z - cell.wz < -.001:
                return abs(s_atom.z - l_atom.z) < .001, i
            elif s_atom.z - cell.wz < -.001:
                return abs(s_atom.z - (l_atom.z / cell.wz - (factor-1))) < .001, i
            elif l_atom.z - cell.wz < -.001:
                return abs(l_atom.z - (s_atom.z / cell.wz - 1)) < .001, i
            else:
                return abs((l_atom.z / cell.wz - (factor-1)) \
                           - (s_atom.z / cell.wz - 1)) < .001, i
    def add_atom(self, atom):
        self.atoms.append(atom)
        return atom
    def get_atom(self, id=0):
        return self.atoms[id]
    def get_par_convert_factor(self, p1, p2):
        "Allows conversion between different specifications of 'par'"
        a0 = 0.52918
        r = [1.0, a0, a0*self.celldm]
        s = ["angstrom", "bohr", "alat"]
        try:
            return r[s.index(p1)]/r[s.index(p2)]
        except Exception as e:
            error("error: couldn't convert cell par '%s' to '%s'" % (p1, p2))
    def get_pos_convert_matrix(self, p1, p2):
        "Allows conversion between different specifications of 'pos'"
        a0 = 0.52918
        s = ["angstrom", "bohr", "alat", "crystal"]
        I = np.matrix(np.eye(3))
        A = self.get_par_convert_factor(self.par, "angstrom") * np.matrix( \
            [[self.ux, self.vx, self.wx], \
             [self.uy, self.vy, self.wy], \
             [self.uz, self.vz, self.wz]])
        R = [I, a0*I, a0*self.celldm*I, A]
        return R[s.index(p2)].I*R[s.index(p1)]
    def set_par(self, par):
        "Converts between different 'par'"
        r = self.get_par_convert_factor(self.par, par)
        self.ux, self.uy, self.uz = r*self.ux, r*self.uy, r*self.uz
        self.vx, self.vy, self.vz = r*self.vx, r*self.vy, r*self.vz
        self.wx, self.wy, self.wz = r*self.wx, r*self.wy, r*self.wz
        self.par = par
    def set_pos(self, pos):
        "Converts between different 'pos'"
        R = self.get_pos_convert_matrix(self.pos, pos)
        for atom in self.atoms:
            atom_xyz = (R*np.matrix([atom.x, atom.y, atom.z]).T).A1
            atom.x, atom.y, atom.z = atom_xyz
        self.pos = pos
    def reduce_cell(self):
        "Reduces an ibrav=7 unit cell to simple vectors"
        self.set_par("angstrom")
        self.set_pos("angstrom")
        #Change unit vectors
        self.ux, self.uy, self.uz = (2*abs(self.ux), 0, 0)
        self.vx, self.vy, self.vz = (0,2*abs(self.vy), 0)
        self.wx, self.wy, self.wz = (abs(self.wx), abs(self.wy), abs(self.wz))
        for atom in self.atoms:
            if abs(atom.z-self.wz) < 0.001:
                atom.z = 0.0
            atom.x, atom.y, atom.z = abs(atom.x), abs(atom.y), abs(atom.z)
            if atom.z/self.wz > 1:
                atom.z = self.wz*((atom.z/self.wz)-1)
                atom.x = .5*self.ux
                atom.y = .5*self.vy
        self.set_pos("crystal")
    def align_forces(self):
        "Used in axsf conversion, rotates horizontal forces to a/b axes"
        #Finds first basis vector
        base = np.matrix([[1],[0]])
        for atom in self.atoms:
            if abs(atom.Fx) > 0.00001 or abs(atom.Fy) > 0.00001:
                base = np.matrix([[atom.Fx],[atom.Fy]])
                break
        base = base/np.linalg.norm(base)
        #Creates second basis vector
        base2 = np.matrix([[0,-1],[1,0]])*base
        #Transformation matrices between a/b basis and force basis
        T_b2e = np.matrix([[base[0,0],base2[0,0]],[base[1,0],base2[1,0]]])
        T_e2b = np.linalg.inv(T_b2e)
        for atom in self.atoms:
            neg = 1
            force = np.matrix([[atom.Fx],[atom.Fy]])
            force = T_e2b*force
            #Checks predominant force direction and assigns force to
            #solely that direction
            if abs(force[0,0]) >= abs(force[1,0]):
                if force[0,0] < 0:
                    neg = -1
                F = neg*np.sqrt(force[0,0]**2 + force[1,0]**2)
                force[0,0],force[1,0] = F,0.0
            else:
                if force[1,0] < 0:
                    neg = -1
                F = neg*np.sqrt(force[0,0]**2 + force[1,0]**2)
                force[0,0],force[1,0] = 0.0,F
            force = T_b2e*force
            atom.Fx,atom.Fy = force[0,0],force[1,0]

################################################################################
class AXSF():
    def __init__(self):
        self.cells = list()
    def __repr__(self):
        n = len(self.cells)
        if n < 1:
            return "ANIMSTEPS  0"
        c = self.cells[0]
        #Specifies .axsf formatting
        sp = lambda x: "  " if np.copysign(1.0,x) < 0 or abs(x) > 10 else "   "
        frames = """\
ANIMSTEPS  %d
CRYSTAL
PRIMVEC
    %6f    %6f    %6f
    %6f    %6f    %6f
    %6f    %6f %s%6f\
""" % (n, c.ux, c.uy, c.uz, c.vx, c.vy, c.vz, c.wx, c.wy, sp(c.wz), c.wz)
        for i in range(n):
            c = self.cells[i]
            frames += "\nPRIMCOORD  %d\n      %d   1" % (i+1, len(c.atoms))
            for atom in c.atoms:
                frames += "\n   %-2s  %s%.6f%s%.6f%s%.6f%s%.6f%s%.6f%s%.6f"\
                    % (atom.name, sp(atom.x), atom.x, sp(atom.y), atom.y, \
                    sp(atom.z), atom.z, sp(atom.Fx), atom.Fx, \
                    sp(atom.Fy), atom.Fy, sp(atom.Fz), atom.Fz)
        return frames
    @classmethod
    def parse_axsf(cls, txt):
        "Reads an axsf file and returns a new AXSF object"
        if os.path.isfile(txt):
            f = open(txt, "r")
            txt = f.read()
            f.close()
        else:
            error("error: file %s not found"%txt)
        ax = cls()
        txt = txt.split("PRIMCOORD")
        uvw = txt[0]
        uvw = uvw.split("PRIMVEC")[1]
        uvw = uvw.strip().split("\n")
        try:
            ux, uy, uz = [float(i) for i in uvw[0].split()]
            vx, vy, vz = [float(i) for i in uvw[1].split()]
            wx, wy, wz = [float(i) for i in uvw[2].split()]
        except Exception as e:
            error("error: couldn't parse cell parameters: '%s'"%txt)
        for block in txt[1:]:
            c = Cell("angstrom", "angstrom", ux, uy, uz, vx, vy, vz, wx, wy, wz)
            block = block.strip().split("\n")[2:]
            for atom_str in block:
                atom = [i.strip() for i in atom_str.strip().split()]
                if len(atom) != 7:
                    error("error: couldn't parse atom: '%s'"%atom_str)
                name = atom[0]
                try:
                   x, y, z, Fx, Fy, Fz = [float(i) for i in atom[1:]]
                except Exception as e:
                    error("error: couldn't convert to floats: '%s'"%atom_str)
                c.add_atom(Atom(name, x, y, z, Fx, Fy, Fz))
            ax.add_cell(c)
        return ax
    def parse_dynmat(self, txt):
        "Reads a .dynmat.out file and resets forces on existing axsf file"
        if os.path.isfile(txt):
            f = open(txt, "r")
            txt = f.read()
            f.close()
        else:
            error("error: file %s not found"%txt)
        txt = txt.split("freq")
        if len(txt)-len(self.cells) != 1:
            error("error: axsf and dynmat.out have different number of cells")
        for i in range(1,len(txt)):
            block = txt[i]
            c = self.cells[i-1]
            block = block.strip().split("\n")[1:]
            if i == len(txt)-1:
                block = block[:-1]
            for j in range(0,len(block)):
                force_str = block[j]
                force = [i.strip() for i in force_str.strip().split()]
                if len(force) != 8:
                    error("error: couldn't parse forces: '%s'"%force_str)
                try:
                   Fx, xx, Fy, yy, Fz, zz = [float(i) for i in force[1:-1]]
                except Exception as e:
                    error("error: couldn't convert to floats: '%s'"%force_str)
                #Corrects issue near x and m pts where forces are in cols 2,4,6
                if abs(xx) > abs(Fx):
                    Fx = xx
                if abs(yy) > abs(Fy):
                    Fy = yy
                if abs(zz) > abs(Fz):
                    Fz = zz
                Fx, Fy, Fz = round(0.1*Fx, 6), round(0.1*Fy, 6), round(0.1*Fz, 6)
                c.atoms[j].Fx, c.atoms[j].Fy, c.atoms[j].Fz = Fx, Fy, Fz
    def add_cell(self, cell):
        if cell.par != "angstrom":
            warn("warning: converting the cell parameters to angstroms")
            cell.set_par("angstrom")
        if cell.pos != "angstrom":
            warn("warning: converting the atom positions to angstroms")
            cell.set_pos("angstrom")
        self.cells.append(cell)
    def map_unit_cell(self, new_cell, atom_map):
        "Copies forces from current atoms and maps to atoms in 'new_cell.'"
        new_cell.set_par("angstrom")
        new_cell.set_pos("angstrom")
        for i in range(len(self.cells)):
            old_cell = self.cells[i]
            cell = Cell.copy(new_cell)
            # copy forces from new_cell to each output in map
            for j in range(len(atom_map)):
                old_atom = old_cell.atoms[j]
                if old_atom.name != atom_map[j]["name"]:
                    error("error: bad map line %d '%s' and '%s' don't match" \
                        % (j+1, old_atom.name, atom_map[j]["name"]))
                for k in atom_map[j]["outputs"]:
                    a = cell.atoms[k-1]
                    if old_atom.name != a.name:
                        error("error: bad map line %d output atom %d not '%s'" \
                            % (j+1, k, old_atom.name))
#                    equiv, num = Cell.equivalent(old_cell, cell)
#                    if not equiv:
#                        error("error: cells not equivalent. Hi. Check atom %d" % (num))
                    a.Fx, a.Fy, a.Fz = old_atom.Fx, old_atom.Fy, old_atom.Fz
            # free(self.cells[i])
            self.cells[i] = cell
    def orthogonalize_forces(self):
        "Determines whether the mode is vertical or horizontal, \
            then assigns each atom's forces to one axis."
        for cell in self.cells:
            #Determine if the mode's largest force is horizontal or vertical
            max_force = 0
            horizontal = True
            for atom in cell.atoms:
                Fh = np.sqrt(atom.Fx**2 + atom.Fy**2)
                if abs(atom.Fz) > Fh:
                    F = atom.Fz
                    atom.Fx, atom.Fy, atom.Fz = (0.0,0.0,F)
                    if abs(atom.Fz) > max_force:
                        max_force = abs(atom.Fz)
                        horizontal = False
                elif Fh > max_force:
                    max_force = Fh
                    horizontal = True
            #If horizontal, ensures forces are aligned with cell vectors
            if horizontal:
                cell.align_forces()
                for atom in cell.atoms:
                    neg = 1
                    if abs(atom.Fx) >= abs(atom.Fy):
                        if atom.Fx < 0:
                            neg = -1
                        F = neg*np.sqrt(atom.Fx**2 + atom.Fy**2)
                        atom.Fx, atom.Fy, atom.Fz = (F,0.0,0.0)
                    elif abs(atom.Fy) > abs(atom.Fx):
                        if atom.Fy < 0:
                            neg = -1
                        F = neg*np.sqrt(atom.Fx**2 + atom.Fy**2)
                        atom.Fx, atom.Fy, atom.Fz = (0.0,F,0.0)

################################################################################
