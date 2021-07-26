#!/usr/bin/env python3
import json
import sys
from pw_parse import Atom, Cell, AXSF

# map.json:
#   must be surjective
#   must be 1-indexed
#   maps input atoms to output atoms

def main(args):
    #Fill in these fields
    #Data folder must end with slash. This is where axsf files are read from
    # and written to.
    data_folder = "/home/phonons/results/"
    cell_type = "simple" #Output cell style. Options are 
                            #'full', 'simple', and (for x,m) 'super'
    compound = "bisco"
    prefix = "g" #Name of Brillouin zone point (only 'x', 'm' influence this code)
    extra = "_15k" #Used when running multiple calculations for same compound/q-point
    dynmat = True #Specifies whether to reset forces by the dynmat.out file
                    #This is recommended, as the axsf file may not properly
                    # distinguish _u and _g modes.
    #make sure file names match these defaults
    axsf_file = prefix + "_" + compound + extra + ".axsf"
    dynmat_file = prefix + "_" + compound + extra + ".dynmat.out"
    pw_in = cell_type + "-" + compound +  ".in"
    map_file = "map_sg-" + cell_type + "-" + compound + ".json"
    
    #Begin conversion -- this code shouldn't need alteration
    sg_axsf = AXSF.parse_axsf(data_folder + axsf_file)
    if dynmat:
        sg_axsf.parse_dynmat(data_folder + dynmat_file)
    sg_axsf.orthogonalize_forces()
    for cell in sg_axsf.cells:
        cell.reduce_cell()

    try:
        f = open(map_file, "r")
        map_str = f.read()
        f.close()
    except Exception as e:
        print("error: map file '%s' not found"%map_file, file=sys.stderr)

    try:
        atom_map = json.loads(map_str)
    except Exception as e:
        print("error: couldn't parse map file '%s'"%map_file, file=sys.stderr)

    new_cell = Cell.parse_pw_in(pw_in)
    sg_axsf.map_unit_cell(new_cell, atom_map)

    if cell_type == "super":
        cells = sg_axsf.cells
        for j in range(len(cells)):
            cell = cells[j]
            base_length = len(cell.atoms)//4
            if prefix == "x":
                s00,s10,s01,s11 = False,True,False,True
            elif prefix == "m":
                s00,s10,s01,s11 = False,True,True,False
            for i in range(len(cell.atoms)):
                atom = cell.atoms[i]
                num = i // base_length
                if s00 and num == 0:
                    atom.Fx, atom.Fy, atom.Fz = -atom.Fx, -atom.Fy, -atom.Fz
                elif s10 and num == 1:
                    atom.Fx, atom.Fy, atom.Fz = -atom.Fx, -atom.Fy, -atom.Fz
                elif s01 and num == 2:
                    atom.Fx, atom.Fy, atom.Fz = -atom.Fx, -atom.Fy, -atom.Fz
                elif s11 and num == 3:
                    atom.Fx, atom.Fy, atom.Fz = -atom.Fx, -atom.Fy, -atom.Fz
    outfile = open(data_folder+prefix+'_'+cell_type+'_'+compound+extra+'.axsf','w')
    outfile.write(repr(sg_axsf))
    outfile.close()

if __name__ == "__main__":
    main(sys.argv[1:])
