############################################################################
#######                  AXSF CELL CONVERSION CODE                   #######
############################################################################

This code transfers force vectors (from Quantum ESPRESSO's AXSF output)
from a Wigner-Seitz cell to a corresponding rectangular cell. Currently 
implemented for La2CuO4, Bi2Sr2CuO6, and Bi2Sr2CaCu2O8, but can be expanded
for use with other cells.


### RUNNING FOR EXISTING COMPOUNDS ###

1. Ensure AXSF and dynmat.out files match the structure used in code (or
modify code to match filenames) sg_conversion.py
2. Edit designated fields in sg_conversion.py to the desired cell type and
compound. 
    KEY: 
        'full'   - Rectangular unit cell containing all atoms
        'simple' - Smaller rectangular cell which eliminates doubled layers
                   found in studied compounds. Useful for reducing the 
                   height of unit cells for representation purposes.
        'super'  - 'simple' unit cell expanded into a 2 x 2 x 1 supercell.
                   For designated 'x' and 'm' q-points, code will alternate
                   direction of forces to represent short wavelength.
                   NOTE: Not implemented for Bi2212, as we did not have the
                   resources to do multiple calculations for this compound.
3. Run sg_conversion.py. File is created in same directory as axsf and 
dynmat.out files.


### USING NEW COMPOUNDS ###

1. Create a pw.x input file. Header not necessary but CELL_PARAMETERS and 
ATOMIC_POSITIONS are. This is what you want the final cell to look like.
2. Create a .json map designating which atoms in the final cell correspond
to which atoms in the Wigner-Seitz Cell. Use existing files as an example.


### FURTHER USES ###
Our pw_parse.py code is capable of parsing, modifying, and outputting pw.x
input files, including converting atomic positions from crystal units to 
angstrom or bohr, and vice versa. This can be useful when creating pw.x
input files when atomic positions are known in only one of these forms.
