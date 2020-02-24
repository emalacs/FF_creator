from read_input_files import *
from functions import *
from output import *
from atomtypes_aa_definitions import *

# This first part is derived from PSCRIPT_ATP_TOP.PY
# which creates the bonds, angles and dihedrals FROM THE PEPTIDE SMOG to paste into the new topology file
# Topology bonds, angles and dihedrals

print('Topology file to paste into .top')
print('Topology atoms, bonds, angles and dihedrals')

write_topology_atoms(read_pep_atoms())

# vimdiff ok

top_bonds = make_topology_bonds(read_pep_bonds())
write_topology_bonds(top_bonds)

# vimdiff ok

top_angles = make_topology_angles(read_pep_angles())
write_topology_angles(top_angles)

# vimdiff ok

top_dihedrals = make_topology_dihedrals(read_pep_dihedrals())
write_topology_dihedrals(top_dihedrals)

# vimdiff ok
# ok up to here after script reorganization

print('Topology files created')
print('')

# This script is derived from PSCRIPT_FF.PY
# which taking the bonds, angles and bonds from the PEPTIDE SMOG and prepare the FF

# atomtype definition

print('atomtype.atp creation')

pep_atoms = read_pep_atoms()
atp, atomtypes, dict_pep_atomtypes = make_atomtypes_and_dict(pep_atoms)
# DICTIONARY MASS
write_atomtypes_atp(atp)

# vimdiff ok, cambiate tutte le masse

print('atomtypes.atp written')
print('')

# ffbonded.itp Requires bonds, angles and dihedrals

print('Preparing the peptide ffbonded.itp')

# Making a dictionary out of it to change the atomnumber to the atomtype
# The dictionary is based on the peptide atoms

# Covalent preparation
pep_bonds = read_pep_bonds()
pep_ff_bonds = ffbonded_bonds(pep_bonds, dict_pep_atomtypes)

pep_angles = read_pep_angles()
pep_ff_angles = ffbonded_angles(pep_angles, dict_pep_atomtypes)

pep_dihedrals = read_pep_dihedrals()
pep_ff_dihedrals = ffbonded_dihedrals(pep_dihedrals, dict_pep_atomtypes)
write_pep_ffbonded(pep_ff_bonds, pep_ff_angles, pep_ff_dihedrals)

print('Peptide ffbonded.itp created')
print('')
# vimdiff ok

# Nonbonded preparation
print('Preparing the peptide ffnonbonded.itp')
# Addition of parameters into atomtypes
# DICTIONARY C12
# Preparing the pep_pairs
pep_pairs = read_pep_pairs()
pep_ff_pairs = ffnonbonded_pep_pairs(pep_pairs, dict_pep_atomtypes)
write_pep_ffnonbonded(atomtypes, pep_ff_pairs)

print('Peptide ffnonbonded.itp created')
print('')
print('Peptide topology to paste and FF ready!!')
print('')

# vimdiff ok, cambiati massa e c12

print('Preparing the fibril ffbonded.itp')

# Making a dictionary out of it to change the atomnumber to the atomtype
# The dictionary is based on the fibril atoms
fib_atoms = read_fib_atoms()
atp, atomtypes, dict_fib_atomtypes = make_atomtypes_and_dict(fib_atoms)

# vimdiff ok, sono lo stesso file che c'e nel peptide

# Covalent preparation
fib_bonds = read_fib_bonds()
fib_ff_bonds = ffbonded_bonds(fib_bonds, dict_fib_atomtypes)

fib_angles = read_fib_angles()
fib_ff_angles = ffbonded_angles(fib_angles, dict_fib_atomtypes)

fib_dihedrals = read_fib_dihedrals()
fib_ff_dihedrals = ffbonded_dihedrals(fib_dihedrals, dict_fib_atomtypes)

write_fib_ffbonded(fib_ff_bonds, fib_ff_angles, fib_ff_dihedrals)

print('Fibril ffbonded.itp created')
print('')
# vimdiff ok

print('Preparation of fibril ffnonbonded.itp')

# Preparing the fib_pairs
# OCCHIO A QUESTA PARTE PERCHE CI SONO DELLE DIFFERENZE NEL PAIRS
fib_pairs = read_fib_pairs()
fib_ff_pairs = ffnonbonded_fib_pairs(fib_pairs, dict_fib_atomtypes)
# DICTIONARY C12
write_fib_ffnonbonded(atomtypes, fib_ff_pairs)

# vimdiff ok

print('Fibril ffnonbonded.itp created')
print('')
print('Fibril FF ready!!')
print('')
print('Merge ffbonded.itp preparation')

merge_dihedrals = ffbonded_merge_dihedrals(pep_dihedrals, fib_dihedrals, dict_pep_atomtypes, dict_fib_atomtypes)
write_merge_ffbonded(pep_ff_bonds, pep_ff_angles, merge_dihedrals)

# vimdiff ok

print('Merge ffbonded.itp created')
print('')

print('Merge ffnonbonded.itp preparation')

merge_pairs = ffnonbonded_merge_pairs(pep_pairs, fib_pairs, dict_pep_atomtypes, dict_fib_atomtypes)
write_merge_ffnonbonded(atomtypes, merge_pairs)

# vimdiff ok

print('Merge ffnonbonded.itp created')
print('')
print('All the three Force Fields are ready. Charlie is super happy!!!')
