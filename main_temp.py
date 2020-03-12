from read_input_files import *
from functions import *
from output import *
from atomtypes_aa_definitions import *

# SAREBBE DA RIFARE LA GESTIONE DELL'INPUT. IN QUESTO MODO VIENE LETTO TUTTO UNA VOLTA, MODIFICATO UNA VOLTA E POI
# LE DIFFERENZE STANNO NELL'OUTPUT

########################################################################################################################
########################################################################################################################
# TEST GROMOS

#print(read_gro_atoms().to_string())

# Dictionaries from gromos topology
#gromos_mass_dict, gromos_res_atom_dict, dict_gro_atomtypes = make_gromos_topology_dictionaries(read_gro_atoms(), read_gro_impropers())

########################################################################################################################
########################################################################################################################


# This first part is derived from PSCRIPT_ATP_TOP.PY
# which creates the bonds, angles and dihedrals FROM THE PEPTIDE SMOG to paste into the new topology file
# Topology bonds, angles and dihedrals

print('Topology file to paste into .top')
print('Topology atoms, bonds, angles and dihedrals')

write_topology_atoms(read_pep_atoms())

top_bonds = make_topology_bonds(read_pep_bonds())
write_topology_bonds(top_bonds)

top_angles = make_topology_angles(read_pep_angles())
write_topology_angles(top_angles)

top_dihedrals = make_topology_dihedrals(read_pep_dihedrals())
write_topology_dihedrals(top_dihedrals)

print('Topology files created')
print('')

# atomtype definition

print('atomtype.atp creation')

#pep_atoms = read_pep_atoms()
atp, atomtypes, dict_pep_atomtypes, dict_pep_aminores = make_atomtypes_and_dict(read_pep_atoms())

write_atomtypes_atp(atp)

print('atomtypes.atp written')
print('')

# ffbonded.itp Requires bonds, angles and dihedrals

print('Preparing the peptide ffbonded.itp')

# Making a dictionary out of it to change the atomnumber to the atomtype
# The dictionary is based on the peptide atoms

# Covalent preparation
#pep_bonds = read_pep_bonds()
pep_ff_bonds = ffbonded_bonds(read_pep_bonds(), dict_pep_atomtypes, dict_pep_aminores)

#pep_angles = read_pep_angles()
pep_ff_angles = ffbonded_angles(read_pep_angles(), dict_pep_atomtypes, dict_pep_aminores)

#pep_dihedrals = read_pep_dihedrals()
pep_ff_dihedrals = ffbonded_dihedrals(read_pep_dihedrals(), dict_pep_atomtypes, dict_pep_aminores)

write_pep_ffbonded(pep_ff_bonds, pep_ff_angles, pep_ff_dihedrals)

print('Peptide ffbonded.itp created')
print('')

# Nonbonded preparation
print('Preparing the peptide ffnonbonded.itp')
# Addition of parameters into atomtypes
# DICTIONARY C12
# Preparing the pep_pairs
#pep_pairs = read_pep_pairs()
pep_ff_pairs = ffnonbonded_pep_pairs(read_pep_pairs(), dict_pep_atomtypes)

write_pep_ffnonbonded(atomtypes, pep_ff_pairs)

print('Peptide ffnonbonded.itp created')
print('')
print('Peptide topology to paste and FF ready!!')
print('')

print('Preparing the fibril ffbonded.itp')

# Making a dictionary out of it to change the atomnumber to the atomtype
# The dictionary is based on the fibril atoms
#fib_atoms = read_fib_atoms()
atp, atomtypes, dict_fib_atomtypes, dict_fib_aminores = make_atomtypes_and_dict(read_fib_atoms())

# vimdiff ok, sono lo stesso file che c'e nel peptide

# Covalent preparation
#fib_bonds = read_fib_bonds()
fib_ff_bonds = ffbonded_bonds(read_fib_bonds(), dict_fib_atomtypes, dict_fib_aminores)

#fib_angles = read_fib_angles()
fib_ff_angles = ffbonded_angles(read_fib_angles(), dict_fib_atomtypes, dict_fib_aminores)

#fib_dihedrals = read_fib_dihedrals()
fib_ff_dihedrals = ffbonded_dihedrals(read_fib_dihedrals(), dict_fib_atomtypes, dict_fib_aminores)

write_fib_ffbonded(fib_ff_bonds, fib_ff_angles, fib_ff_dihedrals)

print('Fibril ffbonded.itp created')
print('')

print('Preparation of fibril ffnonbonded.itp')

# Preparing the fib_pairs
# OCCHIO A QUESTA PARTE PERCHE CI SONO DELLE DIFFERENZE NEL PAIRS
#fib_pairs = read_fib_pairs()
fib_ff_pairs = ffnonbonded_fib_pairs(read_fib_pairs(), dict_fib_atomtypes)

write_fib_ffnonbonded(atomtypes, fib_ff_pairs)

print('Fibril ffnonbonded.itp created')
print('')
print('Fibril FF ready!!')
print('')
print('Merge ffbonded.itp preparation')

merge_dihedrals = ffbonded_merge_dihedrals(read_pep_dihedrals(), read_fib_dihedrals(), dict_pep_atomtypes,
                                           dict_fib_atomtypes)
write_merge_ffbonded(pep_ff_bonds, pep_ff_angles, merge_dihedrals)

print('Merge ffbonded.itp created')
print('')

print('Merge ffnonbonded.itp preparation')

merge_pairs = ffnonbonded_merge_pairs(read_pep_pairs(), read_fib_pairs(), dict_pep_atomtypes, dict_fib_atomtypes)
write_merge_ffnonbonded(atomtypes, merge_pairs)

print('Merge ffnonbonded.itp created')
print('')
print('All the three Force Fields are ready. Charlie is super happy!!!')
