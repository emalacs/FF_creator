from read_input_files import *
from functions import *

pep_atoms = read_pep_atoms()
dict_pep_atomtypes = make_atp_dict(pep_atoms)

pep_pairs = read_pep_pairs()
fib_pairs = read_fib_pairs()

fib_atoms = read_fib_atoms()
dict_fib_atomtypes = make_atp_dict(fib_atoms)

merge_pairs = ffnonbonded_merge_pairs(pep_pairs, fib_pairs, dict_pep_atomtypes, dict_fib_atomtypes)

print(merge_pairs)

exit()


pep_atoms = read_pep_atoms()
dict_pep_atomtypes = make_atp_dict(pep_atoms)

pep_merge_dihedrals = read_pep_dihedrals()

fib_atoms = read_fib_atoms()
dict_fib_atomtypes = make_atp_dict(fib_atoms)

fib_merge_dihedrals = read_fib_dihedrals()

ffbonded_merge_dihedrals(pep_merge_dihedrals, fib_merge_dihedrals, dict_pep_atomtypes, dict_fib_atomtypes)

#print(fib_dihedrals)