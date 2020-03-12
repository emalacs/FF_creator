import numpy as np
import pandas as pd

# INPUTS FROM SMOG PEPTIDE

def read_pep_atoms():
    pep_atoms = pd.read_csv('input/pep_atoms', sep = '\s+', header = None)
    # header=0 because it counts ; as a column. Therefore, I don't care about the header and I'll recreate one from
    # scratch.
    pep_atoms["charge"] = ""
    # addition of a charge empty column
    pep_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]
    # I've created the header with the correct column names
    # The following function is necessary to associate the resID with the residue number.
    # Changing the atom type column -> patoms column 5 + 3
    pep_atoms["type"] = pep_atoms["atom"].apply(str) + '_' + pep_atoms["resnr"].apply(str)
    # Afterwards the function make_pep_atp_dict will create a dictionary based on this column
    return pep_atoms


def read_pep_bonds():
    # Reading the peptide bonds
    pep_bonds = pd.read_csv('input/pep_bonds', sep = "\s+", header = None)
    pep_bonds.columns = [";ai", "aj", "func", "r0", "kb"]
    return pep_bonds


def read_pep_angles():
    # Reading the peptide angles
    pep_angles = pd.read_csv('input/pep_angles', sep = "\s+", header = None)
    pep_angles.columns = [";ai", "aj", "ak", "func", "th0", "Ka"]
    return pep_angles


def read_pep_dihedrals():
    # Reading the peptide dihedrals
    pep_dihedrals = pd.read_csv('input/pep_dihedrals', sep = "\s+", header = None)
    pep_dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
    pep_dihedrals['mult'] = pep_dihedrals['mult'].fillna(value = '')
    # dihedrals = dihedrals.replace(np.nan, '', regex=True)
    # NaN are replaced with and empty line.
    return pep_dihedrals


def read_pep_pairs():
    # Reading the peptide pairs
    pep_pairs = pd.read_csv('input/pep_pairs', sep = "\s+", header = None)
    pep_pairs.columns = [";ai", "aj", "type", "A", "B"]
    return pep_pairs


# INPUTS FROM SMOG FIBRIL

def read_fib_atoms():
    fib_atoms = pd.read_csv('input/fib_atoms', sep = '\s+', header = None)
    # header=0 because it counts ; as a column.
    # Therefore, I don't care about the header and I'll recreate one from scratch.
    fib_atoms["charge"] = ""
    # addition of a charge empty column
    fib_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]
    fib_atoms['resnr'] = fib_atoms['resnr'] % 11
    fib_atoms['resnr'] = fib_atoms['resnr'].replace(0, 11)
    # I've created the header with the correct column names
    fib_atoms["type"] = fib_atoms["atom"].apply(str) + '_' + fib_atoms["resnr"].apply(str)
    return fib_atoms


def read_fib_bonds():
    # Reading the fib_atomstide bonds
    fib_bonds = pd.read_csv('input/fib_bonds', sep = "\s+", header = None)
    fib_bonds.columns = [";ai", "aj", "func", "r0", "kb"]
    return fib_bonds


def read_fib_angles():
    # Reading the fib_atomstide angles
    fib_angles = pd.read_csv('input/fib_angles', sep = "\s+", header = None)
    fib_angles.columns = [";ai", "aj", "ak", "func", "th0", "Ka"]
    return fib_angles


def read_fib_dihedrals():
    # Reading the fib_atomstide dihedrals
    fib_dihedrals = pd.read_csv('input/fib_dihedrals', sep = "\s+", header = None)
    fib_dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
    fib_dihedrals['mult'] = fib_dihedrals['mult'].fillna(value = '')
    # NaN are replaced with and empty line.
    return fib_dihedrals


def read_fib_pairs():
    # Reading the fib_atomstide pairs
    fib_pairs = pd.read_csv('input/fib_pairs', sep = "\s+", header = None)
    fib_pairs.columns = [";ai", "aj", "type", "A", "B"]
    return fib_pairs

########################################################################################################################
########################################################################################################################

# GROMOS TEST


def read_gro_atoms():
    # Reading the atoms section from gromos topology
    gro_atoms = pd.read_csv('input/pep_gro_atoms', sep = "\s+", header = None)
    gro_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", 'charge', 'mass']
    gro_atoms = gro_atoms.drop(['charge', 'cgnr'], axis = 1)
    gro_atoms["atom_nmr"] = gro_atoms["atom"].apply(str) + '_' + gro_atoms["resnr"].apply(str)
    gro_atoms['res_atom'] = gro_atoms['residue'] + '_' + gro_atoms['atom']
    return gro_atoms


def read_gro_impropers():
    # Reading the impropers obtained from gromos topology
    gro_impropers = pd.read_csv('input/pep_gro_impropers', sep="\s+", header=None)
    gro_impropers.columns = ["; ai", "aj", "ak", "al", "funct", "define"]
    return gro_impropers