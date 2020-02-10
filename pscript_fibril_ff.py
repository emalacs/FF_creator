#! /usr/bin/python

import numpy as np
import pandas as pd
import prettytable as pt
import copy
import re

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#inputs:
#atoms
#bonds
#angles
#dihedrals

## atoms = pd.read_csv('atoms', sep="\s+", header=None)   # header=0 ci piace ma conta il ; come colonna.
                                                        # provare a toglierlo durante l'importazione dei dati e aggiungerlo poi
                                                        # al momento si toglie tutto e vengono definite le colonne da me
## atoms["charge"] = ""                                    # addition of a charge empty column
## atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]

# changing the informations about the residue number. All chains must be the same.

## atoms["resnr"] = atoms["resnr"] % 11
## atoms["resnr"] = atoms["resnr"].replace(0, 11)

#changing the atom type column -> patoms column 5 + 3

## atoms["type"] = atoms["atom"].apply(str) + '_' + atoms["resnr"].apply(str)


## bonds = pd.read_csv('bonds', sep="\s+", header=None)    # stessa cosa di prima, bisogna togliere il nome delle colonne e reinserirle
## bonds.columns = [";ai", "aj", "func", "r0", "kb"]       # in attesa di un metodo piu intelligente

## angles = pd.read_csv('angles', sep="\s+", header=None)
## angles.columns = [";ai", "aj", "ak", "func", "th0", "Ka"]

## dihedrals = pd.read_csv('dihedrals', sep="\s+", header=None)
## dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
## dihedrals = dihedrals.replace(np.nan, '', regex=True)   # ci sono dei NaN che se tolti vengono messi i numeri con int normali

## dict_atomtypes = atoms.set_index("; nr")["type"].to_dict() # associazione di ["type"] a ("nr"). Per ogni nr c'e un type
                                                           # farei atom+resnr coso del percento 11
                                                           # contatore delle chains per cancellare atoms dal top e 
                                                           # per selezionare progressivamente i bonds e altro

## bonds[";ai"].replace(dict_atomtypes, inplace=True)
## bonds["aj"].replace(dict_atomtypes, inplace=True)
## bonds.to_string(index=False)
## kb_notation = bonds["kb"].map(lambda x: '{:.9e}'.format(x))
## r0_notation = bonds["r0"].map(lambda x: '{:.9e}'.format(x))
## bonds = bonds.assign(kb=kb_notation)
## bonds = bonds.assign(r0=r0_notation)

## angles[";ai"].replace(dict_atomtypes, inplace=True)
## angles["aj"].replace(dict_atomtypes, inplace=True)
## angles["ak"].replace(dict_atomtypes, inplace=True)
## angles.to_string(index=False)
## th0_notation = angles["th0"].map(lambda x: '{:.9e}'.format(x))
## ka_notation = angles["Ka"].map(lambda x: '{:.9e}'.format(x))
## angles = angles.assign(th0=th0_notation)
## angles = angles.assign(Ka=ka_notation)

## dihedrals[";ai"].replace(dict_atomtypes, inplace=True)
## dihedrals["aj"].replace(dict_atomtypes, inplace=True)
## dihedrals["ak"].replace(dict_atomtypes, inplace=True)
## dihedrals["al"].replace(dict_atomtypes, inplace=True)
## dihedrals.to_string(index=False)
## phi0_notation = dihedrals["phi0"].map(lambda x: '{:.9e}'.format(x))
## kd_notation = dihedrals["Kd"].map(lambda x: '{:.9e}'.format(x))
## dihedrals = dihedrals.assign(phi0=phi0_notation)
## dihedrals = dihedrals.assign(Kd=kd_notation)
## dihedrals["func"] = dihedrals["func"].replace(1, 9)

## bonds.columns = ["; i", "j", "func", "b0", "kb"]   #nuovi nomi delle colonne per usare direttamente il file come FF
#bonds.align = '1'              # questi due li ho messi non ricordo perche ma a quanto pare non servono. 
#bonds.border = False
## angles.columns = [";    i", "j", "k", "func", "th0", "Ka"]
## dihedrals.columns = [";  i", "j", "k", "l", "func", "phi", "kd", "mult"]

## file = open("ffbonded.itp", "w")
## file.write("[ bondtypes ]\n")
## file.write(str(bonds.to_string(index=False)))
## file.write("\n")
## file.write("\n")
## file.write("[ angletypes ]\n")
## file.write(str(angles.to_string(index=False)))
## file.write("\n")
## file.write("\n")
## file.write("[ dihedraltypes ]\n")
## file.write(str(dihedrals.to_string(index=False)))
## file.close()

## atomtypes = pd.DataFrame(atoms, columns = ['type'])
## atomtypes.columns = ["; type"]
## atomtypes["mass"] = '1.0000'
## atomtypes.insert(1, 'at.num', 1)
## atomtypes["charge"] = '0.000000'
## atomtypes["ptype"] = 'A'
## atomtypes["sigma"] = '0.00000e+00'
## atomtypes["epsilon"] = '5.96046e-09'

## pairs = pd.read_csv('pairs', sep="\s+", header=None)
## pairs.columns = [";ai", "aj", "type", "A", "B"]

## pairs[";ai"].replace(dict_atomtypes, inplace=True)
## pairs["aj"].replace(dict_atomtypes, inplace=True)
## pairs.to_string(index=False)
## A_notation = pairs["A"].map(lambda x: '{:.9e}'.format(x))
## B_notation = pairs["B"].map(lambda x: '{:.9e}'.format(x))
## pairs = pairs.assign(A=A_notation)
## pairs = pairs.assign(B=B_notation)
## pairs.columns = ["ai", "aj", "type", "A", "B"]

#print(pairs)

#file = open("temp_pairs", "w")
#file.write(str(pairs.to_string(index=False)))
#file.close()

# make a copy of the pairs with ai and aj inverted
## inv_pairs = pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
## inv_pairs.columns = ["ai", "aj", "type", "A", "B"]


#print(inv_pairs)

#file = open("inv_pairs", "w")
#file.write(str(inv_pairs.to_string(index=False)))
#file.close()


# append inv_pairs below pairs
## pairs_full = pairs.append(inv_pairs, sort=False, ignore_index=True)

#print(pairs_full)

#file = open("pairs_full", "w")
#file.write(str(pairs_full.to_string(index=False)))
#file.close()


# select the number of the residue and copy into two new rows
## n_ai = pairs_full.ai.str.extract('(\d+)')
## n_aj = pairs_full.aj.str.extract('(\d+)')

## pairs_full['n_ai'] = n_ai
## pairs_full['n_aj'] = n_aj

#print(pairs_full)


#file = open("pairs_full_aa", "w")
#file.write(str(pairs_full.to_string(index=False)))
#file.close()

# cleaning the inverted duplicates
# if ai is greater than aj then put in cond ai or else nan
## pairs_full['cond'] = np.where((pairs_full['ai'] >= pairs_full['aj']), pairs_full['ai'], np.nan)
## pairs_full = pairs_full.dropna()

#print(pairs_full)

#file = open("pairs_full_nan", "w")
#file.write(str(pairs_full.to_string(index=False)))
#file.close()

# deleting extra columns
## pairs_full = pairs_full.drop(["cond", "n_ai", "n_aj"], axis=1)

#print(pairs_full)

# sort based on ai
## pairs_full.sort_values(by=['ai', 'aj', 'A'], inplace=True)

#print(pairs_full)

# merge columns in order to drop the duplicates
## pairs_temp = pairs_full.copy()
## pairs_temp["ai"] = pairs_full["ai"].apply(str) + ':' + pairs_full["aj"].apply(str)

#print(pairs_temp)

# clean the remaining duplicates
## temp_clean = pairs_temp.drop_duplicates(subset='ai', keep="first")

#print(temp_clean)

#file = open("temp_clean", "w")
#file.write(str(temp_clean.to_string(index=False)))
#file.close()


# clean duplicates also on atomtypes
## atomtypes = atomtypes.drop_duplicates(subset='; type', keep="first")


# separate columns
## ai_aj = temp_clean["ai"].str.split(":", n = 1, expand = True)
## temp_clean["ai"] = ai_aj[0]
## temp_clean["aj"] = ai_aj[1]
## temp_clean.columns = [";ai", "aj", "type", "A", "B"]

## file = open("ffnonbonded.itp", "w")
## file.write("[ atomtypes ]\n")
## file.write(str(atomtypes.to_string(index=False)))
## file.write("\n")
## file.write("\n")
## file.write("[ nonbond_params ]\n")
## file.write(str(temp_clean.to_string(index=False)))
## file.close()