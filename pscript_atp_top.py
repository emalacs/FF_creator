#! /usr/bin/python

import numpy as np
import pandas as pd
import prettytable as pt
import copy

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

#changing the atom type column -> patoms column 5 + 3 

## atoms["type"] = atoms["atom"].apply(str) + '_' + atoms["resnr"].apply(str)

## file = open("topology_atoms", "w")
## file.write("[ atoms ]")
## file.write("\n")
## file.write(str(atoms.to_string(index=False)))
## file.close()

## atomtypes = pd.DataFrame(atoms, columns = ['type'])
## atomtypes.columns = ["; type"]
## atomtypes["mass"] = '1.0000'
##file = open("atomtypes.atp", "w")
##file.write("[ atomtypes ]")
##file.write("\n")
##file.write(str(atomtypes.to_string(index=False, header=False)))
##file.close()

## bonds = pd.read_csv('bonds', sep="\s+", header=None)    # stessa cosa di prima, bisogna togliere il nome delle colonne e reinserirle
## bonds.columns = [";ai", "aj", "func", "r0", "kb"]       # in attesa di un metodo piu intelligente

## angles = pd.read_csv('angles', sep="\s+", header=None)
## angles.columns = [";ai", "aj", "ak", "func", "th0", "Ka"]

## dihedrals = pd.read_csv('dihedrals', sep="\s+", header=None)
## dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
## dihedrals = dihedrals.replace(np.nan, '', regex=True)   # ci sono dei NaN che se tolti vengono messi i numeri con int normali

## atomtypes = pd.DataFrame(atoms, columns = ['type'])

## top_bonds = bonds.copy()
## top_angles = angles.copy()
## top_dihedrals = dihedrals.copy()

## top_bonds["r0"] = ""
## top_bonds["kb"] = ""
## top_angles["th0"] = ""
## top_angles["Ka"] = ""
## top_dihedrals["phi0"] = ""
## top_dihedrals["Kd"] = ""
## top_dihedrals["mult"] = ""

#file = open("topology_bonds", "w")
#file.write("[ bonds ]")
#file.write("\n")
#file.write(str(top_bonds.to_string(index=False)))
#file.close()

#file = open("topology_angles", "w")
#file.write("[ angles ]")
#file.write("\n")
#file.write(str(top_angles.to_string(index=False)))
#file.close()

#file = open("topology_dihedrals", "w")
#file.write("[ dihedrals ]")
#file.write("\n")
#file.write(str(top_dihedrals.to_string(index=False)))
#file.close()
