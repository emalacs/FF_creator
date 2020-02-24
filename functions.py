import numpy as np
import pandas as pd
from atomtypes_aa_definitions import *


# Dictionary definitions for both peptide and fibril

# atomtypes.atp created using the peptide and used also with the peptide

def make_atomtypes_and_dict(atomtypes):  # qui si mette l'output di read_*_atoms
    # print ("[ atomtypes ] of ffnonbonded")
    # This function prepare the file for ffnonbonded of the peptide
    # Insertion of information in atomtypes
    dict_atomtypes = atomtypes.set_index("; nr")["type"].to_dict()
    atomtypes['at.group'] = atomtypes['residue'] + '_' + atomtypes['atom']
    atomtypes['at.group'].replace(gromos_aa, inplace = True)
    # E QUI CI SIAMO
    atomtypes.insert(3, 'at.num', 4)
    atomtypes['at.num'] = atomtypes['at.group'].map(gromos_atp['at.num'])
    atomtypes.insert(4, 'mass', 5)
    atomtypes['mass'] = atomtypes['at.group'].map(gromos_atp['mass'])
    atomtypes["charge"] = '0.000000'
    atomtypes.insert(9, 'ptype', 10)
    atomtypes["ptype"] = 'A'
    atomtypes['c6'] = '0.00000e+00'
    atomtypes['c12'] = atomtypes['at.group'].map(gromos_atp['c12'])
    c12_notation = atomtypes["c12"].map(lambda x:'{:.6e}'.format(x))
    atomtypes = atomtypes.assign(c12 = c12_notation)
    atomtypes.drop(columns = ['; nr', 'resnr', 'residue', 'atom', 'cgnr', 'at.group'], inplace = True)
    atomtypes.rename(columns = {'type':'; type'}, inplace = True)
    # As expected this drop duplicates does not bother the peptide FF
    # Anche se a dire il vero gia ho fatto un drop duplicates
    # OCCHIO QUI

    atomtypes = atomtypes.drop_duplicates(subset = '; type', keep = 'first')

    # This part returns the atomtypes.atp
    atp = pd.DataFrame(atomtypes, columns = ['; type', 'mass'])
    return atp, atomtypes, dict_atomtypes


# Function for topology

def make_topology_bonds(top_bonds):
    # Emptying those two columns which are not needed into the topology
    top_bonds["r0"] = ""
    top_bonds["kb"] = ""
    # IN CASO CI FOSSERO DUPLICATI
    #  top_dihedrals = top_dihedrals.drop_duplicates(subset = [';ai', 'aj', 'ak', 'al'], keep = 'first')
    return top_bonds


def make_topology_angles(top_angles):
    # Emptying those two columns which are not needed into the topology
    top_angles["th0"] = ""
    top_angles["Ka"] = ""
    # IN CASO CI FOSSERO DUPLICATI
    #  top_dihedrals = top_dihedrals.drop_duplicates(subset = [';ai', 'aj', 'ak', 'al'], keep = 'first')
    return top_angles


def make_topology_dihedrals(top_dihedrals):
    # Emptying those three columns which are not needed into the topology
    top_dihedrals["phi0"] = ""
    top_dihedrals["Kd"] = ""
    top_dihedrals["mult"] = ""
    # the double dihedrals MUST be removed also in the topology
    top_dihedrals = top_dihedrals.drop_duplicates(subset = [';ai', 'aj', 'ak', 'al'], keep = 'first')
    return top_dihedrals


# Functions to prepare the bonds, angles and dihedrals for the ffbonded.itp creation

def ffbonded_bonds(bonds, dict_atomtypes):
    # Changing the atomnumber with the atomtype defined in the dictionary
    bonds[";ai"].replace(dict_atomtypes, inplace = True)
    bonds["aj"].replace(dict_atomtypes, inplace = True)
    # Cose per farlo funzionare
    bonds.to_string(index = False)
    #print(f'bonds\n{bonds}')
    bonds['kb'] = bonds['kb'] * (300 / 70)
    #print(f'bonds\n{bonds}')
    # Separazione delle colonne con dei numeri per la notazione scientifica
    kb_notation = bonds["kb"].map(lambda x:'{:.9e}'.format(x))
    r0_notation = bonds["r0"].map(lambda x:'{:.9e}'.format(x))
    # Sostituzione delle colonne all'interno del dataframe
    bonds = bonds.assign(kb = kb_notation)
    bonds = bonds.assign(r0 = r0_notation)
    # Ennesima definizione delle colonne
    bonds.columns = ["; i", "j", "func", "b0", "kb"]
    return bonds


def ffbonded_angles(angles, dict_atomtypes):
    # Changing the atomnumber with the atomtype defined in the dictionary
    angles[";ai"].replace(dict_atomtypes, inplace = True)
    angles["aj"].replace(dict_atomtypes, inplace = True)
    angles["ak"].replace(dict_atomtypes, inplace = True)
    # Cose per farlo funzionare
    angles.to_string(index = False)
    #print(f'angles\n{angles}')
    angles['Ka'] = angles['Ka'] * (300 / 70)
    #print(f'angles\n{angles}')
    # Separazione delle colonne con dei numeri per la notazione scientifica
    th0_notation = angles["th0"].map(lambda x:'{:.9e}'.format(x))
    ka_notation = angles["Ka"].map(lambda x:'{:.9e}'.format(x))
    # Sostituzione delle colonne all'interno del dataframe
    angles = angles.assign(th0 = th0_notation)
    angles = angles.assign(Ka = ka_notation)
    # Ennesima definizione delle colonne
    angles.columns = [";    i", "j", "k", "func", "th0", "Ka"]
    return angles


def ffbonded_dihedrals(dihedrals, dict_atomtypes):
    # Changing the atomnumber with the atomtype defined in the dictionary
    dihedrals[";ai"].replace(dict_atomtypes, inplace = True)
    dihedrals["aj"].replace(dict_atomtypes, inplace = True)
    dihedrals["ak"].replace(dict_atomtypes, inplace = True)
    dihedrals["al"].replace(dict_atomtypes, inplace = True)
    # Cose per farlo funzionare
    dihedrals.to_string(index = False)
    dihedrals['Kd'] = dihedrals['Kd'] * (300 / 70)
    # Separazione delle colonne con dei numeri per la notazione scientifica
    phi0_notation = dihedrals["phi0"].map(lambda x:'{:.9e}'.format(x))
    kd_notation = dihedrals["Kd"].map(lambda x:'{:.9e}'.format(x))
    # Sostituzione delle colonne all'interno del dataframe
    dihedrals = dihedrals.assign(phi0 = phi0_notation)
    dihedrals = dihedrals.assign(Kd = kd_notation)
    # Nei diedri e necessario sostituire l'1 con il 9 perche alcuni diedri vengono definiti due volte
    dihedrals["func"] = dihedrals["func"].replace(1, 9)
    # Ennesima definizione delle colonne
    dihedrals.columns = [";  i", "j", "k", "l", "func", "phi", "kd", "mult"]
    return dihedrals


# This is for the merge peptide + fibril, which is quite similar but there are a few lines more


def ffbonded_merge_dihedrals(pep_dihedrals, fib_dihedrals, dict_pep_atomtypes, dict_fib_atomtypes):
    # This function is used to calculate dihedrals when merging the FF of peptide and fibril
    # Changing the atomnumber with the atomtype defined in the dictionary
    # Peptide input handling
    pep_dihedrals[";ai"].replace(dict_pep_atomtypes, inplace = True)
    pep_dihedrals["aj"].replace(dict_pep_atomtypes, inplace = True)
    pep_dihedrals["ak"].replace(dict_pep_atomtypes, inplace = True)
    pep_dihedrals["al"].replace(dict_pep_atomtypes, inplace = True)
    # Fibril input handling

    fib_dihedrals[";ai"].replace(dict_fib_atomtypes, inplace = True)
    fib_dihedrals["aj"].replace(dict_fib_atomtypes, inplace = True)
    fib_dihedrals["ak"].replace(dict_fib_atomtypes, inplace = True)
    fib_dihedrals["al"].replace(dict_fib_atomtypes, inplace = True)

    # Here starts the difference between the normal dihedral and the merge one
    # Some dihedrals will be duplicates for peptide and fibril and we want to keep both
    # Discrimination of type 9 (which now are 1) and the type 2
    # QUESTO PASSAGGIO TRIGGERA IL COPY WARNING MA AL MOMENTO E COSI CHE VOGLIO FARE
    pep_kd_1 = pep_dihedrals.loc[pep_dihedrals['func'] == 1]
    pep_kd_2 = pep_dihedrals.loc[pep_dihedrals['func'] == 2]
    fib_kd_1 = fib_dihedrals.loc[fib_dihedrals['func'] == 1]
    # Type 9 (1) of peptide and fibril will be divided by the half
    # and the type 2 will be kept only the one from the peptide
    # QUESTI DUE COMANDI RISOLVONO IL COPY WARNING
    pep_kd_1.is_copy = False
    fib_kd_1.is_copy = False
    pep_kd_1.loc[:, 'Kd'] = pep_kd_1.loc[:, 'Kd'].divide(2)
    fib_kd_1.loc[:, 'Kd'] = fib_kd_1.loc[:, 'Kd'].divide(2)
    # QUESTI DUE COMANDI SONO QUELLI ORIGINALI PRIMA DELL'OTTIMIZZAZIONE COPY WARNING
    # pep_kd_1['Kd'] = pep_kd_1['Kd'].divide(2)
    # fib_kd_1['Kd'] = fib_kd_1['Kd'].divide(2)
    # Having the half kd_1 of peptide and fibril appending peptide kd_1, kd_2 and fibril kd_1
    merge_pep_dihedrals = pep_kd_1.append(pep_kd_2, sort = False, ignore_index = True)
    merge_dihedrals = merge_pep_dihedrals.append(fib_kd_1, sort = False, ignore_index = True)
    # Here is the common part of the dihedrals for peptide and fibril
    merge_dihedrals.to_string(index = False)
    merge_dihedrals['Kd'] = merge_dihedrals['Kd'] * (300 / 70)
    phi0_notation = merge_dihedrals["phi0"].map(lambda x:'{:.9e}'.format(x))
    kd_notation = merge_dihedrals["Kd"].map(lambda x:'{:.9e}'.format(x))
    merge_dihedrals = merge_dihedrals.assign(phi0 = phi0_notation)
    merge_dihedrals = merge_dihedrals.assign(Kd = kd_notation)
    merge_dihedrals["func"] = merge_dihedrals["func"].replace(1, 9)
    # This step is required since gromacs wants the dihedrals sorted by type
    merge_dihedrals.sort_values(by = [';ai', 'aj', 'ak', 'al'], inplace = True)
    # Ennesima definizione delle colonne
    merge_dihedrals.columns = [";  i", "j", "k", "l", "func", "phi", "kd", "mult"]
    return merge_dihedrals


def ffnonbonded_pep_pairs(pep_pairs, dict_pep_atomtypes):
    pep_pairs[";ai"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs["aj"].replace(dict_pep_atomtypes, inplace = True)
    # Cose per farlo funzionare
    pep_pairs.to_string(index = False)
    # Ennesima definizione delle colonne
    pep_pairs.columns = [";ai", "aj", "type", "A", "B"]

    # This part allow to reweight all the C6 and C12 of the pairs. So it is possible for us to use higher temperatures
    # such as 300 K and have a proper viscosity during the simulation.
    # Basically, the intramolecular are ok at 70K but the molecules moves like they're in water at 70K.
    # print(pep_pairs)
    # pep_old_epsilon = (pep_pairs['A'] ** 2) / (4 * (pep_pairs['B']))
    # print(f'\tc6^2/(4xc12)'
    #      f'\n'
    #      f'\tOld Peptide Epsilon = {pep_old_epsilon.iloc[0]}\n')
    # t_check = pep_old_epsilon / (8.314462618 * (10 ** -3))
    # print(f'\toldEpsilon/R in KJ'
    #      f'\n'
    #      f'\tTemperature = {t_check.iloc[0]}\n')
    # pep_new_epsilon = pep_old_epsilon * (300 / t_check)  # qui sarebbe carino mettere un input della temperatura
    # print(f'\toldEpsilon/(300/{t_check.iloc[0]})'
    #      f'\n'
    #      f'\tNew Peptide Epsilon = {pep_new_epsilon.iloc[0]}\n')
    # In questo caso stai comunque moltiplicando per due i valori di C6 e C12 e dunque aumentano.
    # fib_pairs_full['A'] = fib_pairs_full['A'] * (300/t_check)
    # fib_pairs_full['B'] = fib_pairs_full['B'] * (300/t_check)

    pep_pairs['A'] = pep_pairs['A'] * (300 / 70)
    pep_pairs['B'] = pep_pairs['B'] * (300 / 70)
    # print(fib_pairs_full)

    # Notazione scientifica
    A_notation = pep_pairs["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = pep_pairs["B"].map(lambda x:'{:.9e}'.format(x))
    # Sostituzione delle colonne all'interno del dataframe
    pep_pairs = pep_pairs.assign(A = A_notation)
    pep_pairs = pep_pairs.assign(B = B_notation)
    # print(pep_pairs)
    return pep_pairs


def ffnonbonded_fib_pairs(fib_pairs, dict_atomtypes):
    fib_pairs[";ai"].replace(dict_atomtypes, inplace = True)
    fib_pairs["aj"].replace(dict_atomtypes, inplace = True)
    # Cose per farlo funzionare
    fib_pairs.to_string(index = False)

    # Ennesima definizione delle colonne
    fib_pairs.columns = ["ai", "aj", "type", "A", "B"]

    # Here there is the specific part for the fibril, to search and destroy duplicates
    # The strategy is based on a copy of the appending of a copy of the pairs with inverted ai and aj
    inv_fib_pairs = fib_pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
    inv_fib_pairs.columns = ['ai', 'aj', 'type', 'A', 'B']
    fib_pairs_full = fib_pairs.append(inv_fib_pairs, sort = False, ignore_index = True)
    # print(fib_pairs_full)
    # Then, the residue number is selected and copied into two new rows
    n_ai = fib_pairs_full.ai.str.extract('(\d+)')
    n_aj = fib_pairs_full.aj.str.extract('(\d+)')
    fib_pairs_full['n_ai'] = n_ai
    fib_pairs_full['n_aj'] = n_aj
    # This third step cleans the inverted duplicates
    # If ai is greater than aj then put in cond ai or else nan
    fib_pairs_full['cond'] = np.where((fib_pairs_full['n_ai'] >= fib_pairs_full['n_aj']), fib_pairs_full['ai'], np.nan)
    fib_pairs_full = fib_pairs_full.dropna()
    # Finally the deletion of the extra columns
    fib_pairs_full = fib_pairs_full.drop(['cond', 'n_ai', 'n_aj'], axis = 1)
    # Sorting the pairs
    fib_pairs_full.sort_values(by = ['ai', 'aj', 'A'], inplace = True)
    # Duplicates removal
    # Merging columns in order to drop the duplicates
    # fib_pairs_full['ai'] = fib_pairs_full['ai'].apply(str) + ':' + fib_pairs_full['aj'].apply(str) # questo come
    #                                                                                                  funzionava prima
    # Cleaning the remaining duplicates
    fib_pairs_full = fib_pairs_full.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Column separation
    ai_aj = fib_pairs_full['ai'].str.split(":", n = 1, expand = True)
    fib_pairs_full.loc[:, 'ai'] = ai_aj.loc[:, 0]
    fib_pairs_full.columns = [';ai', 'aj', 'type', 'A', 'B']

    # This part allow to reweight all the C6 and C12 of the pairs. So it is possible for us to use higher temperatures
    # such as 300 K and have a proper viscosity during the simulation.
    # Basically, the intramolecular are ok at 70K but the molecules moves like they're in water at 70K.
    # print(fib_pairs_full)
    # fib_old_epsilon = (fib_pairs_full['A'] ** 2) / (4 * (fib_pairs_full['B']))
    # print(f'\tc6^2/(4xc12)'
    #      f'\n'
    #      f'\tOld Fibril Epsilon = {fib_old_epsilon.iloc[0]}\n')
    # t_check = fib_old_epsilon / (8.314462618 * (10 ** -3))
    # print(f'\toldEpsilon/R in KJ'
    #      f'\n'
    #      f'\tTemperature = {t_check.iloc[0]}\n')
    # fib_new_epsilon = fib_old_epsilon * (300/t_check) # qui sarebbe carino mettere un input della temperatura
    # print(f'\toldEpsilon/(300/{t_check.iloc[0]})'
    #      f'\n'
    #      f'\tNew Fibril Epsilon = {fib_new_epsilon.iloc[0]}\n')
    # In questo caso stai comunque moltiplicando per due i valori di C6 e C12 e dunque aumentano.
    # fib_pairs_full['A'] = fib_pairs_full['A'] * (300/t_check)
    # fib_pairs_full['B'] = fib_pairs_full['B'] * (300/t_check)

    fib_pairs_full['A'] = fib_pairs_full['A'] * (300 / 70)
    fib_pairs_full['B'] = fib_pairs_full['B'] * (300 / 70)
    # print(fib_pairs_full)

    # Notazione scientifica
    A_notation = fib_pairs_full["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = fib_pairs_full["B"].map(lambda x:'{:.9e}'.format(x))
    # Sostituzione delle colonne all'interno del dataframe
    fib_pairs_full = fib_pairs_full.assign(A = A_notation)
    fib_pairs_full = fib_pairs_full.assign(B = B_notation)
    # print(fib_pairs_full)
    return fib_pairs_full


# This function is similar to the previous one
# however it appends also the values for the merged ff pairs

def ffnonbonded_merge_pairs(pep_pairs, fib_pairs, dict_pep_atomtypes, dict_fib_atomtypes):
    # This script allow to merge the pairs of peptide and fibril.
    # The main difference between the other two pairs function is that peptide C6 and C12 are reweighted.
    # This is because SMOG normalize the LJ potential based on the total number of contacts.
    # Since the peptide has less contacts, the LJ potential is stronger than the fibril
    # Peptide input handling
    pep_pairs[";ai"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs["aj"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs.to_string(index = False)
    pep_pairs.columns = ["ai", "aj", "type", "A", "B"]

    pep_pairs['A'] = pep_pairs['A'] * (300 / 70)
    pep_pairs['B'] = pep_pairs['B'] * (300 / 70)

    # Fibril input handling
    fib_pairs['ai'].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs["aj"].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs.to_string(index = False)
    fib_pairs.columns = ["ai", "aj", "type", "A", "B"]

    fib_pairs['A'] = fib_pairs['A'] * (300 / 70)
    fib_pairs['B'] = fib_pairs['B'] * (300 / 70)

    # Calcolo di epsilon per peptide e fibrilla
    pep_epsilon = (pep_pairs['A'] ** 2) / (4 * (pep_pairs['B']))
    fib_epsilon = (fib_pairs['A'] ** 2) / (4 * (fib_pairs['B']))
    ratio = pep_epsilon[0] / fib_epsilon[0]
    # QUESTO PRINT CI PIACE MOLTO MA E' SOLO PER PYTHON 3
    print(f'\n'
          f'\tPeptide epsilon: {pep_epsilon[0]}\n'
          f'\tFibril epsilon: {fib_epsilon[0]}\n'
          f'\tRatio: {ratio}'
          f'\n')

    # Reweight peptide LJ
    pep_pairs['A'] = pep_pairs['A'] / ratio
    pep_pairs['B'] = pep_pairs['B'] / ratio

    # From now the function behaves like the others
    A_notation = pep_pairs["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = pep_pairs["B"].map(lambda x:'{:.9e}'.format(x))
    pep_pairs = pep_pairs.assign(A = A_notation)
    pep_pairs = pep_pairs.assign(B = B_notation)
    A_notation = fib_pairs["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = fib_pairs["B"].map(lambda x:'{:.9e}'.format(x))
    fib_pairs = fib_pairs.assign(A = A_notation)
    fib_pairs = fib_pairs.assign(B = B_notation)

    # One last step about merging the pairs
    pairs = pep_pairs.append(fib_pairs, sort = False, ignore_index = True)

    # Cleaning the duplicates (the logic has already been explained above
    inv_pairs = pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
    inv_pairs.columns = ['ai', 'aj', 'type', 'A', 'B']
    pairs_full = pairs.append(inv_pairs, sort = False, ignore_index = True)
    n_ai = pairs_full.ai.str.extract('(\d+)')
    n_aj = pairs_full.aj.str.extract('(\d+)')
    pairs_full['n_ai'] = n_ai
    pairs_full['n_aj'] = n_aj
    pairs_full['cond'] = np.where((pairs_full['n_ai'] >= pairs_full['n_aj']), pairs_full['ai'], np.nan)
    pairs_full = pairs_full.dropna()
    pairs_full = pairs_full.drop(['cond', 'n_ai', 'n_aj'], axis = 1)
    # Sorting the pairs
    pairs_full.sort_values(by = ['ai', 'aj', 'A'], inplace = True)
    # Duplicates removal
    # Merging columns in order to drop the duplicates
    # pairs_full['ai'] = pairs_full['ai'].apply(str) + ':' + pairs_full['aj'].apply(str) # questo come funzionava prima
    # Cleaning the remaining duplicates
    pairs_full = pairs_full.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Column separation
    ai_aj = pairs_full['ai'].str.split(":", n = 1, expand = True)
    # QUESTI DUE COMANDI SONO DOPO IL SET COPY WARNING
    pairs_full.loc[:, 'ai'] = ai_aj.loc[:, 0]

    # QUESTI DUE COMANDI SONO QUELLI PRIMA DEL COPY WARNING
    # clean_pairs['ai'] = ai_aj[0]
    # clean_pairs['aj'] = ai_aj[1]
    pairs_full.columns = [';ai', 'aj', 'type', 'A', 'B']
    return pairs_full


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# This function is a basic pairs analyzing tool which can come handy and to develop in the future

def pairs_separator(pairs):
    # Poi commenta un po' quello che serve che neanche ora ti ricordi come funziona questo script

    ai_in = pairs['ai'].str.split('_', n = 1, expand = True)
    aj_jn = pairs['aj'].str.split('_', n = 1, expand = True)
    pairs['ain'] = ai_in[1]
    pairs['air'] = ai_in[0]
    pairs['ajn'] = aj_jn[1]
    pairs['ajr'] = aj_jn[0]
    pairs.astype(str)

    # Here we start to select the backbone and the sidechains
    # The mixed pairs between the two are not included atm
    back_atomtype = ['CA', 'N', 'O', 'C', 'OXT']
    # Backbone
    drop_sidechains = pairs[pairs['air'].isin(back_atomtype)]
    backbone = drop_sidechains[drop_sidechains['ajr'].isin(back_atomtype)]
    # Sidechains
    drop_backbone = pairs[~ pairs['air'].isin(back_atomtype)]
    sidechains = drop_backbone[~ drop_backbone['ajr'].isin(back_atomtype)]
    # mixed #con questo metodo non funziona
    mix_air = pairs[pairs['air'].isin(back_atomtype)]
    mix_ajr = pairs[pairs['ajr'].isin(back_atomtype)]
    mix_full = mix_air.append(mix_ajr, sort = False, ignore_index = True)
    # mixed poi ci penso
    pymol_backbone = backbone.copy()
    pymol_sidechains = sidechains.copy()
    pymol_mixed = mix_full.copy()

    # cleaning the dataframes
    backbone = backbone.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)
    sidechains = sidechains.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)
    mix_full = mix_full.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)

    file = open("pymol/backbone", "w")
    file.write(str(backbone.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/sidechains", "w")
    file.write(str(sidechains.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/mixed", "w")
    file.write(str(mix_full.to_string(index = False, header = False)))
    file.close()

    # backbone
    pymol_backbone = pymol_backbone.drop(['type', 'A', 'B'], axis = 1)
    pymol_backbone['interaction'] = np.where((pymol_backbone['ai'] == pymol_backbone['aj']), 'v', 'h')
    # PROVA AD UNTILIZZARE .JOIN
    name = pymol_backbone['interaction'].apply(str) + '_' + pymol_backbone['ai'].apply(str) + ':' + pymol_backbone[
        'aj'].apply(str)
    pymol_backbone = pymol_backbone.drop(['ai', 'aj'], axis = 1)
    pymol_backbone.insert(0, 'distance', 'distance')
    pymol_backbone.insert(1, 'name', name)
    pymol_backbone.insert(2, 'resi', ', resi')
    pymol_backbone.insert(4, 'and', 'and name')
    pymol_backbone.insert(6, 'resi2', ', resi')
    pymol_backbone.insert(8, 'and2', 'and name')
    pymol_backbone.insert(10, 'cutoff', ', 6')
    pymol_backbone = pymol_backbone.sort_values(by = 'interaction')
    pymol_backbone = pymol_backbone.drop(['interaction'], axis = 1)

    # sidechains
    pymol_sidechains = pymol_sidechains.drop(['type', 'A', 'B'], axis = 1)
    pymol_sidechains['interaction'] = np.where((pymol_sidechains['ai'] == pymol_sidechains['aj']), 'v', 'h')
    name = pymol_sidechains['interaction'].apply(str) + '_' + pymol_sidechains['ai'].apply(str) + ':' + \
           pymol_sidechains['aj'].apply(str)
    pymol_sidechains = pymol_sidechains.drop(['ai', 'aj'], axis = 1)
    pymol_sidechains.insert(0, 'distance', 'distance')
    pymol_sidechains.insert(1, 'name', name)
    pymol_sidechains.insert(2, 'resi', ', resi')
    pymol_sidechains.insert(4, 'and', 'and name')
    pymol_sidechains.insert(6, 'resi2', ', resi')
    pymol_sidechains.insert(8, 'and2', 'and name')
    pymol_sidechains.insert(10, 'cutoff', ', 6')
    pymol_sidechains = pymol_sidechains.sort_values(by = 'interaction')
    pymol_sidechains = pymol_sidechains.drop(['interaction'], axis = 1)

    # mixed
    pymol_mixed = pymol_mixed.drop(['type', 'A', 'B'], axis = 1)
    pymol_mixed['interaction'] = np.where((pymol_mixed['ai'] == pymol_mixed['aj']), 'v', 'h')
    name = pymol_mixed['interaction'].apply(str) + '_' + pymol_mixed['ai'].apply(str) + ':' + pymol_mixed['aj'].apply(
        str)
    pymol_mixed = pymol_mixed.drop(['ai', 'aj'], axis = 1)
    pymol_mixed.insert(0, 'distance', 'distance')
    pymol_mixed.insert(1, 'name', name)
    pymol_mixed.insert(2, 'resi', ', resi')
    pymol_mixed.insert(4, 'and', 'and name')
    pymol_mixed.insert(6, 'resi2', ', resi')
    pymol_mixed.insert(8, 'and2', 'and name')
    pymol_mixed.insert(10, 'cutoff', ', 6')
    pymol_mixed = pymol_mixed.sort_values(by = 'interaction')
    pymol_mixed = pymol_mixed.drop(['interaction'], axis = 1)

    file = open("pymol/pymol_backbone", "w")
    file.write(str(pymol_backbone.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/pymol_sidechains", "w")
    file.write(str(pymol_sidechains.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/pymol_mixed", "w")
    file.write(str(pymol_mixed.to_string(index = False, header = False)))
    file.close()
