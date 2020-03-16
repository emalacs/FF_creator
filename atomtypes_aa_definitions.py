import pandas as pd
from read_input_files import *

# Script containing the definitions of the atomtypes of every aminoacid used in here. Taken from gromos/aminoacids.rtp
# along with the mass and c12 for every atomtype contained in Gromos

# Gromos atomtypes definition
# Abbiamo tenuto solo le info che ci servivano. Da aggiungere le definizioni dei terminali e quelli che serviranno
# per i prossimi aminoacidi. Abbiamo aggiunto gli H nella massa (e nel c12) - OA ed N -
# per aggiungere C12 "regola del 20".

# GROMOS 53a6

# Gromos dictionaries defintions of:
    # dict_gro_atomtypes -> number:residue_atom 1:TYR_N
    # gromos_mass -> atom_type:mass N:150067
    # res_atom_dict -> residue_atom:chemical type TYR_N:N
gro_atoms = read_gro_atoms()
dict_gro_atomtypes = gro_atoms.set_index('; nr')['res_atom'].to_dict()

gromos_mass = gro_atoms[['type', 'mass']].drop_duplicates(subset=['type'], keep='first').copy()
gromos_mass_dict = gromos_mass.set_index('type')['mass'].to_dict() # This one will be used in make_atomtypes

print('RICORDATI DI CAMBIARE LE MASSE NEL FF, QUI GLI H SONO TUTTI ESPLICITI!!!!')
gromos_res_atom_dict = gro_atoms.set_index('res_atom')['type'].to_dict() # This one will be used in make_atomtypes
# This one will be used for proper dihedrals from SMOG to past to GROMOS topology
gromos_resatom_nmr = gro_atoms[['res_atom', 'resnr', '; nr']].copy()
gromos_resatom_nmr['resnr_2'] = gromos_resatom_nmr['resnr'].astype(str)
gromos_resatom_nmr['res_atom'] = gromos_resatom_nmr['res_atom'] + '_' + gromos_resatom_nmr['resnr_2']
gromos_resatom_nmr_dict = gromos_resatom_nmr.set_index('res_atom')['; nr'].to_dict()


print(gromos_resatom_nmr_dict)
#gromos_resatom_nmr_dict = gro_atoms.set_index('res_atom')['; nr']#.to_dict()


print(len(gromos_resatom_nmr_dict))

# Gromos impropers dictionary residue_atom to : define
gro_impropers = read_gro_impropers()
gro_impropers['; ai'].replace(dict_gro_atomtypes, inplace = True)
gro_impropers['aj'].replace(dict_gro_atomtypes, inplace = True)
gro_impropers['ak'].replace(dict_gro_atomtypes, inplace = True)
gro_impropers['al'].replace(dict_gro_atomtypes, inplace = True)
gro_impropers['dihedrals'] = gro_impropers['; ai'] + '+' + gro_impropers['aj'] + '+' +\
                             gro_impropers['ak'] + '+' + gro_impropers['al']
dict_gro_impropers = gro_impropers.set_index('dihedrals')['define'].to_dict()

# THE C12 RATIO CHANGED a little bit.
gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 'CH2', 'CH3', 'CH2r', 'NT'],
     'mass': [16, 17, 15, 12, 13, 14, 15, 14, 17],
     'at.num': [8, 8, 7, 6, 6, 6, 6, 6, 7],
     'c12': [1e-06, 3.011e-05, 4.639e-05, 4.937284e-06, 9.70225e-05, 3.3965584e-05,
             2.6646244e-05, 2.8058209e-05, 1.2e-05]
     }
)

gromos_atp.to_dict()
gromos_atp.set_index('name', inplace=True)

# Bonds gromos definition
# What I quicly saw is that gromos can have up to 100 fold increase of the bond strength.
gromos_bonds = pd.DataFrame(
    {'def': ['gb_5', 'gb_10', 'gb_13', 'gb_16', 'gb_18', 'gb_21', 'gb_27'],
     'force': ['1.6600e+07', '1.1800e+07', '1.0200e+07', '1.0800e+07', '8.1800e+06', '8.7100e+06', '7.1500+e06'],
     'len': ['0.1230', '0.1330', '0.1360', '0.1390', '0.1430', '0.1470', '0.1530']
     }
)

# Dictionaries used to replace the SMOG definitions during the bonds creation.
dict_gromos_bonds_force = gromos_bonds.set_index('def')['force'].to_dict()
dict_gromos_bonds_len = gromos_bonds.set_index('def')['len'].to_dict()

gromos_angles = pd.DataFrame(
    {'def': ['ga_13', 'ga_15', 'ga_19', 'ga_27', 'ga_30', 'ga_31', 'ga_33', 'ga_38', 'ga_22'],
     'force': ['520.00', '530.00', '610.00', '560.00', '685.00', '700.00', '730.00', '770.00', '635.00'],
     'angle': ['109.50', '111.00', '115.00', '120.00', '121.00', '122.00', '124.00', '126.00', '117.00']
     }
)

dict_gromos_angles_force = gromos_angles.set_index('def')['force'].to_dict()
dict_gromos_angles_angle = gromos_angles.set_index('def')['angle'].to_dict()

gromos_dihedrals = pd.DataFrame(
    {'def': [''],
     'force': [''],
     'len': ['']
     }
)

gromos_impropers = pd.DataFrame(
    {'def': ['gi_1', 'gi_2'],
     'force': ['167.42309', '334.84617'],
     'dihe': ['0.0', '35.26439']
     }
)

dict_gromos_impropers_force = gromos_impropers.set_index('def')['force'].to_dict()
dict_gromos_impropers_dihe = gromos_impropers.set_index('def')['dihe'].to_dict()

# TTR sequence YTIAALLSPYS

# Alanine bonds definition.
bALA = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

bALA['ai'] = 'ALA_' + bALA['ai'].astype(str)
bALA['aj'] = 'ALA_' + bALA['aj'].astype(str)
bALA.insert(2, 'bond', 3)
bALA['bond'] = bALA['ai'] + '+' + bALA['aj']

bALA_dict = bALA.set_index('bond')['def'].to_dict()

# Alananine angles definition.
aALA = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

# print(aALA)
aALA['ai'] = 'ALA_' + aALA['ai'].astype(str)
aALA['aj'] = 'ALA_' + aALA['aj'].astype(str)
aALA['ak'] = 'ALA_' + aALA['ak'].astype(str)
aALA.insert(3, 'angle', 4)
aALA['angle'] = aALA['ai'] + '+' + aALA['aj'] + '+' + aALA['ak']
# print(aALA)
aALA_dict = aALA.set_index('angle')['def'].to_dict()
# print(aALA_dict)

# Alanine proper dihedrals definition.
dALA = pd.DataFrame({'ai': ['CA', 'C', 'N'],
                     'aj': ['C', 'N', 'CA'],
                     'ak': ['N', 'CA', 'C'],
                     'al': ['CA', 'C', 'N'],
                     'def': ['gd_14', 'gd_39', 'gd_40']}
                    )

# Alananine improper dihedrals definition.
iALA = pd.DataFrame({'ai': ['CA', 'C'],
                     'aj': ['N', 'CA'],
                     'ak': ['C', 'N'],
                     'al': ['CB', 'O'],
                     'def': ['gi_2', 'gi_1']}
                    )

# print(iALA)
iALA['ai'] = 'ALA_' + iALA['ai'].astype(str)
iALA['aj'] = 'ALA_' + iALA['aj'].astype(str)
iALA['ak'] = 'ALA_' + iALA['ak'].astype(str)
iALA['al'] = 'ALA_' + iALA['al'].astype(str)
iALA.insert(4, 'improper', 5)
iALA['improper'] = iALA['ai'] + '+' + iALA['aj'] + '+' + iALA['ak'] + '+' + iALA['al']
# print(iALA)
iALA_dict = iALA.set_index('improper')['def'].to_dict()
# print(iALA_dict)

# Threonine definition.
bTHR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG1', 'CG2', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_27', 'gb_5', 'gb_10']}
                    )

# print(bTHR)
bTHR['ai'] = 'THR_' + bTHR['ai'].astype(str)
bTHR['aj'] = 'THR_' + bTHR['aj'].astype(str)
bTHR.insert(2, 'bond', 3)
bTHR['bond'] = bTHR['ai'] + '+' + bTHR['aj']
# print(bTHR)
bTHR_dict = bTHR.set_index('bond')['def'].to_dict()
# print(bTHR_dict)

aTHR = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CA', 'OG1', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'OG1', 'CG2', 'CG2', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_30', 'ga_19', 'ga_33']}
                    )

# print(aTHR)
aTHR['ai'] = 'THR_' + aTHR['ai'].astype(str)
aTHR['aj'] = 'THR_' + aTHR['aj'].astype(str)
aTHR['ak'] = 'THR_' + aTHR['ak'].astype(str)
aTHR.insert(3, 'angle', 4)
aTHR['angle'] = aTHR['ai'] + '+' + aTHR['aj'] + '+' + aTHR['ak']
# print(aTHR)
aTHR_dict = aTHR.set_index('angle')['def'].to_dict()
# print(aTHR_dict)

# dTHR = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iTHR = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'OG1', 'CA'],
                     'ak': ['C', 'CG2', 'N'],
                     'al': ['CB', 'CA', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

# print(iTHR)
iTHR['ai'] = 'THR_' + iTHR['ai'].astype(str)
iTHR['aj'] = 'THR_' + iTHR['aj'].astype(str)
iTHR['ak'] = 'THR_' + iTHR['ak'].astype(str)
iTHR['al'] = 'THR_' + iTHR['al'].astype(str)
iTHR.insert(4, 'improper', 5)
iTHR['improper'] = iTHR['ai'] + '+' + iTHR['aj'] + '+' + iTHR['ak'] + '+' + iTHR['al']
# print(iTHR)
iTHR_dict = iTHR.set_index('improper')['def'].to_dict()
# print(iTHR_dict)


# Tyrosine definition.
bTYR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'OH', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16',
                             'gb_13', 'gb_5', 'gb_10']}
                    )

# print(bTYR)
bTYR['ai'] = 'TYR_' + bTYR['ai'].astype(str)
bTYR['aj'] = 'TYR_' + bTYR['aj'].astype(str)
bTYR.insert(2, 'bond', 3)
bTYR['bond'] = bTYR['ai'] + '+' + bTYR['aj']
# print(bTYR)
bTYR_dict = bTYR.set_index('bond')['def'].to_dict()
# print(bTYR_dict)

aTYR = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CB', 'CB', 'CD1', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE1',
                            'CE2', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ',
                            'CZ', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'CG', 'CD1', 'CD2', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'CE2', 'OH',
                            'OH', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27',
                             'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_30', 'ga_19', 'ga_33']}
                    )

aTYR['ai'] = 'TYR_' + aTYR['ai'].astype(str)
aTYR['aj'] = 'TYR_' + aTYR['aj'].astype(str)
aTYR['ak'] = 'TYR_' + aTYR['ak'].astype(str)
aTYR.insert(3, 'angle', 4)
aTYR['angle'] = aTYR['ai'] + '+' + aTYR['aj'] + '+' + aTYR['ak']
aTYR_dict = aTYR.set_index('angle')['def'].to_dict()

# dTYR = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iTYR = pd.DataFrame({'ai': ['CA', 'CG', 'CG', 'CG', 'CD1', 'CD1', 'CD2', 'CD2', 'CZ', 'C'],
                     'aj': ['N', 'CD1', 'CD1', 'CD2', 'CG', 'CE1', 'CG', 'CE2', 'CE1', 'CA'],
                     'ak': ['C', 'CD2', 'CE1', 'CE2', 'CD2', 'CZ', 'CD1', 'CZ', 'CE2', 'N'],
                     'al': ['CB', 'CB', 'CZ', 'CZ', 'CE2', 'CE2', 'CE1', 'CE1', 'OH', 'O'],
                     'def': ['gi_2', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1']}
                    )

# print(iTYR)
iTYR['ai'] = 'TYR_' + iTYR['ai'].astype(str)
iTYR['aj'] = 'TYR_' + iTYR['aj'].astype(str)
iTYR['ak'] = 'TYR_' + iTYR['ak'].astype(str)
iTYR['al'] = 'TYR_' + iTYR['al'].astype(str)
iTYR.insert(4, 'improper', 5)
iTYR['improper'] = iTYR['ai'] + '+' + iTYR['aj'] + '+' + iTYR['ak'] + '+' + iTYR['al']
# print(iTYR)
iTYR_dict = iTYR.set_index('improper')['def'].to_dict()
# print(iTYR_dict)


# Isoleucine definition.
bILE = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'CG1', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG1', 'CG2', 'CD1', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

# print(bILE)
bILE['ai'] = 'ILE_' + bILE['ai'].astype(str)
bILE['aj'] = 'ILE_' + bILE['aj'].astype(str)
bILE.insert(2, 'bond', 3)
bILE['bond'] = bILE['ai'] + '+' + bILE['aj']
# print(bILE)
bILE_dict = bILE.set_index('bond')['def'].to_dict()
# print(bILE_dict)

aILE = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CA', 'CG1', 'CB', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'CG1', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'CG1', 'CG2', 'CG2', 'CD1', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_15', 'ga_15', 'ga_30',
                             'ga_19', 'ga_33']}
                    )

# print(aILE)
aILE['ai'] = 'ILE_' + aILE['ai'].astype(str)
aILE['aj'] = 'ILE_' + aILE['aj'].astype(str)
aILE['ak'] = 'ILE_' + aILE['ak'].astype(str)
aILE.insert(3, 'angle', 4)
aILE['angle'] = aILE['ai'] + '+' + aILE['aj'] + '+' + aILE['ak']
# print(aILE)
aILE_dict = aILE.set_index('angle')['def'].to_dict()
# print(aILE_dict)

# dILE = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iILE = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'CG1', 'CA'],
                     'ak': ['C', 'CG2', 'N'],
                     'al': ['CB', 'CA', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

# print(iILE)
iILE['ai'] = 'ILE_' + iILE['ai'].astype(str)
iILE['aj'] = 'ILE_' + iILE['aj'].astype(str)
iILE['ak'] = 'ILE_' + iILE['ak'].astype(str)
iILE['al'] = 'ILE_' + iILE['al'].astype(str)
iILE.insert(4, 'improper', 5)
iILE['improper'] = iILE['ai'] + '+' + iILE['aj'] + '+' + iILE['ak'] + '+' + iILE['al']
# print(iILE)
iILE_dict = iILE.set_index('improper')['def'].to_dict()
# print(iILE_dict)

# Leucine definition.
bLEU = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CG', 'CG', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG', 'CD1', 'CD2', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

bLEU['ai'] = 'LEU_' + bLEU['ai'].astype(str)
bLEU['aj'] = 'LEU_' + bLEU['aj'].astype(str)
bLEU.insert(2, 'bond', 3)
bLEU['bond'] = bLEU['ai'] + '+' + bLEU['aj']
bLEU_dict = bLEU.set_index('bond')['def'].to_dict()

aLEU = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CB', 'CB', 'CD1', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CG', 'CG', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'CG', 'CD1', 'CD2', 'CD2', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_15', 'ga_15', 'ga_30', 'ga_19',
                             'ga_33']}
                    )

aLEU['ai'] = 'LEU_' + aLEU['ai'].astype(str)
aLEU['aj'] = 'LEU_' + aLEU['aj'].astype(str)
aLEU['ak'] = 'LEU_' + aLEU['ak'].astype(str)
aLEU.insert(3, 'angle', 4)
aLEU['angle'] = aLEU['ai'] + '+' + aLEU['aj'] + '+' + aLEU['ak']
aLEU_dict = aLEU.set_index('angle')['def'].to_dict()

# dLEU = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iLEU = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'CD1', 'CA'],
                     'ak': ['C', 'CD2', 'N'],
                     'al': ['CB', 'CG', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

iLEU['ai'] = 'LEU_' + iLEU['ai'].astype(str)
iLEU['aj'] = 'LEU_' + iLEU['aj'].astype(str)
iLEU['ak'] = 'LEU_' + iLEU['ak'].astype(str)
iLEU['al'] = 'LEU_' + iLEU['al'].astype(str)
iLEU.insert(4, 'improper', 5)
iLEU['improper'] = iLEU['ai'] + '+' + iLEU['aj'] + '+' + iLEU['ak'] + '+' + iLEU['al']
iLEU_dict = iLEU.set_index('improper')['def'].to_dict()

# Serine definition.
bSER = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG', 'O', 'N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_5', 'gb_10']}
                    )

bSER['ai'] = 'SER_' + bSER['ai'].astype(str)
bSER['aj'] = 'SER_' + bSER['aj'].astype(str)
bSER.insert(2, 'bond', 3)
bSER['bond'] = bSER['ai'] + '+' + bSER['aj']
bSER_dict = bSER.set_index('bond')['def'].to_dict()

aSER = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'CB', 'OG', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

# print(aSER)
aSER['ai'] = 'SER_' + aSER['ai'].astype(str)
aSER['aj'] = 'SER_' + aSER['aj'].astype(str)
aSER['ak'] = 'SER_' + aSER['ak'].astype(str)
aSER.insert(3, 'angle', 4)
aSER['angle'] = aSER['ai'] + '+' + aSER['aj'] + '+' + aSER['ak']
# print(aSER)
aSER_dict = aSER.set_index('angle')['def'].to_dict()
# print(aSER_dict)

# dSER = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSER = pd.DataFrame({'ai': ['CA', 'C'],
                     'aj': ['N', 'CA'],
                     'ak': ['C', 'N'],
                     'al': ['CB', 'O'],
                     'def': ['gi_2', 'gi_1']}
                    )

# print(iSER)
iSER['ai'] = 'SER_' + iSER['ai'].astype(str)
iSER['aj'] = 'SER_' + iSER['aj'].astype(str)
iSER['ak'] = 'SER_' + iSER['ak'].astype(str)
iSER['al'] = 'SER_' + iSER['al'].astype(str)
iSER.insert(4, 'improper', 5)
iSER['improper'] = iSER['ai'] + '+' + iSER['aj'] + '+' + iSER['ak'] + '+' + iSER['al']
# print(iSER)
iSER_dict = iSER.set_index('improper')['def'].to_dict()
# print(iSER_dict)

# Serine terminal definition
bSERT = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'C', 'C', 'C'],
                      'aj': ['CA', 'CB', 'C', 'OG', 'O', 'N', 'OXT'],
                      'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_5', 'gb_10', 'gb_5']}
                     )

# print(bSERT)
bSERT['ai'] = 'SERT_' + bSERT['ai'].astype(str)
bSERT['aj'] = 'SERT_' + bSERT['aj'].astype(str)
bSERT.insert(2, 'bond', 3)
bSERT['bond'] = bSERT['ai'] + '+' + bSERT['aj']
# print(bSERT)
bSERT_dict = bSERT.set_index('bond')['def'].to_dict()
# print(bSERT_dict)

aSERT = pd.DataFrame({'ai': ['C', 'N', 'N', 'C', 'CA', 'CA', 'CA', 'O', 'CA', 'O'],
                      'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'C', 'C', 'C', 'C', 'C'],
                      'ak': ['CA', 'CB', 'C', 'CB', 'OG', 'O', 'N', 'N', 'OXT', 'OXT'],
                      'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_22', 'ga_22', 'ga_38']}
                     )

# print(aSERT)
aSERT['ai'] = 'SERT_' + aSERT['ai'].astype(str)
aSERT['aj'] = 'SERT_' + aSERT['aj'].astype(str)
aSERT['ak'] = 'SERT_' + aSERT['ak'].astype(str)
aSERT.insert(3, 'angle', 4)
aSERT['angle'] = aSERT['ai'] + '+' + aSERT['aj'] + '+' + aSERT['ak']
# print(aSERT)
aSERT_dict = aSERT.set_index('angle')['def'].to_dict()
# print(aSERT_dict)

# dSERT = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSERT = pd.DataFrame({'ai': ['CA', 'C'],
                      'aj': ['N', 'CA'],
                      'ak': ['C', 'N'],
                      'al': ['CB', 'O'],
                      'def': ['gi_2', 'gi_1']}
                     )

# print(iSERT)
iSERT['ai'] = 'SERT_' + iSERT['ai'].astype(str)
iSERT['aj'] = 'SERT_' + iSERT['aj'].astype(str)
iSERT['ak'] = 'SERT_' + iSERT['ak'].astype(str)
iSERT['al'] = 'SERT_' + iSERT['al'].astype(str)
iSERT.insert(4, 'improper', 5)
iSERT['improper'] = iSERT['ai'] + '+' + iSERT['aj'] + '+' + iSERT['ak'] + '+' + iSERT['al']
# print(iSERT)
iSERT_dict = iSERT.set_index('improper')['def'].to_dict()
# print(iSERT_dict)

# Proline definition.
bPRO = pd.DataFrame({'ai': ['N', 'CD', 'CA', 'CA', 'CB', 'CG', 'C', 'C'],
                     'aj': ['CA', 'N', 'CB', 'C', 'CG', 'CD', 'O', 'N'],
                     'def': ['gb_21', 'gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

bPRO['ai'] = 'PRO_' + bPRO['ai'].astype(str)
bPRO['aj'] = 'PRO_' + bPRO['aj'].astype(str)
bPRO.insert(2, 'bond', 3)
bPRO['bond'] = bPRO['ai'] + '+' + bPRO['aj']
bPRO_dict = bPRO.set_index('bond')['def'].to_dict()

aPRO = pd.DataFrame({'ai': ['C', 'C', 'CA', 'N', 'N', 'C', 'CA', 'CB', 'N', 'CA', 'CA', 'O'],
                     'aj': ['N', 'N', 'N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CD', 'C', 'C', 'C'],
                     'ak': ['CA', 'CD', 'CD', 'CB', 'C', 'CB', 'CG', 'CD', 'CG', 'O', 'N', 'N'],
                     'def': ['ga_31', 'ga_31', 'ga_21', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30',
                             'ga_19', 'ga_33']}
                    )

# print(aPRO)
aPRO['ai'] = 'PRO_' + aPRO['ai'].astype(str)
aPRO['aj'] = 'PRO_' + aPRO['aj'].astype(str)
aPRO['ak'] = 'PRO_' + aPRO['ak'].astype(str)
aPRO.insert(3, 'angle', 4)
aPRO['angle'] = aPRO['ai'] + '+' + aPRO['aj'] + '+' + aPRO['ak']
# print(aPRO)
aPRO_dict = aPRO.set_index('angle')['def'].to_dict()
# print(aPRO_dict)

# dPRO = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iPRO = pd.DataFrame({'ai': ['N', 'CA', 'C'],
                     'aj': ['C', 'N', 'CA'],
                     'ak': ['CA', 'C', 'N'],
                     'al': ['CD', 'CB', 'O'],
                     'def': ['gi_1', 'gi_2', 'gi_1']}
                    )

iPRO['ai'] = 'PRO_' + iPRO['ai'].astype(str)
iPRO['aj'] = 'PRO_' + iPRO['aj'].astype(str)
iPRO['ak'] = 'PRO_' + iPRO['ak'].astype(str)
iPRO['al'] = 'PRO_' + iPRO['al'].astype(str)
iPRO.insert(4, 'improper', 5)
iPRO['improper'] = iPRO['ai'] + '+' + iPRO['aj'] + '+' + iPRO['ak'] + '+' + iPRO['al']
iPRO_dict = iPRO.set_index('improper')['def'].to_dict()
#print(iPRO_dict)


# Creation of a dictionaries for bonds, angles and impropers between two aminoacids
b_aa = pd.DataFrame({'ai': ['TYR_C', 'THR_C', 'ILE_C', 'ALA_C', 'LEU_C', 'SER_C', 'PRO_C', 'TYR_C'],
                     'aj': ['THR_N', 'ILE_N', 'ALA_N', 'LEU_N', 'SER_N', 'PRO_N', 'TYR_N', 'SERT_N'],
                     'def': ['gb_5', 'gb_5', 'gb_5', 'gb_5', 'gb_5', 'gb_5', 'gb_5', 'gb_5']}
                    )

b_aa['bond'] = b_aa['ai'] + '+' + b_aa['aj']
b_aa_dict = b_aa.set_index('bond')['def'].to_dict()

a_aa = pd.DataFrame({'ai': ['TYR_CA', 'TYR_C', 'TYR_O', 'THR_CA', 'THR_C', 'THR_O', 'ILE_CA', 'ILE_C', 'ILE_O',
                            'ALA_CA', 'ALA_C', 'ALA_O', 'LEU_CA', 'LEU_C', 'LEU_O', 'SER_CA', 'SER_C', 'SER_C', 'SER_O',
                            'PRO_CA', 'PRO_C', 'PRO_O', 'TYR_CA', 'TYR_C', 'TYR_O'],
                     'aj': ['TYR_C', 'THR_N', 'TYR_C', 'THR_C', 'ILE_N', 'THR_C', 'ILE_C', 'ALA_N', 'ILE_C', 'ALA_C',
                            'LEU_N', 'ALA_C', 'LEU_C', 'SER_N', 'LEU_C', 'SER_C', 'PRO_N', 'PRO_N', 'SER_C', 'PRO_C',
                            'TYR_N', 'PRO_C', 'TYR_C', 'SERT_N', 'TYR_C'],
                     'ak': ['THR_N', 'THR_CA', 'THR_N', 'ILE_N', 'ILE_CA', 'ILE_N', 'ALA_N', 'ALA_CA', 'ALA_N', 'LEU_N',
                            'LEU_CA', 'LEU_N', 'SER_N', 'SER_CA', 'SER_N', 'PRO_N', 'PRO_CA', 'PRO_CD', 'PRO_N', 'TYR_N',
                            'TYR_CA', 'TYR_N', 'SERT_N', 'SERT_CA', 'SERT_N'],
                     'def': ['ga_31', 'ga_31', 'ga_33', 'ga_31', 'ga_31', 'ga_33', 'ga_31', 'ga_31', 'ga_33', 'ga_31',
                             'ga_31', 'ga_33', 'ga_31', 'ga_31', 'ga_33', 'ga_31', 'ga_31', 'ga_31', 'ga_33', 'ga_31',
                             'ga_31', 'ga_33', 'ga_31', 'ga_31', 'ga_33']}
                    )

a_aa['angle'] = a_aa['ai'] + '+' + a_aa['aj'] + '+' + a_aa['ak']
a_aa_dict = a_aa.set_index('angle')['def'].to_dict()



# Merging all the dictionary into one to replace in functions.py
aa_bonds = {**bALA_dict, **bTYR_dict, **bILE_dict, **bLEU_dict, **bPRO_dict, **bSER_dict, **bSERT_dict, **bTHR_dict,
            **b_aa_dict}
aa_angles = {**aALA_dict, **aTYR_dict, **aILE_dict, **aLEU_dict, **aPRO_dict, **aSER_dict, **aSERT_dict, **aTHR_dict,
             **a_aa_dict}
aa_impropers = {**iALA_dict, **iTYR_dict, **iILE_dict, **iLEU_dict, **iPRO_dict, **iSER_dict, **iSERT_dict, **iTHR_dict}
#print(aa_impropers)

# bALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'def': ['']}
#                    )

# aALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'def': ['']}
#                    )

# aALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )
