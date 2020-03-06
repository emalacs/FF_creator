import pandas as pd

# Script containing the definitions of the atomtypes of every aminoacid used in here. Taken from gromos/aminoacids.rtp
# along with the mass and c12 for every atomtype contained in Gromos

# Gromos atomtypes definition
# Abbiamo tenuto solo le info che ci servivano. Da aggiungere le definizioni dei terminali e quelli che serviranno
# per i prossimi aminoacidi. Abbiamo aggiunto gli H nella massa (e nel c12) - OA ed N -
# per aggiungere C12 "regola del 20".

# GROMOS 53a6

# THE C12 RATIO CHANGED a little bit.
gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 'CH2', 'CH3', 'CH2r', 'NT'],
     'mass': [16, 17, 15, 12, 13, 14, 15, 14, 17],
     'at.num': [8, 8, 7, 6, 6, 6, 6, 6, 7],
     'c12': [1e-06, 3.011e-05, 4.639e-05, 4.937284e-06, 9.70225e-05, 3.3965584e-05,
             2.6646244e-05, 2.8058209e-05, 1.2e-05]
     }
)

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
    {'def': ['ga_13', 'ga_15', 'ga_19', 'ga_27', 'ga_30', 'ga_31', 'ga_33'],
     'force': ['520.00', '530.00', '610.00', '560.00', '685.00', '700.00', '730.00'],
     'len': ['109.50', '111.00', '115.00', '120.00', '121.00', '122.00', '124.00']
    }
)

gromos_dihedrals = pd.DataFrame(
    {'def': [''],
     'force': [''],
     'len': ['']
    }
)

gromos_impropers = pd.DataFrame(
    {'def': ['gi_1', 'gi_2'],
     'force': ['167.42309', '334.84617'],
     'len': ['0.0', '35.26439']
    }
)

gromos_atp.to_dict()
gromos_atp.set_index('name', inplace=True)

# TTR sequence YTIAALLSPYS

# Alanine atomtypes definition.
ALA = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH3', 'C', 'O']}
                   )
ALA.to_dict()

# Alanine bonds definition.
bALA = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'O', 'N+'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )
print(bALA)
bALA['ai'] = 'ALA_' + bALA['ai'].astype(str)
bALA['aj'] = 'ALA_' + bALA['aj'].astype(str)
bALA.insert(2, 'bond', 3)
bALA['bond'] = bALA['ai'] + '+' + bALA['aj']
print(bALA)
bALA_dict = bALA.set_index('bond')['def'].to_dict()
print(bALA_dict)

# Alananine angles definition.
aALA = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

print(aALA)
aALA['ai'] = 'ALA_' + aALA['ai'].astype(str)
aALA['aj'] = 'ALA_' + aALA['aj'].astype(str)
aALA['ak'] = 'ALA_' + aALA['ak'].astype(str)
aALA.insert(3, 'angle', 4)
aALA['angle'] = aALA['ai'] + '+' + aALA['aj'] + '+' + aALA['ak']
print(aALA)
aALA_dict = aALA.set_index('angle')['def'].to_dict()
print(aALA_dict)

# Alanine proper dihedrals definition.
dALA = pd.DataFrame({'ai': ['-CA', '-C', 'N'],
                     'aj': ['-C', 'N', 'CA'],
                     'ak': ['N', 'CA', 'C'],
                     'al': ['CA', 'C', '+N'],
                     'def': ['gd_14', 'gd_39', 'gd_40']}
                    )

# Alananine improper dihedrals definition.
iALA = pd.DataFrame({'ai': ['CA', 'C'],
                     'aj': ['N', 'CA'],
                     'ak': ['C', '+N'],
                     'al': ['CB', 'O'],
                     'def': ['gi_2', 'gi_1']}
                    )

print(iALA)
iALA['ai'] = 'ALA_' + iALA['ai'].astype(str)
iALA['aj'] = 'ALA_' + iALA['aj'].astype(str)
iALA['ak'] = 'ALA_' + iALA['ak'].astype(str)
iALA['al'] = 'ALA_' + iALA['al'].astype(str)
iALA.insert(4, 'improper', 5)
iALA['improper'] = iALA['ai'] + '+' + iALA['aj'] + '+' + iALA['ak'] + '+' + iALA['al']
print(iALA)
iALA_dict = iALA.set_index('improper')['def'].to_dict()
print(iALA_dict)

# Threonine definition.
THR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG1', 'HG1', 'CG2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'OA', 'H', 'CH3', 'C', 'O']}
                   )

THR.to_dict()

bTHR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG1', 'CG2', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_27', 'gb_5', 'gb_10']}
                    )

print(bTHR)
bTHR['ai'] = 'THR_' + bTHR['ai'].astype(str)
bTHR['aj'] = 'THR_' + bTHR['aj'].astype(str)
bTHR.insert(2, 'bond', 3)
bTHR['bond'] = bTHR['ai'] + '+' + bTHR['aj']
print(bTHR)
bTHR_dict = bTHR.set_index('bond')['def'].to_dict()
print(bTHR_dict)

aTHR = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'OG1', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'OG1', 'CG2', 'CG2', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_30', 'ga_19', 'ga_33']}
                    )

print(aTHR)
aTHR['ai'] = 'THR_' + aTHR['ai'].astype(str)
aTHR['aj'] = 'THR_' + aTHR['aj'].astype(str)
aTHR['ak'] = 'THR_' + aTHR['ak'].astype(str)
aTHR.insert(3, 'angle', 4)
aTHR['angle'] = aTHR['ai'] + '+' + aTHR['aj'] + '+' + aTHR['ak']
print(aTHR)
aTHR_dict = aTHR.set_index('angle')['def'].to_dict()
print(aTHR_dict)

#dTHR = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iTHR = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'OG1', 'CA'],
                     'ak': ['C', 'CG2', '+N'],
                     'al': ['CB', 'CA', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

print(iTHR)
iTHR['ai'] = 'THR_' + iTHR['ai'].astype(str)
iTHR['aj'] = 'THR_' + iTHR['aj'].astype(str)
iTHR['ak'] = 'THR_' + iTHR['ak'].astype(str)
iTHR['al'] = 'THR_' + iTHR['al'].astype(str)
iTHR.insert(4, 'improper', 5)
iTHR['improper'] = iTHR['ai'] + '+' + iTHR['aj'] + '+' + iTHR['ak'] + '+' + iTHR['al']
print(iTHR)
iTHR_dict = iTHR.set_index('improper')['def'].to_dict()
print(iTHR_dict)


# Tyrosine definition.
TYR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ',
                             'OH', 'HH', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'C', 'C', 'HC', 'C', 'HC', 'C', 'HC', 'C',
                             'HC', 'C', 'OA', 'H', 'C', 'O']}
                   )

TYR.to_dict()

bTYR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'OH', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16',
                             'gb_13', 'gb_5', 'gb_10']}
                    )

print(bTYR)
bTYR['ai'] = 'TYR_' + bTYR['ai'].astype(str)
bTYR['aj'] = 'TYR_' + bTYR['aj'].astype(str)
bTYR.insert(2, 'bond', 3)
bTYR['bond'] = bTYR['ai'] + '+' + bTYR['aj']
print(bTYR)
bTYR_dict = bTYR.set_index('bond')['def'].to_dict()
print(bTYR_dict)

aTYR = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CB', 'CB', 'CD1', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE1',
                            'CE2', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ',
                            'CZ', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'CG', 'CD1', 'CD2', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'CE2', 'OH',
                            'OH', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27',
                             'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_30', 'ga_19', 'ga_33']}
                    )

print(aTYR)
aTYR['ai'] = 'TYR_' + aTYR['ai'].astype(str)
aTYR['aj'] = 'TYR_' + aTYR['aj'].astype(str)
aTYR['ak'] = 'TYR_' + aTYR['ak'].astype(str)
aTYR.insert(3, 'angle', 4)
aTYR['angle'] = aTYR['ai'] + '+' + aTYR['aj'] + '+' + aTYR['ak']
print(aTYR)
aTYR_dict = aTYR.set_index('angle')['def'].to_dict()
print(aTYR_dict)

#dTYR = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iTYR = pd.DataFrame({'ai': ['CA', 'CG', 'CG', 'CG', 'CD1', 'CD1', 'CD2', 'CD2', 'CZ', 'C'],
                     'aj': ['N', 'CD1', 'CD1', 'CD2', 'CG', 'CE1', 'CG', 'CE2', 'CE1', 'CA'],
                     'ak': ['C', 'CD2', 'CE1', 'CE2', 'CD2', 'CZ', 'CD1', 'CZ', 'CE2', '+N'],
                     'al': ['CB', 'CB', 'CZ', 'CZ', 'CE2', 'CE2', 'CE1', 'CE1', 'OH', 'O'],
                     'def': ['gi_2', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1', 'gi_1']}
                    )

print(iTYR)
iTYR['ai'] = 'TYR_' + iTYR['ai'].astype(str)
iTYR['aj'] = 'TYR_' + iTYR['aj'].astype(str)
iTYR['ak'] = 'TYR_' + iTYR['ak'].astype(str)
iTYR['al'] = 'TYR_' + iTYR['al'].astype(str)
iTYR.insert(4, 'improper', 5)
iTYR['improper'] = iTYR['ai'] + '+' + iTYR['aj'] + '+' + iTYR['ak'] + '+' + iTYR['al']
print(iTYR)
iTYR_dict = iTYR.set_index('improper')['def'].to_dict()
print(iTYR_dict)


# Isoleucine definition.
ILE = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'CH2', 'CH3', 'CH3', 'C', 'O']}
                   )

ILE.to_dict()

bILE = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'CG1', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG1', 'CG2', 'CD', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_27',  'gb_5',  'gb_10']}
                    )

print(bILE)
bILE['ai'] = 'ILE_' + bILE['ai'].astype(str)
bILE['aj'] = 'ILE_' + bILE['aj'].astype(str)
bILE.insert(2, 'bond', 3)
bILE['bond'] = bILE['ai'] + '+' + bILE['aj']
print(bILE)
bILE_dict = bILE.set_index('bond')['def'].to_dict()
print(bILE_dict)

aILE = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'CG1', 'CB', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'CG1', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'CG1', 'CG2', 'CG2', 'CD', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_15', 'ga_15', 'ga_30',
                             'ga_19', 'ga_33']}
                    )

print(aILE)
aILE['ai'] = 'ILE_' + aILE['ai'].astype(str)
aILE['aj'] = 'ILE_' + aILE['aj'].astype(str)
aILE['ak'] = 'ILE_' + aILE['ak'].astype(str)
aILE.insert(3, 'angle', 4)
aILE['angle'] = aILE['ai'] + '+' + aILE['aj'] + '+' + aILE['ak']
print(aILE)
aILE_dict = aILE.set_index('angle')['def'].to_dict()
print(aILE_dict)

#dILE = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iILE = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'CG1', 'CA'],
                     'ak': ['C', 'CG2', '+N'],
                     'al': ['CB', 'CA', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

print(iILE)
iILE['ai'] = 'ILE_' + iILE['ai'].astype(str)
iILE['aj'] = 'ILE_' + iILE['aj'].astype(str)
iILE['ak'] = 'ILE_' + iILE['ak'].astype(str)
iILE['al'] = 'ILE_' + iILE['al'].astype(str)
iILE.insert(4, 'improper', 5)
iILE['improper'] = iILE['ai'] + '+' + iILE['aj'] + '+' + iILE['ak'] + '+' + iILE['al']
print(iILE)
iILE_dict = iILE.set_index('improper')['def'].to_dict()
print(iILE_dict)

# Leucine definition.
LEU = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'CH1', 'CH3', 'CH3', 'C', 'O']}
                   )

LEU.to_dict()

bLEU = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CG', 'CG', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG', 'CD1', 'CD2', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

print(bLEU)
bLEU['ai'] = 'LEU_' + bLEU['ai'].astype(str)
bLEU['aj'] = 'LEU_' + bLEU['aj'].astype(str)
bLEU.insert(2, 'bond', 3)
bLEU['bond'] = bLEU['ai'] + '+' + bLEU['aj']
print(bLEU)
bLEU_dict = bLEU.set_index('bond')['def'].to_dict()
print(bLEU_dict)

aLEU = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CB', 'CB', 'CD1', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CG', 'CG', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'CG', 'CD1', 'CD2', 'CD2', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_15', 'ga_15', 'ga_30', 'ga_19',
                             'ga_33']}
                    )

print(aLEU)
aLEU['ai'] = 'LEU_' + aLEU['ai'].astype(str)
aLEU['aj'] = 'LEU_' + aLEU['aj'].astype(str)
aLEU['ak'] = 'LEU_' + aLEU['ak'].astype(str)
aLEU.insert(3, 'angle', 4)
aLEU['angle'] = aLEU['ai'] + '+' + aLEU['aj'] + '+' + aLEU['ak']
print(aLEU)
aLEU_dict = aLEU.set_index('angle')['def'].to_dict()
print(aLEU_dict)

#dLEU = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iLEU = pd.DataFrame({'ai': ['CA', 'CB', 'C'],
                     'aj': ['N', 'CD1', 'CA'],
                     'ak': ['C', 'CD2', '+N'],
                     'al': ['CB', 'CG', 'O'],
                     'def': ['gi_2', 'gi_2', 'gi_1']}
                    )

print(iLEU)
iLEU['ai'] = 'LEU_' + iLEU['ai'].astype(str)
iLEU['aj'] = 'LEU_' + iLEU['aj'].astype(str)
iLEU['ak'] = 'LEU_' + iLEU['ak'].astype(str)
iLEU['al'] = 'LEU_' + iLEU['al'].astype(str)
iLEU.insert(4, 'improper', 5)
iLEU['improper'] = iLEU['ai'] + '+' + iLEU['aj'] + '+' + iLEU['ak'] + '+' + iLEU['al']
print(iLEU)
iLEU_dict = iLEU.set_index('improper')['def'].to_dict()
print(iLEU_dict)

# Serine definition.
SER = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O']}
                   )

SER.to_dict()

bSER = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_5', 'gb_10']}
                    )

print(bSER)
bSER['ai'] = 'SER_' + bSER['ai'].astype(str)
bSER['aj'] = 'SER_' + bSER['aj'].astype(str)
bSER.insert(2, 'bond', 3)
bSER['bond'] = bSER['ai'] + '+' + bSER['aj']
print(bSER)
bSER_dict = bSER.set_index('bond')['def'].to_dict()
print(bSER_dict)

aSER = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'OG', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

print(aSER)
aSER['ai'] = 'SER_' + aSER['ai'].astype(str)
aSER['aj'] = 'SER_' + aSER['aj'].astype(str)
aSER['ak'] = 'SER_' + aSER['ak'].astype(str)
aSER.insert(3, 'angle', 4)
aSER['angle'] = aSER['ai'] + '+' + aSER['aj'] + '+' + aSER['ak']
print(aSER)
aSER_dict = aSER.set_index('angle')['def'].to_dict()
print(aSER_dict)

#dSER = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSER = pd.DataFrame({'ai': ['CA', 'C'],
                     'aj': ['N', 'CA'],
                     'ak': ['C', '+N'],
                     'al': ['CB', 'O'],
                     'def': ['gi_2', 'gi_1']}
                    )

print(iSER)
iSER['ai'] = 'SER_' + iSER['ai'].astype(str)
iSER['aj'] = 'SER_' + iSER['aj'].astype(str)
iSER['ak'] = 'SER_' + iSER['ak'].astype(str)
iSER['al'] = 'SER_' + iSER['al'].astype(str)
iSER.insert(4, 'improper', 5)
iSER['improper'] = iSER['ai'] + '+' + iSER['aj'] + '+' + iSER['ak'] + '+' + iSER['al']
print(iSER)
iSER_dict = iSER.set_index('improper')['def'].to_dict()
print(iSER_dict)

# Serine terminal definition
SERT = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O', 'OXT'],
                     'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O', 'O']}
                    )

SERT.to_dict()

bSERT = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_5', 'gb_10']}
                    )

print(bSERT)
bSERT['ai'] = 'SERT_' + bSERT['ai'].astype(str)
bSERT['aj'] = 'SERT_' + bSERT['aj'].astype(str)
bSERT.insert(2, 'bond', 3)
bSERT['bond'] = bSERT['ai'] + '+' + bSERT['aj']
print(bSERT)
bSERT_dict = bSERT.set_index('bond')['def'].to_dict()
print(bSERT_dict)

aSERT = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'OG', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

print(aSERT)
aSERT['ai'] = 'SERT_' + aSERT['ai'].astype(str)
aSERT['aj'] = 'SERT_' + aSERT['aj'].astype(str)
aSERT['ak'] = 'SERT_' + aSERT['ak'].astype(str)
aSERT.insert(3, 'angle', 4)
aSERT['angle'] = aSERT['ai'] + '+' + aSERT['aj'] + '+' + aSERT['ak']
print(aSERT)
aSERT_dict = aSERT.set_index('angle')['def'].to_dict()
print(aSERT_dict)

#dSERT = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSERT = pd.DataFrame({'ai': ['CA', 'C'],
                     'aj': ['N', 'CA'],
                     'ak': ['C', '+N'],
                     'al': ['CB', 'O'],
                     'def': ['gi_2', 'gi_1']}
                    )

print(iSERT)
iSERT['ai'] = 'SERT_' + iSERT['ai'].astype(str)
iSERT['aj'] = 'SERT_' + iSERT['aj'].astype(str)
iSERT['ak'] = 'SERT_' + iSERT['ak'].astype(str)
iSERT['al'] = 'SERT_' + iSERT['al'].astype(str)
iSERT.insert(4, 'improper', 5)
iSERT['improper'] = iSERT['ai'] + '+' + iSERT['aj'] + '+' + iSERT['ak'] + '+' + iSERT['al']
print(iSERT)
iSERT_dict = iSERT.set_index('improper')['def'].to_dict()
print(iSERT_dict)

# Proline definition.
PRO = pd.DataFrame({'type': ['N', 'CA', 'CB', 'CG', 'CD', 'C', 'O'],
                    'name': ['N', 'CH1', 'CH2r', 'CH2r', 'CH2r', 'C', 'O']}
                   )

PRO.to_dict()

bPRO = pd.DataFrame({'ai': ['N', 'N', 'CA', 'CA', 'CB', 'CG', 'C', 'C'],
                     'aj': ['CA', 'CD', 'CB', 'C', 'CG', 'CD', 'O', '+N'],
                     'def': ['gb_21', 'gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_5', 'gb_10']}
                    )

print(bPRO)
bPRO['ai'] = 'PRO_' + bPRO['ai'].astype(str)
bPRO['aj'] = 'PRO_' + bPRO['aj'].astype(str)
bPRO.insert(2, 'bond', 3)
bPRO['bond'] = bPRO['ai'] + '+' + bPRO['aj']
print(bPRO)
bPRO_dict = bPRO.set_index('bond')['def'].to_dict()
print(bPRO_dict)

aPRO = pd.DataFrame({'ai': ['-C', '-C', 'CA', 'N', 'N', 'CB', 'CA', 'CB', 'N', 'CA', 'CA', 'O'],
                     'aj': ['N', 'N', 'N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CD', 'C', 'C', 'C'],
                     'ak': ['CA', 'CD', 'CD', 'CB', 'C', 'C', 'CG', 'CD', 'CG', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_31', 'ga_21', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_30',
                             'ga_19', 'ga_33']}
                    )

print(aPRO)
aPRO['ai'] = 'PRO_' + aPRO['ai'].astype(str)
aPRO['aj'] = 'PRO_' + aPRO['aj'].astype(str)
aPRO['ak'] = 'PRO_' + aPRO['ak'].astype(str)
aPRO.insert(3, 'angle', 4)
aPRO['angle'] = aPRO['ai'] + '+' + aPRO['aj'] + '+' + aPRO['ak']
print(aPRO)
aPRO_dict = aPRO.set_index('angle')['def'].to_dict()
print(aPRO_dict)

#dPRO = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iPRO = pd.DataFrame({'ai': ['N', 'CA', 'C'],
                     'aj': ['-C', 'N', 'CA'],
                     'ak': ['CA', 'C', '+N'],
                     'al': ['CD', 'CB', 'O'],
                     'def': ['gi_1', 'gi_2', 'gi_1']}
                    )

print(iPRO)
iPRO['ai'] = 'PRO_' + iPRO['ai'].astype(str)
iPRO['aj'] = 'PRO_' + iPRO['aj'].astype(str)
iPRO['ak'] = 'PRO_' + iPRO['ak'].astype(str)
iPRO['al'] = 'PRO_' + iPRO['al'].astype(str)
iPRO.insert(4, 'improper', 5)
iPRO['improper'] = iPRO['ai'] + '+' + iPRO['aj'] + '+' + iPRO['ak'] + '+' + iPRO['al']
print(iPRO)
iPRO_dict = iPRO.set_index('improper')['def'].to_dict()
print(iPRO_dict)


aa_bonds = {**bALA_dict, **bTYR_dict, **bILE_dict, **bLEU_dict, **bPRO_dict, **bSER_dict, **bSERT_dict, **bTHR_dict}
print(aa_bonds)



gromos_aa = {}
reslist = [ALA, TYR, ILE, LEU, PRO, SER, SERT, THR]
reslistname = ['ALA', 'TYR', 'ILE', 'LEU', 'PRO', 'SER', 'SERT', 'THR']
for res, resn in zip(reslist, reslistname):
    for i in range(len(res['type'])):
        str = resn + '_' + res.iloc[i, 0]
        gromos_aa[str] = res.iloc[i, 1]
# print(gromos_aa)

#bALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'def': ['']}
#                    )

#aALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'def': ['']}
#                    )

#aALA = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )