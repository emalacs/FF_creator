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

# Alananine angles definition.
aALA = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_30', 'ga_19', 'ga_33']}
                    )

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

# Threonine definition.
THR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG1', 'HG1', 'CG2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'OA', 'H', 'CH3', 'C', 'O']}
                   )

bTHR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'OG1', 'CG2', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_18', 'gb_27', 'gb_5', 'gb_10']}
                    )

aTHR = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'OG1', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'OG1', 'CG2', 'CG2', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_30', 'ga_19', 'ga_33']}
                    )

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

THR.to_dict()

# Tyrosine definition.
TYR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ',
                             'OH', 'HH', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'C', 'C', 'HC', 'C', 'HC', 'C', 'HC', 'C',
                             'HC', 'C', 'OA', 'H', 'C', 'O']}
                   )

bTYR = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'OH', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16', 'gb_16',
                             'gb_13', 'gb_5', 'gb_10']}
                    )

aTYR = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CB', 'CB', 'CD1', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE1',
                            'CE2', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CG', 'CG', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ',
                            'CZ', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'CG', 'CD1', 'CD2', 'CD2', 'CE1', 'CE2', 'CZ', 'CZ', 'CE2', 'OH',
                            'OH', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27',
                             'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_27', 'ga_19', 'ga_33']}
                    )

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

TYR.to_dict()

# Isoleucine definition.
ILE = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'CH2', 'CH3', 'CH3', 'C', 'O']}
                   )

bILE = pd.DataFrame({'ai': ['N', 'CA', 'CA', 'CB', 'CB', 'CG1', 'C', 'C'],
                     'aj': ['CA', 'CB', 'C', 'CG1', 'CG2', 'CD', 'O', '+N'],
                     'def': ['gb_21', 'gb_27', 'gb_27', 'gb_27', 'gb_27', 'gb_27',  'gb_5',  'gb_10']}
                    )

aILE = pd.DataFrame({'ai': ['-C', 'N', 'N', 'CB', 'CA', 'CA', 'CG1', 'CB', 'CA', 'CA', 'O'],
                     'aj': ['N', 'CA', 'CA', 'CA', 'CB', 'CB', 'CB', 'CG1', 'C', 'C', 'C'],
                     'ak': ['CA', 'CB', 'C', 'C', 'CG1', 'CG2', 'CG2', 'CD', 'O', '+N', '+N'],
                     'def': ['ga_31', 'ga_13', 'ga_13', 'ga_13', 'ga_15', 'ga_15', 'ga_15', 'ga_15', 'ga_30',
                             'ga_19', 'ga_33']}
                    )

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

ILE.to_dict()

# Leucine definition.
LEU = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'CH1', 'CH3', 'CH3', 'C', 'O']}
                   )

bLEU = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'def': ['']}
                    )

aLEU = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'def': ['']}
                    )

#dLEU = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iLEU = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'al': [''],
                     'def': ['']}
                    )

LEU.to_dict()

# Serine definition.
SER = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O']}
                   )

bSER = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'def': ['']}
                    )

aSER = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'def': ['']}
                    )

#dSER = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSER = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'al': [''],
                     'def': ['']}
                    )

SER.to_dict()

# Serine terminal definition
SERT = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O', 'OXT'],
                     'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O', 'O']}
                    )

bSERT = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'def': ['']}
                    )

aSERT = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'def': ['']}
                    )

#dSERT = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iSERT = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'al': [''],
                     'def': ['']}
                    )

SERT.to_dict()

# Proline definition.
PRO = pd.DataFrame({'type': ['N', 'CA', 'CB', 'CG', 'CD', 'C', 'O'],
                    'name': ['N', 'CH1', 'CH2r', 'CH2r', 'CH2r', 'C', 'O']}
                   )

bPRO = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'def': ['']}
                    )

aPRO = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'def': ['']}
                    )

#dPRO = pd.DataFrame({'ai': [''],
#                     'aj': [''],
#                     'ak': [''],
#                     'al': [''],
#                     'def': ['']}
#                    )

iPRO = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'al': [''],
                     'def': ['']}
                    )

PRO.to_dict()

gromos_aa = {}
reslist = [ALA, TYR, ILE, LEU, PRO, SER, SERT, THR]
reslistname = ['ALA', 'TYR', 'ILE', 'LEU', 'PRO', 'SER', 'SERT', 'THR']
for res, resn in zip(reslist, reslistname):
    for i in range(len(res['type'])):
        str = resn + '_' + res.iloc[i, 0]
        gromos_aa[str] = res.iloc[i, 1]

bALA = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'def': ['']}
                    )

aALA = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'def': ['']}
                    )

aALA = pd.DataFrame({'ai': [''],
                     'aj': [''],
                     'ak': [''],
                     'al': [''],
                     'def': ['']}
                    )