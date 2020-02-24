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
            2.6646244e-05, 2.8058209e-05, 1.2e-05],
     'bonds': [],
     'angles': [],
     'improper': []
     # da aggiungere le definizioni di bonded, angles e improper
     }
)

gromos_atp.to_dict()
gromos_atp.set_index('name', inplace = True)

# TTR sequence YTIAALLSPYS
# Alanine definition.
ALA = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH3', 'C', 'O']}
                   )
ALA.to_dict()

# Threonine definition
THR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG1', 'HG1', 'CG2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'OA', 'H', 'CH3', 'C', 'O']}
                   )
THR.to_dict()

# Tyrosine definition.
TYR = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ',
                            'OH', 'HH', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'C', 'C', 'HC', 'C', 'HC', 'C', 'HC', 'C',
                            'HC', 'C', 'OA', 'H', 'C', 'O']}
                   )

TYR.to_dict()

# Isoleucine definition.
ILE = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH1', 'CH2', 'CH3', 'CH3', 'C', 'O']}
                   )
ILE.to_dict()

# Leucine definition.
LEU = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'CH1', 'CH3', 'CH3', 'C', 'O']}
                   )
LEU.to_dict()

# Serine definition.
SER = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O'],
                    'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O']}
                   )
SER.to_dict()

# Serine terminal definition
SERT = pd.DataFrame({'type': ['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O', 'OXT'],
                     'name': ['N', 'H', 'CH1', 'CH2', 'OA', 'H', 'C', 'O', 'O']}
                    )
SERT.to_dict()

# Proline definition.
PRO = pd.DataFrame({'type': ['N', 'CA', 'CB', 'CG', 'CD', 'C', 'O'],
                    'name': ['N', 'CH1', 'CH2r', 'CH2r', 'CH2r', 'C', 'O']}
                   )
PRO.to_dict()

gromos_aa = {}
reslist = [ALA, TYR, ILE, LEU, PRO, SER, SERT, THR]
reslistname = ['ALA', 'TYR', 'ILE', 'LEU', 'PRO', 'SER', 'SERT', 'THR']
for res, resn in zip(reslist, reslistname):
    for i in range(len(res['type'])):
        str = resn + '_' + res.iloc[i, 0]
        gromos_aa[str] = res.iloc[i, 1]
