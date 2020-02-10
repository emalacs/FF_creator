#! /usr/bin/python

import numpy as np
import pandas as pd
import prettytable as pt
import copy
import re

pairs = pd.read_csv('ffnonbonded', sep="\s+", header=None)
#sistemare la notazine scientifica
pairs.columns = ["ai", "aj", "type", "A", "B"]
ai_in = pairs["ai"].str.split("_", n = 1, expand = True)
aj_jn = pairs["aj"].str.split("_", n = 1, expand = True)
pairs["ain"] = ai_in[1]
pairs["air"] = ai_in[0]
pairs["ajn"] = aj_jn[1]
pairs["ajr"] = aj_jn[0]
#pairs["merge"] = pairs["air"].apply(str) + '_' + pairs["ajr"].apply(str)
pairs.astype(str)

print("Selection of backbone, sidechains and mixed")

back_atomtype = ["CA", "N", "O", "C", "OXT"]

# backbone
drop_sidechains = pairs[pairs['air'].isin(back_atomtype)]
backbone = drop_sidechains[drop_sidechains['ajr'].isin(back_atomtype)]

# sidechains
drop_backbone = pairs[~ pairs['air'].isin(back_atomtype)]
sidechains = drop_backbone[~ drop_backbone['ajr'].isin(back_atomtype)]

# mixed #con questo metodo non funziona
mix_air = pairs[pairs['air'].isin(back_atomtype)]
mix_ajr = pairs[pairs['ajr'].isin(back_atomtype)]
mix_full = mix_air.append(mix_ajr, sort=False, ignore_index=True)

# mixed poi ci penso


pymol_backbone = backbone.copy()
pymol_sidechains = sidechains.copy()
pymol_mixed = mix_full.copy()

# cleaning the dataframes
backbone = backbone.drop(['air', 'ajr', 'ain', 'ajn'], axis=1)
sidechains = sidechains.drop(['air', 'ajr', 'ain', 'ajn'], axis=1)
mix_full = mix_full.drop(['air', 'ajr', 'ain', 'ajn'], axis=1)


print("Writing backbone, sidechains and mixed files")

file = open("backbone", "w")
file.write(str(backbone.to_string(index=False, header=False)))
file.close()

file = open("sidechains", "w")
file.write(str(sidechains.to_string(index=False, header=False)))
file.close()

file = open("mixed", "w")
file.write(str(mix_full.to_string(index=False, header=False)))
file.close()

print ("Backbone, sidechains and mixed wrote")
print ("Creating three pymol file to paste")
# distance prova, resi 10 and name OH, resi 3 and name CG1, 6 -> comando python con cutoff 8

# backbone
pymol_backbone = pymol_backbone.drop(['type', 'A', 'B'], axis=1)
pymol_backbone['interaction'] = np.where((pymol_backbone['ai'] == pymol_backbone['aj']), 'v', 'h')
name =  pymol_backbone['interaction'].apply(str) + '_' + pymol_backbone['ai'].apply(str) + ':' + pymol_backbone['aj'].apply(str)
pymol_backbone = pymol_backbone.drop(['ai', 'aj'], axis=1)
pymol_backbone.insert(0, 'distance', 'distance')
pymol_backbone.insert(1, 'name', name)
pymol_backbone.insert(2, 'resi', ', resi')
pymol_backbone.insert(4, 'and', 'and name')
pymol_backbone.insert(6, 'resi2', ', resi')
pymol_backbone.insert(8, 'and2', 'and name')
pymol_backbone.insert(10, 'cutoff', ', 6')
pymol_backbone = pymol_backbone.sort_values(by='interaction')
pymol_backbone = pymol_backbone.drop(['interaction'], axis=1)

# sidechains
pymol_sidechains = pymol_sidechains.drop(['type', 'A', 'B'], axis=1)
pymol_sidechains['interaction'] = np.where((pymol_sidechains['ai'] == pymol_sidechains['aj']), 'v', 'h')
name =  pymol_sidechains['interaction'].apply(str) + '_' + pymol_sidechains['ai'].apply(str) + ':' + pymol_sidechains['aj'].apply(str)
pymol_sidechains = pymol_sidechains.drop(['ai', 'aj'], axis=1)
pymol_sidechains.insert(0, 'distance', 'distance')
pymol_sidechains.insert(1, 'name', name)
pymol_sidechains.insert(2, 'resi', ', resi')
pymol_sidechains.insert(4, 'and', 'and name')
pymol_sidechains.insert(6, 'resi2', ', resi')
pymol_sidechains.insert(8, 'and2', 'and name')
pymol_sidechains.insert(10, 'cutoff', ', 6')
pymol_sidechains = pymol_sidechains.sort_values(by='interaction')
pymol_sidechains = pymol_sidechains.drop(['interaction'], axis=1)

# mixed
pymol_mixed = pymol_mixed.drop(['type', 'A', 'B'], axis=1)
pymol_mixed['interaction'] = np.where((pymol_mixed['ai'] == pymol_mixed['aj']), 'v', 'h')
name =  pymol_mixed['interaction'].apply(str) + '_' + pymol_mixed['ai'].apply(str) + ':' + pymol_mixed['aj'].apply(str)
pymol_mixed = pymol_mixed.drop(['ai', 'aj'], axis=1)
pymol_mixed.insert(0, 'distance', 'distance')
pymol_mixed.insert(1, 'name', name)
pymol_mixed.insert(2, 'resi', ', resi')
pymol_mixed.insert(4, 'and', 'and name')
pymol_mixed.insert(6, 'resi2', ', resi')
pymol_mixed.insert(8, 'and2', 'and name')
pymol_mixed.insert(10, 'cutoff', ', 6')
pymol_mixed = pymol_mixed.sort_values(by='interaction')
pymol_mixed = pymol_mixed.drop(['interaction'], axis=1)

print("Writing the two pymol files")

file = open("pymol_backbone", "w")
file.write(str(pymol_backbone.to_string(index=False, header=False)))
file.close()


file = open("pymol_sidechains", "w")
file.write(str(pymol_sidechains.to_string(index=False, header=False)))
file.close()

file = open("pymol_mixed", "w")
file.write(str(pymol_mixed.to_string(index=False, header=False)))
file.close()

