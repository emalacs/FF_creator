import pandas as pd

pep_atoms = pd.read_csv('input/pep_atoms', sep = '\s+', header = None)

dat = pd.read_csv('analysis/histo-time.dat', sep = '\s+')

print(dat.to_string())
