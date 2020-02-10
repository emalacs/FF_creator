import pandas as pd

#funzioni con qualcosa di chiaro in ingresso e qualcosa di chiaro in uscita


#for f in os.listdir(args.input_folder):
    #if
#    if f.endswith('.top'):
        #Call my function


df = pd.read_csv('fibril8.top', header=None, names=range(11))
table_names = ['[ defaults ]',
               '[ atomtypes ]',
                '[ moleculetype ]',
                '[ atoms ]',
                '[ bonds ]',
                '[ angles ]',
                '[ dihedrals ]',
                '[ pairs ]',
                '[ exclusions ]',
                '[ system ]',
                '[ molecules ]']
groups = df[0].isin(table_names).cumsum()
tables = {g.iloc[0,0]: g.iloc[1:] for k,g in df.groupby(groups)}

atoms = tables['[ atoms ]']
print(atoms)
atoms["charge"] = ""
atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]
print(atoms)