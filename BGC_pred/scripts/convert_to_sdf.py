from rdkit import Chem
from IPython.display import SVG
from rdkit.Chem import AllChem, Draw, Descriptors, PandasTools
import os
import sys

smiles_data = sys.argv[1]
lpath = sys.argv[2]

#smiles_data = 'CC1=CC(=NN1CC(=O)N2CCC(CC2)C3=NC(=CS3)C4=NOC(C4)C5=C(C=CC=C5F)F)C(F)(F)F'
#sdf_path = '/home/seika_oiwa_3590/notebooks/Mydata_analysis/Antismash/Docking/ligand/tmp.sdf'

sdf_path = f'{lpath}/tmp.sdf'

mh = Chem.AddHs(Chem.MolFromSmiles(str(smiles_data)))
AllChem.EmbedMolecule(mh, AllChem.ETKDGv2()) 
writer = Chem.SDWriter(sdf_path)
writer.write(mh)