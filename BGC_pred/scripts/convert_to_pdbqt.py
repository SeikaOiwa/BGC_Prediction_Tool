import oddt
from oddt import toolkit
from oddt import docking
import os
import sys

lpath = sys.argv[1]

sdf_path = f'{lpath}/tmp.sdf'

ligands = next(oddt.toolkit.readfile('sdf',sdf_path))
oddt.docking.AutodockVina.write_vina_pdbqt(ligands,lpath,name_id='ligand')
os.remove(sdf_path)