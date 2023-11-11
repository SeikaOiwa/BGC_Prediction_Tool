from vina import Vina
import os
import sys

denzyme_path = sys.argv[1]
dligand_path = sys.argv[2]
binding_cite_ = sys.argv[3].split(',')
binding_cite = [float(num) for num in binding_cite_]
search_area_ = sys.argv[4].split(',')
search_area = [float(num2) for num2 in search_area_]
exhaustiveness = sys.argv[5]
data_num = sys.argv[6]
opath = sys.argv[7]

v = Vina(sf_name='vina')
v.set_receptor(denzyme_path)
v.set_ligand_from_file(dligand_path)
v.compute_vina_maps(center=binding_cite, box_size=search_area)
v.dock(exhaustiveness=int(exhaustiveness),n_poses=int(data_num))
v.write_poses(f'{opath}/tmp.pdbqt',n_poses=5,overwrite=True)