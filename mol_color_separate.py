import mdtraj as md
import numpy as np
import json
import os
import argparse

parser = argparse.ArgumentParser(description='このプログラムの説明（なくてもよい）')
parser.add_argument('arg1', help='この引数の説明（なくてもよい）')
args = parser.parse_args()
input_path = args.arg1

input_dir_path = os.path.dirname(os.path.abspath(input_path))
base_fname = os.path.splitext(os.path.basename(input_path))[0]
save_pdb_fname = input_dir_path + "/" + base_fname + "_MCS.pdb"

script_fname = os.path.basename(__file__)
script_dir_path = os.path.dirname(os.path.abspath(__file__))
json_fname = script_dir_path + "/mol_color_separate.json"

gro = md.load(input_path)
top = gro.topology

table, bonds = top.to_dataframe()
color_list = np.column_stack((table['name'].to_numpy(), np.zeros(top.n_atoms, dtype=float)))

input_resName = table['resName'][0]

with open(json_fname, 'r') as f:
    json_data = json.load(f)
res_json_data = json_data.get(input_resName, {})
atom_names = res_json_data.get('atoms', [])

for i in range(color_list.shape[0]):
    atom_name = color_list[i, 0]
    if atom_name in atom_names:
        color_list[i, 1] += 1  # 同じ原子名があった場合は +1
    else:
        color_list[i, 1] -= 1  # なかった場合は -1

b_factor_list = color_list[:,1]
gro.save_pdb(save_pdb_fname, bfactors=b_factor_list)
