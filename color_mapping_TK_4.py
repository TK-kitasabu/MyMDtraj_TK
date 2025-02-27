import mdtraj as md
import numpy as np
import json
import os
import sys
import argparse
import concurrent.futures
import gc
import sys

#input values
#gro_fname = "UUZU-2N_1920mols_373k_r.gro"
#trr_fname = "UUZU-2N_1920mols_373k_2_whole.xtc"
#skip = 100 #ps

args = sys.argv
gro_fname = args[1]
trr_fname = args[2]
skip = args[3]
output_dir_path = args[4]
json_file_path = args[5]
output_dir_name = args[6]

gro = md.load(gro_fname)
top = gro.topology
top_table, top_bonds = top.to_dataframe()
input_resName = top_table['resName'][0]

#input_json
#script_fname = os.path.basename(__file__)
#script_dir_path = os.path.dirname(os.path.abspath(__file__))
#json_fname = script_dir_path + "/mol_color_separate.json"
json_fname = f"{json_file_path}/color_mapping_TK.json"
with open(json_fname, 'r') as f:
    json_data = json.load(f)
    
res_json_data = json_data.get(input_resName, {})
op_atom_list = res_json_data.get('select_op_long', [])
atom_name_A = op_atom_list[0]
atom_name_B = op_atom_list[1]



atom_list_A = top.select(" name {} ".format(atom_name_A)) # target atom_name
atom_list_B = top.select(" name {} ".format(atom_name_B))
atom_index = np.union1d(atom_list_A, atom_list_B)
xtc = md.load_xtc(trr_fname, top=gro_fname, stride=skip, atom_indices=atom_index) # stride : skip


# Basic info
steps = xtc.n_frames
nmol = xtc.n_atoms / 2
xyz_axis = 3 # The trajectory is in three dimensions, XYZ

# Calculation of order parameter
xtc_index_A = np.where(np.isin(atom_index, atom_list_A))[0]
xtc_index_B = np.where(np.isin(atom_index, atom_list_B))[0]
xr = xtc.atom_slice(xtc_index_A).xyz - xtc.atom_slice(xtc_index_B).xyz
#del xtc
e = xr / np.linalg.norm(xr, axis=2, keepdims=True) # Normalised
mol_tensor = np.einsum('...jk,...jl->...kl', e, e) # Tensor cal.
Q = (mol_tensor * ( 3 / nmol ) - np.identity ( 3 )) / 2 # Q tensor
eigenvalues, eigenvectors = np.linalg.eig( Q ) # Eig cal.
P_2 = np.max(eigenvalues, axis=1).reshape(-1, 1) # sort of eig ( <P_2> )

#formatted_output = "\n".join(["{:.8f}".format(x) for x in P_2.flatten()]) #output format
#print(formatted_output) #output

'''
color_mapping
'''
print("color_mapping")
trajectory = md.load_xtc(trr_fname, top=gro_fname, stride=skip) #traj_load
# 総フレーム数
num_steps = trajectory.n_frames  # 601 フレームあるなら num_steps=601

eign_index = np.argmax(eigenvalues, axis=1)
director = eigenvectors[np.arange(eigenvectors.shape[0]), :, eign_index] #ダイレクター
cos_theta = np.sum(director[-1, np.newaxis, :] * e, axis=2) #cos = n dot e

sign_list = np.sign(np.mean(cos_theta, axis=1)) #正負のlist
sign_list = sign_list*sign_list
P_1 = sign_list[:, np.newaxis] * cos_theta #修正したOP_1のリスト
#directorは修正していない

mol_atoms_num = trajectory.n_atoms / nmol # 1分子当たりの原子数
b_factor_list = np.repeat(P_1, mol_atoms_num, axis=1)

# 出力ディレクトリ
output_dir = f"{output_dir_path}/color_mapping_strcture{output_dir_name}"
os.makedirs(output_dir, exist_ok=True)

# 各ステップごとに処理

def save_pdb(step):
    save_pdb_fname = f"{output_dir}/step_{step}.pdb"
    trajectory[step].save_pdb(save_pdb_fname, bfactors=b_factor_list[step, :])
    print(f"Saved: {save_pdb_fname}")  # 進行状況を表示

def main():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(save_pdb, range(steps), chunksize=max(1, steps // (os.cpu_count() * 2)))
        
if __name__ == "__main__":
    main()

'''
ToDoList
    読み込み、保存ファイル名の変更機能
    保存先pathへの移動
    出力間隔の調節機能
    pymolへの接続機能
    pymolへのデータ渡し機能
    データ加工機能
    pdbファイルへの並列化処理機能
'''


'''
補足
xtcやtrrファイルなどのバイナリファイルから直接、座標データをloadしていることからoctaveの出力結果と異なる場合ある
ただ、inputファイルの桁数による誤差であることを検証済み
(このスクリプトを使ってバイナリファイルよりも精度の低いgroファイルから諸計算を行なったところ、スクリプトのOrderParameter_v2.mによるP1と値は一致した。符号が全てずれていた)

補助資料
>>> eigenvectors.shape
(601, 3, 3)
>>> eign_index.shape
(601,)
eigenvectors[0,:,eign_index[0]]
eigenvectors[1,:,eign_index[1]]
eigenvectors[2,:,eign_index[2]]
...

関数：np.newaxis
使い方：https://qiita.com/rtok/items/10f803a226892a760d75

OK : director[0,:] * e[0,0,:]
NG : director[0,:] * e[0,:,:]
OK : director[0,np.newaxis,:] * e[0,:,:]
このように、サイズが異なっている場合、計算できない
'''
