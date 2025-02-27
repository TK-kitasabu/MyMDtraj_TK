#/Applications/PyMOL.app/Contents/bin/python3.7 pymol_multiprocessing_test.py -in_dir directory_path -out_dir directory_path
import os
import glob
import re
import pymol
import concurrent.futures
import gc
import sys
import json

#pymolのpythonを開く必要がある
#import sys
#print(sys.executable)

args = sys.argv
input_dir_path = args[1]
output_dir_path = args[2]
json_file_path = args[3]
input_resName = args[4]
output_dir_name = args[5]

# PDBファイルが保存されているディレクトリ
structure_dir = f"{input_dir_path}/color_mapping_strcture{output_dir_name}"
output_dir = f"{output_dir_path}/color_mapping_images{output_dir_name}"
os.makedirs(output_dir, exist_ok=True)


json_fname = f"{json_file_path}/color_mapping_TK.json"
with open(json_fname, 'r') as f:
    json_data = json.load(f)
res_json_data = json_data.get(input_resName, {})
custom_view = res_json_data.get('custom_view', [])

# 画像の出力設定
IMAGE_WIDTH = 1920   # 画像の幅
IMAGE_HEIGHT = 1080  # 画像の高さ
IMAGE_DPI = 150

# PDB ファイルの検索
pdb_files = glob.glob(f"{structure_dir}/step_*.pdb")
# ソート
pdb_files = sorted(pdb_files, key=lambda x: int(re.search(r'step_(\d+)\.pdb', x).group(1)))  # ファイル名が step_0.pdb, step_1.pdb などの前提
# PyMOL を起動 (GUI 非表示)
def render_pdb(pdb_file):
    pymol.finish_launching(['pymol', '-cq'])
    pymol.cmd.reinitialize()  # PyMOL をリセット
    pymol.cmd.load(pdb_file)  # PDB ファイルを読み込み
    pymol.cmd.set("orthoscopic", "true")      # 正射投影モード
    pymol.cmd.set("ray_shadow", "off")        # 影をオフにする
    pymol.cmd.set("cgo_line_width", 3)
    pymol.cmd.set("stick_h_scale", 1.0)
    pymol.cmd.hide("everything")
    pymol.cmd.show("stick")
    pymol.cmd.show("cell")
    pymol.cmd.spectrum("b", minimum=-1, maximum=1)
    # 指定の視点を適用
    pymol.cmd.set_view(custom_view)
    # 画像保存
    image_filename = os.path.join(output_dir, os.path.splitext(os.path.basename(pdb_file))[0] + ".png")
    pymol.cmd.png(image_filename, width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=IMAGE_DPI, ray=1)
    print(f"Saved: {image_filename}")
    pymol.cmd.delete("all")

for pdb_file in pdb_files:
    render_pdb(pdb_file)
