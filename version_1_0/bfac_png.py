import os
import sys
import json
from pymol import cmd

pdb_file = sys.argv[1]
# ファイルを読み込む
cmd.load(pdb_file)
# モデル内のすべてのオブジェクト名を取得する
objects_name_list = cmd.get_object_list()
# 最後に読み込まれたオブジェクト名
last_loaded_object_name = objects_name_list[-1] if objects_name_list else None
# モデル内のすべての残基名を取得する
residue_list = set()
cmd.iterate(last_loaded_object_name, 'residue_list.add(resn)', space=locals())
#print(residue_list)
# リストを表示
for residue_list in residue_list:
    print("")




json_open = open(sys.argv[2], 'r')
json_load = json.load(json_open)

b_factor_min = json_load[residue_list]['b_factor_min']
b_factor_max = json_load[residue_list]['b_factor_max']
turn_angle_x = json_load[residue_list]['turn_angle_x']
turn_angle_y = json_load[residue_list]['turn_angle_y']
turn_angle_z = json_load[residue_list]['turn_angle_z']

# Bファクタースペクトルを設定
cmd.spectrum("b", minimum=b_factor_min, maximum=b_factor_max)
# X軸を中心に90度回転
cmd.turn("x", turn_angle_x)
# Y軸を中心に90度回転
cmd.turn("y", turn_angle_y)
# Z軸を中心に90度回転
cmd.turn("z", turn_angle_z)
# 正射投影を設定
cmd.set("orthoscopic", 1)
# 全ての原子をスティック表示
cmd.show("sticks", "all")
# 水素原子を非表示
cmd.hide("everything", "hydrogens")
#描画:300DPI
cmd.ray(width=2559,height=2143)
#pngにて指定ファイル名で保存
cmd.png(last_loaded_object_name + ".png")






#iterate (resi 1), print(name + " %.2f" % b)
