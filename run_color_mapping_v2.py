import argparse
import subprocess
import os
import mdtraj as md

def run_script_pdb(script, args):
    """ 指定されたスクリプトをサブプロセスとして実行 """
    cmd = ["python", script] + args
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
def run_script_png(script, args):
    """ 指定されたスクリプトをサブプロセスとして実行 """
    cmd = ["/Applications/PyMOL.app/Contents/bin/python3.7", script] + args
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="MD解析のスクリプトチェーン実行")
    
    parser.add_argument("-gro", required=True, help="トポロジー用のgroファイル")
    parser.add_argument("-traj", required=True, help="トラジェクトリファイル (xtc, trr, gro, pdb)")
    parser.add_argument("-skip", required=True, help="出力間隔 (1, 10, 100)")
    parser.add_argument("-output_dir_path", default=None, help="出力フォルダのパス (デフォルト: trajファイルと同じ場所)")
    parser.add_argument("-output_dir_name", default=None, help="出力フォルダの名前 (デフォルト: None)")

    args = parser.parse_args()

    # 出力フォルダが指定されていない場合、trajと同じ場所に作成
    if args.output_dir_path is None:
        args.output_dir_path = os.path.dirname(os.path.abspath(args.traj))
    if args.output_dir_name is not None:
        args.output_dir_name = "_" + args.output_dir_name
        
    script_fname = os.path.basename(__file__)
    script_dir_path = os.path.dirname(os.path.abspath(__file__))
    make_pdb_fname = script_dir_path + "/color_mapping_TK_4.py"
    make_png_fname = script_dir_path + "/pymol_load_TK_6.py"
    
    gro_fname = args.gro
    gro = md.load(gro_fname)
    top = gro.topology
    top_table, top_bonds = top.to_dataframe()
    input_resName = top_table['resName'][0]

    # mdtraj.py の実行
    mdtraj_args = [args.gro, args.traj, args.skip, args.output_dir_path, script_dir_path, args.output_dir_name]
    print(mdtraj_args)
    run_script_pdb(make_pdb_fname, mdtraj_args)

    # pymol.py の実行
    pymol_args = [args.output_dir_path, args.output_dir_path, script_dir_path, input_resName, args.output_dir_name]
    print(pymol_args)
    run_script_png(make_png_fname, pymol_args)

if __name__ == "__main__":
    main()
