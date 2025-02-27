#!/bin/bash
#オーダーパラメータを作成するスクリプトです。
#指示はmanual.txtに記入
#function
function exit_trap() {
    if [ $? != 0 ]; then
        echo "\n\n"
        echo "Command [$BASH_COMMAND] is failed"
        eval echo "Command [$BASH_COMMAND] is failed"
        echo "\n\n"
        cd ../;
        #rm -r -f "${dirName}";
        exit 1
    fi
}
function path_cd() {
    local fname_f=${1:-noPrameter}
    if [[ "$fname_f" == noPrameter ]]; then
        read -p "drag & drop send folder or file >> " fname_f #ドラッグ＆ドロップによりフルパスの取得
    fi
    #ファイルかフォルダーか判断
    if [[ "$(basename $fname_f)" == *.* ]]; then
        input_type="file"
        cd $(dirname $fname_f);#入力ファイルがあるフォルダーに移動
        echo ""
        echo "  "$(pwd) "に移動しました。 input_type: "${input_type}
        echo "  -------------------------------------------------------------"
        echo ""
    else
        input_type="folder"
        cd $fname_f; #入力フォルダに移動
        echo ""
        echo "  +++sh/scp/path_cd"
        echo "  "$(pwd) "に移動しました。 input_type: "${input_type}
        echo "  -------------------------------------------------------------"
        echo ""
    fi
}
function Make_idx() {
    #echo "どの原子でindexを作りますか？"
    read -p "atom_a: " atom_a
    read -p "atom_b: " atom_b
    echo "a ${atom_a} | a ${atom_b}\nq" | gmx_d make_ndx -f ${targetFname_tpr}.tpr -o target_atoms_idx0.ndx
    cat target_atoms_idx0.ndx | awk -v pattern="${atom_a}_${atom_b}" '$0 ~ pattern { found=1 } found { print; }' > target_atoms_idx.ndx  #目的のインデックス以外削除, [a-b]を残す。
   
    for i in `cat target_atoms_idx.ndx | tail -n +2 | xargs`; do echo $i; done > idx.txt
    echo "index selected Atom" > idx_atom.txt
    echo "atom_a : ${atom_a}\natom_b : ${atom_b}" >> idx_atom.txt
}
trap exit_trap ERR
run_path=$(cd $(dirname $0); cd ../../; pwd)

opm_fname=`pwd`"/orderparameter_ver-8-0.m"
fpwd=`pwd`
path_cd #changeDir
echo " m-version: 8-0"
#opDirName=`date +"%Y%m%d"`
#cheakDIR=`ls ./analysis | awk -v dname="op_"${opDirName} '{if($0==dname){print "yes"}else{print "no"}}'`
#ls ./analysis | awk -v dname="op_"${opDirName} 'BEGIN{vname=dname"v[0-9]"}{if($0 ~ vname){a=0}else if($0==dname){print "yes"}else{print "no"};}'
#if [ ${cheakDIR} = "yes" ];then cheakDIR="error";exit; elif [ ${cheakDIR} = "no" ];then mkdir ./analysis/op_${opDirName};cd ./analysis/op_${opDirName};fi;
#mkdir -p ./analysis/op_${opDirName};cd ./analysis/op_${opDirName}
titleName="op"
dirName=${titleName}"_"`date +"%Y%m%d"`
if [ -d "./analysis/$dirName" ];then #dirがあるかどうか
    dirNum=$(ls -d ./analysis/${titleName}_* | sort | awk -F "/" '{print $NF}' | tail -n 1 | awk -F "_" '{print $NF}')
    if [ -d "./analysis/${dirName}_${dirNum}" ];then #dir_数字があるかどうか
        dirNum_new=$(($dirNum+1))
        dirName+="_"$dirNum_new
    else
        dirName+="_1"
        echo "hogehoge"
    fi
else
    dirNum=$(ls -d ./analysis/${titleName}* | sort | awk -F "/" '{print $NF}' | tail -n 1 | awk -F "_" '{print $NF}')
    if [ -d "./analysis/${dirName}_${dirNum}" ];then
        dirNum_new=$(($dirNum+1))
        dirName+="_"$dirNum_new
    fi
fi
mkdir -p ./analysis/${dirName};cd ./analysis/${dirName}
pwd

targetFname=`ls ../../*.trr`
targetFname=${targetFname%.*}
echo ${targetFname}
targetFname_tpr=${targetFname}
if [ ! -e "${targetFname}.tpr" ];then #dirがあるかどうか
    targetFname_tpr="$(dirname ${targetFname})/$(basename -s "_2" ${targetFname})$(basename $(cd ../../ && pwd) | sed -e 's/.*k//' )"
fi
g96fname="trj0.g96"
read -p "skip: " skip
read -p "start_time: " start_time
read -p "end_time: " end_time

Make_idx >& log.txt
echo "      start trjconv"
temp=$(basename `cd ../../;pwd` | tr -d md)
#gmx_d trjconv -f ${targetFname}.trr -s ${targetFname}.tpr -o ${g96fname} -n target_atoms_idx.ndx -novel -pbc whole -skip ${skip} -tu ns -b ${start_time} -e ${end_time} >& log.txt
#time gmx_d trjconv -f ${targetFname}.trr -s ../../*${temp}.tpr -o ${g96fname} -n target_atoms_idx.ndx -novel -pbc whole -skip ${skip} >& log.txt
time gmx_d trjconv -f ${targetFname}.trr -s ../../*${temp}.tpr -o ${g96fname} -n target_atoms_idx.ndx -novel -pbc whole -skip ${skip} >& log.txt #全時間

#gmx_d trjconv -f ${targetFname}.xtc -s ${targetFname}.tpr -o ${g96fname} -n target_atoms_idx.ndx -novel -pbc whole -skip ${skip} -tu ns -b ${start_time} -e ${end_time} >& log.txt
echo "      end   trjconv"
read -p "this_dir: " this_dir
echo "      start c-object"
time ${this_dir}/op #C言語での実行ファイル
echo "      end   c-object"
echo "      start octave"
#cat $g96fname | gsed -n "/TIMESTEP/,/END/p" | awk '{if(NF==2)print $2/1000}' > op_time.xvg
cat "$g96fname" | sed -n '/TIMESTEP/,/END/p' | awk '{if(NF==2)print $2/1000}' > op_time.xvg
time octave -qf ${opm_fname} > result_op${temp}.xvg
cp "${fpwd}/P2_t_inst_NULL.xlsx" ./
echo "      end   octave"
paste -d "," op_time.xvg result_op${temp}.xvg > op.csv
paste -d " " op_time.xvg result_op${temp}.xvg > op.xvg


#gnuplot
#xvgファイルをpng出力

. ${run_path}/xvg_gnuplot.sh op.xvg "\"op\"" "\"time (ns)\"" "\"Order parameter\""
. ${run_path}/expanded.sh $(basename $0)
echo "      end   total"
