#!/bin/bash -eu

nopram=$(ls *_bfac.pdb)
inputfile_1=${1:-$nopram}
inputfile_2=$(dirname $0)"/bfac_png.json"

#変数確認
echo ""
echo "*** "$(basename $0)" : ""inputfile_1 : "$inputfile_1

/Applications/PyMOL.app/Contents/MacOS/python $(dirname $0)/bfac_png.py $inputfile_1 $inputfile_2

echo "*** end script ***"
