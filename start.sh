#!/bin/bash..

script_dir=$(cd "$(dirname "$0")" && pwd)
cd $script_dir
config_name=$script_dir/configs/$1
printf $config_name

printf '\n\n\n\n\n1.create SLC\n'
#python createSLC.py $config_name
printf '\n\n\n\n\n2.create PDG\n'
python createPDG.py $config_name
printf '\n\n\n\n\n3.proc PDG\n'
python procPDG.py $config_name
printf '\n\n\n\n\n4.select PDG\n'
python selectPDG.py $config_name
printf '\n\n\n\n\n5.search OI\n'
python searchOI.py $config_name

