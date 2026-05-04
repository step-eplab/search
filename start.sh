#!/bin/bash..

cd "$(dirname "$0")"

printf '\n\n\n\n\n1.create SLC\n'
#python createSLC.py
printf '\n\n\n\n\n2.create PDG\n'
python createPDG.py
printf '\n\n\n\n\n3.proc PDG\n'
python procPDG.py
printf '\n\n\n\n\n4.select PDG\n'
python selectPDG.py
printf '\n\n\n\n\n5.search OI\n'
python searchOI.py

