#!/bin/bash
cd scf || exit
vacuum=$(echo -e "426\n3\n"|vaspkit|grep "^ Vacuum-Level"|cut -d":" -f2|awk '$1=$1')
cd ../band || exit
gap=$(cat BAND_GAP |grep "Band Gap"|cut -d: -f2|awk '$1=$1')
fermi=$(cat BAND_GAP |grep "Fermi Energy"|cut -d: -f2|awk '$1=$1')
vbm=$(echo $(echo "$vacuum-($fermi)")|bc)
cbm=$(echo $(echo "$vbm-$gap")|bc)
echo -e "vbm:-$vbm\ncbm:-$cbm\n"


