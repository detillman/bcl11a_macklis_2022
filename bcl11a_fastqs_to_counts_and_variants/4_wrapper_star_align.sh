#!/bin/bash

cd /n/holyscratch01/macklis_lab/dtillman/ambient/working

for FILE in *R1.fq.gz; do
TRIMMED=$(echo ${FILE} | cut -f 1 -d 'R')
echo ${TRIMMED}
sbatch /n/holyscratch01/macklis_lab/dtillman/ambient/scripts_ambient/4_script_star_align.run ${TRIMMED}
sleep 1 
done
