#!/bin/bash

cd /n/holyscratch01/macklis_lab/dtillman/ambient/working

for FILE in *_Aligned.sortedByCoord.out.bam; do
TRIMMED=$(echo ${FILE} | cut -f 1 -d 'A')
echo ${TRIMMED}
sbatch /n/holyscratch01/macklis_lab/dtillman/ambient/scripts_ambient/5_script_picard_markDuplicates.run ${TRIMMED}
sleep 1 
done
