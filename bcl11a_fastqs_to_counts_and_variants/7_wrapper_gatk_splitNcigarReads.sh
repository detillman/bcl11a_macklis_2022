#!/bin/bash

cd /n/holyscratch01/macklis_lab/dtillman/ambient/working

for FILE in *_markedDupes.Aligned.sortedByCoord.out.bam; do
TRIMMED=$(echo ${FILE} | cut -f 1 -d '_')
echo ${TRIMMED}
sbatch /n/holyscratch01/macklis_lab/dtillman/ambient/scripts_ambient/7_script_gatk_splitNcigarReads.run ${TRIMMED}
sleep 1 
done
