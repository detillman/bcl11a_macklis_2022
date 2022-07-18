#!/bin/bash

cd /n/holyscratch01/macklis_lab/dtillman/ambient/working

for FILE in *_filtered.haplotypes.g.vcf.gz; do
TRIMMED=$(echo ${FILE} | cut -f 1 -d '_')
echo ${TRIMMED}
sbatch /n/holyscratch01/macklis_lab/dtillman/ambient/scripts_ambient/14_script_gatk_variantsToTable.run ${TRIMMED}
sleep 1 
done
