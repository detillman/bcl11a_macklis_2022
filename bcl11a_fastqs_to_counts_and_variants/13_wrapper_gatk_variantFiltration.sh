#!/bin/bash

cd /n/holyscratch01/macklis_lab/dtillman/ambient/working

for FILE in *_haplotypes.g.vcf.gz; do
TRIMMED=$(echo ${FILE} | cut -f 1 -d '_')
echo ${TRIMMED}
sbatch /n/holyscratch01/macklis_lab/dtillman/ambient/scripts_ambient/13_script_gatk_variantFiltration.run ${TRIMMED}
sleep 1 
done
