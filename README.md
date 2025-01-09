# ASD gene Bcl11a regulates subcellular RNA localization, associative circuitry, and social behavior
RNA-seq and variant calling of CPN somata and GCs that are wild-type, Bcl11a-het, or Bcl11a-null.

GEO contains all raw (fastq inputs), processed (HaplotypeCaller/FeatureCount outputs), and metadata files: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211405

To turn raw files into processed files, use the scripts in bcl11a_fastqs_to_counts_and_variants, following the instructions in workflow.txt.
Expected Run Time: 1-2 days

An R Markdown and miscellaneous files for data analysis and visualization of processed files is in bcl11a_analysis_and_viz.
Expected Run Time: 1-2 hours

Refer to our preprint for additional information: https://doi.org/10.1101/2022.10.06.511159

Software Dependencies:
Trimmomatic (> v0.39)
STAR (> v2.70e)
FeatureCounts (> v1.5.1)
GATK Haplotype Caller (> v4.2.5.0)
DESeq2 (> v1.38.0)
ClusterProfiler (> v4.6.0)
