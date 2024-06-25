Snakemake pipeline for conducting ATACseq analysis for BF528: Applications in Translational Bioinformatics.

# Initial Set-up
This pipeline is run in a conda environment. After installing miniconda, running this command will create the environment containing snakemake and some additional python libraries:
```
conda env create -f base_env.yml
```

This environment can then be activated using the following command:
```
conda activate ATACseq
```

# Running the Pipeline
To run this pipeline on an SGE computer cluster, the following command can be used:
```
snakemake -s atacseq.smk --sdm conda --executor cluster-generic --cluster-generic-submit-cmd "qsub -P bf528 -pe omp {threads}" --jobs 8
```
