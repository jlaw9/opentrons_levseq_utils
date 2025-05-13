# opentrons_levseq_utils
Repo containing helpful code for processing screening data from opentrons robots along with the sequencing of each well using levseq barcoding

## Running LevSeq
I installed levseq in a conda environment on kestrel. I had to make a couple minor changes to the code to get it to run.

Running levseq:
- You'll need two files:
  1. The reference csv (containing the references sequence)
  2. The sequencing reads fastq file
```
# First copy the files to kestrel
scp <ref.csv> <reads.fastq> kd2.hpc.nrel.gov:/projects/bpms/jlaw/projects/other/levseq/250507_run/
```

On kestrel:
```
cd /projects/bpms/jlaw/projects/other/levseq
conda activate /projects/bpms/jlaw/envs/levseq

levseq 250507_run_output/ ./250507_run 250507_run/2025-04-16_GsAdh_ref_levSeq.csv
```

If everything goes well, you'll have a variants.csv file in the output folder with the DNA mutations in from each plate and well

## Mapping to AA mutations and combining with screening data
Run the notebook `250425_map_variants.ipynb`
