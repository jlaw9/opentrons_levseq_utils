# opentrons_levseq_utils
Repo containing helpful code for processing screening data from opentrons robots along with the sequencing of each well using levseq barcoding

# LevSeq
## Installing
First request a [kestrel account](https://www.nrel.gov/hpc/user-accounts) and connect to [kestrel](https://www.nrel.gov/hpc/system-connection) ([kestrel documentation](https://nrel.github.io/HPC/Documentation/Systems/Kestrel/))
(I use `ssh jlaw@kestrel.hpc.nrel.gov`, but you can also use putty for windows for example). 

Then run the following commands:
``` bash
module load conda mamba
cd /projects/bpms
# not sure if this step is needed
mkdir $USER  
cd $USER

# create your own conda environment
# you can either install it in /projects/bpms (slower startup) or in your home folder
# option 1: install in /projects/bpms/$USER/envs
mkdir envs
# optional: download the packages to /scratch (the home folder has 50GB space per user)
# export CONDA_PKGS_DIRS=/scratch/$USER/conda_pkgs
mamba create --yes --prefix ./envs/levseq python=3.12
conda activate ./envs/levseq

# option 2: install in you home folder (faster startup)
mamba create --yes --name levseq python=3.12
conda activate levseq

# install prerequisites
mamba install --yes -c bioconda -c conda-forge samtools minimap2
mamba install --yes -c conda-forge jupyterlab

# install levseq
# export PIP_CACHE_DIR=/scratch/$USER/pip_cache
pip install levseq

# make a levseq folder to hold your inputs and outputs
mkdir levseq
cd levseq
```

### Now you can test your installation
``` bash
# For some reason, some of the drivers don't get mapped correctly.
# Run this command to update the LD_LIBRARY_PATH variable, and add it to your bash profile (so it's automatically loaded when you login to kestrel):
echo "export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:/home/$USER/.conda-envs/levseq/lib\"" >> ~/.bashrc && source ~/.bashrc
# check that the installation works
levseq --help
# try running an example
mkdir runs
cp -r /projects/bpms/jlaw/projects/other/levseq_example/runs/250416_run ./
cd 250416_run
# run levseq: 
# levseq <output_directory> <input_folder_with_fastq_files> <path/to/ref.csv>
levseq outputs inputs inputs/2025-04-16_GsAdh_ref_levSeq.csv
# after it's finished:
ls outputs
# show the top couple lines of the variants file
head -n 3 outputs/variants.csv
```

## Running LevSeq
After the nanopore sequencing data comes back, we need to run levseq to demultiplex the plates and wells

Running levseq:
- You'll need two files:
  1. The reference csv (containing the references sequence)
  2. The sequencing reads fastq file

Copy the files to kestrel. There are other programs you can use for this. I use `scp` on the terminal like this:
```
scp <ref.csv> <reads.fastq> kestrel.hpc.nrel.gov:/projects/bpms/<user>/levseq/<date>_<run_name>/
```

ssh to kestrel (use a DAV node e.g., kd[1-6] so you can run levseq without submitting a job)
```
ssh kd2.hpc.nrel.gov
module load conda
```

Run levseq on kestrel:
```
cd /projects/bpms/$USER/levseq
conda activate levseq  # or /projects/bpms/$USER/envs/levseq

levseq 250507_run_output/ ./250507_run 250507_run/2025-04-16_GsAdh_ref_levSeq.csv
```

If everything goes well, you'll have a variants.csv file in the output folder with the DNA mutations in each plate and well


## Map the variants to amino acid mutations
``` bash
cd /projects/bpms/$USER/levseq
git clone https://github.com/jlaw9/opentrons_levseq_utils.git
# run this python script to load the levseq outputs and map to AA mutations
python opentrons_levseq_utils/map_variants.py --var_file runs/250416_run/outputs/variants.csv
```

Now copy the files back to you computer
