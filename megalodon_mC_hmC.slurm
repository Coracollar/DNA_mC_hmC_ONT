#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=5-00:00:00
#SBATCH --account=punim1048
#SBATCH --qos=gpgpuresplat
#SBATCH --partition=gpgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

source ~/anaconda3/etc/profile.d/conda.sh
conda activate megalodonenv 

fast5folder=""
genomammi="esembleGRCm38.mmi"
megalodondirectory=""

source ~/anaconda3/etc/profile.d/conda.sh
conda activate megalodonenv 

megalodon $fast5folder \
--guppy-params "-d /home/coracollar/Programs/rerio/basecall_models/basecall_models/" \
--guppy-config res_dna_r941_min_modbases_5mC_5hmC_v001.cfg \
--outputs basecalls mappings mods \
--reference $genomammi --mod-motif m C 0 \
--devices 0 --processes 24 --guppy-server-path ~/Programs/guppy_GPU_4.2.2/ont-guppy/bin/guppy_basecall_server \
--output-directory $megalodondirectory 
