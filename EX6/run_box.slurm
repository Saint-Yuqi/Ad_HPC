#!/bin/bash
#SBATCH -C mc
#SBATCH --job-name=box_fft
#SBATCH --output=box_fft_%j.out
#SBATCH --error=box_fft_%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=uzh8

# Load required modules
module load cpeGNU cray-fftw Boost

# Set number of threads for OpenMP and FFTW
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Paths
EXEC=./assign
DATADIR=$SCRATCH/boxdata
RESOLUTIONS=("B100" "B200" "B500" "B1000")

# Loop through all resolutions
for RES in "${RESOLUTIONS[@]}"; do
    echo "Processing $RES..."
    srun $EXEC $DATADIR/$RES 100 > ${RES}.out
done
