#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --output=Dsp_Test_2# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=nsher012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="refined_dsp_test_2_circular"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load singularity
export SINGULARITY_NV=1
module load centos

centos.sh "module load extra; module load GCC/4.9.3-2.25; module load cuda/9.1; ./virus-model -dt=0.001 Data_structure_circular_refined.xml"
