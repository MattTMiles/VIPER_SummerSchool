#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --job-name=single_source_vary
#SBATCH --mem=4gb
#SBATCH --tmp=4gb



touch "varying_mass_${1}"

#Initialising cluster modules
module purge
module load anaconda3/5.0.1
module use /apps/users/pulsar/common/modulefiles
module use /apps/users/pulsar/skylake/modulefiles

module load tempo2/18e1bf6-gcc-9.2.0

#Explicitly define the python environment
source activate viper

srun --cpus-per-task=1 --time=01:00:00 --mem=2gb --tmp=4gb ~/.conda/envs/viper/bin/python /fred/oz002/users/mmiles/VIPER_SummerSchool/combine_run.py -data /fred/oz002/users/mmiles/VIPER_SummerSchool/mdc2 -orf hd -results /fred/oz002/users/mmiles/VIPER_SummerSchool/results_varyMass/mass_${1} -mass ${1} -freq 2e-8 -distance 60 -realisations 2

rm -f "varying_mass_${1}"


echo done
