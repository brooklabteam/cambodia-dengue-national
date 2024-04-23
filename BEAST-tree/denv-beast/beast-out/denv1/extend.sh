#!/bin/bash
#SBATCH --job-name=denv1
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbrook@rcc.uchicago.edu

module load midway2
module load java/1.8
module load cuda/8.0
module load beagle/trunk
module load beast/2.6.2

/home/cbrook/beast/bin/beast -beagle_CPU -seed 777 -resume /home/cbrook/denv-beast/denv1/denv1.xml

