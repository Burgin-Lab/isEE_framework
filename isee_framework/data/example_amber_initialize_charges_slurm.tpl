#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes={{ nodes }}
#SBATCH --gres=gpu
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node={{ ppn }}
#SBATCH --time={{ walltime }}

{{ solver }} -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }} -ref {{ inpcrd }}