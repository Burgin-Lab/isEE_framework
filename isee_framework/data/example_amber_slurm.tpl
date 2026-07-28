#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes={{ nodes }}
#SBATCH --gres=gpu
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node={{ ppn }}
#SBATCH --time={{ walltime }}

cp {{ prmtop }} {{ prmtop }}{{ degen }}; if [ -f {{ prmtop }}{{ degen }}.bak ];  then mv {{ prmtop }}{{ degen }}.bak {{ prmtop }}{{ degen }}; fi; if [ -f {{ inpcrd }}{{ degen }}_min.rst7 ];  then rm {{ inpcrd }}{{ degen }}_min.rst7; fi; {{ solver }} -O -i {{min_inp}} -o {{ name }}_min.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }} -r {{ inpcrd }}{{ degen }}_min.rst7 -ref {{ inpcrd }}; if [ ! -f {{ inpcrd }}{{ degen }}_min.rst7 ];  then exit 1; fi; cp {{ prmtop }}{{ degen }} {{ prmtop }}{{ degen }}.bak; DIFF=$(python {{isee_titrate}} {{ inpcrd }}{{ degen }}_min.rst7 {{ prmtop }}{{ degen }} | tail -1); if [ "$DIFF" != "False" ] && [ "$DIFF" != "True" ];  then exit 1; elif [ "$DIFF" != "False" ];  then mv {{ inpcrd }}{{ degen }}_min.rst7 {{ inpcrd }}{{ degen }}_min_titrated.rst7;  {{ solver }} -O -i {{min_inp}} -o {{ name }}_min.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_min_titrated.rst7 -r {{ inpcrd }}{{ degen }}_min.rst7 -ref {{ inpcrd }}{{ degen }}_min_titrated.rst7;  rm {{ inpcrd }}{{ degen }}_min_titrated.rst7; fi; {{ solver }} -O -i {{heat_inp}} -o {{ name }}_heat.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_min.rst7 -r {{ inpcrd }}{{ degen }}_heat.rst7 -x {{ inpcrd }}{{ degen }}_heat.nc -ref {{ inpcrd }}{{ degen }}_min.rst7; cp {{ prmtop }}{{ degen }} {{ prmtop }}{{ degen }}_${SLURM_JOB_ID}; {{ solver }} -O -i {{ inp }} -o {{ out }} -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_heat.rst7 -r {{ rst }} -x {{ nc }} -ref {{ inpcrd }}{{ degen }}_heat.rst7; if [ -f {{ rst }} ];  then rm {{ prmtop }}{{ degen }}.bak; fi