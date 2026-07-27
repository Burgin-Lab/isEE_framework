#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --gres=gpu
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --time={{ walltime }}

set -e
export CUDA_MPS_PIPE_DIRECTORY={{mps_dir}}/nvidia-mps_{{ mps }}
export CUDA_MPS_LOG_DIRECTORY={{mps_dir}}/nvidia-log_{{ mps }}
nvidia-cuda-mps-control -d

export CUDA_MPS_PIPE_DIRECTORY={{mps_dir}}/nvidia-mps_{{ mps }}

export CUDA_VISIBLE_DEVICES=0

(cp {{ prmtop }} {{ prmtop }}{{ degen }}; if [ -f {{ prmtop }}{{ degen }}.bak ];  then mv {{ prmtop }}{{ degen }}.bak {{ prmtop }}{{ degen }}; fi; if [ -f {{ inpcrd }}{{ degen }}_min.rst7 ];  then rm {{ inpcrd }}{{ degen }}_min.rst7; fi; pmemd -O -i {{min_inp}} -o {{ name }}_min.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }} -r {{ inpcrd }}{{ degen }}_min.rst7 -ref {{ inpcrd }}; if [ ! -f {{ inpcrd }}{{ degen }}_min.rst7 ];  then exit 1; fi; cp {{ prmtop }}{{ degen }} {{ prmtop }}{{ degen }}.bak; DIFF=$(python {{isee_titrate}} {{ inpcrd }}{{ degen }}_min.rst7 {{ prmtop }}{{ degen }} | tail -1); if [ "$DIFF" != "False" ] && [ "$DIFF" != "True" ];  then exit 1; elif [ "$DIFF" != "False" ];  then mv {{ inpcrd }}{{ degen }}_min.rst7 {{ inpcrd }}{{ degen }}_min_titrated.rst7;  pmemd -O -i {{min_inp}} -o {{ name }}_min.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_min_titrated.rst7 -r {{ inpcrd }}{{ degen }}_min.rst7 -ref {{ inpcrd }}{{ degen }}_min_titrated.rst7;  rm {{ inpcrd }}{{ degen }}_min_titrated.rst7; fi; pmemd.cuda -O -i {{heat_inp}} -o {{ name }}_heat.out -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_min.rst7 -r {{ inpcrd }}{{ degen }}_heat.rst7 -x {{ inpcrd }}{{ degen }}_heat.nc -ref {{ inpcrd }}{{ degen }}_min.rst7; cp {{ prmtop }}{{ degen }} {{ prmtop }}{{ degen }}_${SLURM_JOB_ID}; pmemd.cuda -O -i {{ inp }} -o {{ out }} -p {{ prmtop }}{{ degen }} -c {{ inpcrd }}{{ degen }}_heat.rst7 -r {{ rst }} -x {{ nc }} -ref {{ inpcrd }}{{ degen }}_heat.rst7; if [ -f {{ rst }} ];  then rm {{ prmtop }}{{ degen }}.bak; fi)

echo quit | nvidia-cuda-mps-control