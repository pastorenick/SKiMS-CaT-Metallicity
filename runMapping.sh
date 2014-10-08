#!/bin/bash
# specify the queue name
#PBS -q sstar
# resource requests
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -l walltime=01:00:00:00
source /home/npastore/.bashrc
# run process
module unload python
module load python/2.7.2
python /lustre/projects/p040_swin/data/Nicola/REDUCTION/SKiMS_CaTmet/runMapping.py

	

