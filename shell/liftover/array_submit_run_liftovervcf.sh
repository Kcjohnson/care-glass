#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=liftovervcf-%a
#SBATCH --chdir=/projects/verhaak-lab/USERS/johnsk/glass4
#SBATCH --output=/projects/verhaak-lab/USERS/johnsk/glass4/logs/mutect2/liftover/liftover-b37tohg38-%a.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kevin.c.johnson@jax.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=00:30:00
#SBATCH --array=1-92
#SBATCH --partition=compute
#SBATCH --qos=batch


# Run merge script.
bash shell/liftover/run_liftovervcf.sh ${SLURM_ARRAY_TASK_ID}

### END ###
