#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=gencode_cov2db
#SBATCH --chdir=/projects/verhaak-lab/USERS/johnsk/glass4
#SBATCH --output=/projects/verhaak-lab/USERS/johnsk/glass4/logs/gencode_cov2db/gencode_cov2db_glss_lx_batch2.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kevin.c.johnson@jax.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=36:00:00
#SBATCH --partition=compute
#SBATCH --qos=batch


# Run script to upload gencode coverage to the database.
Rscript R/misc/gencode-coverage2db.R

### END ###
