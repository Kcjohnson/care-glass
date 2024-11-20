#!/bin/bash
#SBATCH --get-user-env
#SBATCH --job-name=norlux_md5sums_%a
#SBATCH --chdir=/projects/verhaak-lab/USERS/johnsk/glass4
#SBATCH --output=/projects/verhaak-lab/USERS/johnsk/glass4/logs/md5/norlux_batch2_aligned_bams_md5_%a.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kevin.c.johnson@jax.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=03:00:00
#SBATCH --array=1-72
#SBATCH --partition=compute
#SBATCH --qos=batch

# Define the input file and retrieve the filename for this task
INPUT_FILE="/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/GLSS_LX_batch2_realigned_filesize_filterme.txt"
FILEPATH=$(awk "NR==$SLURM_ARRAY_TASK_ID {print \$2}" $INPUT_FILE)
FILENAME=$(basename $FILEPATH)

STARTTIME=`date`
echo $STARTTIME
echo "Analyzing $FILEPATH"
echo ""

# Calculate the md5sum for the file
md5sum $FILEPATH > /projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/aligned_bam_md5/${FILENAME}.md5sum

### END ###
