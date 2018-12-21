#!/bin/bash

# This is the slurm version of the remote submission script for bunch a psipred tasks
#
#
#


LOCATION=`realpath $1`
NJOBS=$2
echo $NJOBS
#cd $LOCATION

#SBATCH --job-name=psiPredArray
#SBATCH -p express-mobi
#SBATCH --qos express-mobi
#SBATCH --output=psiPredArray%A_%a.out
#SBATCH --error=psiPredArray%A_%a.err
#SBATCH --array=1-$NJOBS
#SBATCH --time=0:05:00
#SBATCH --workdir=$LOCATION
#SBATCH --ntasks=1
#SBATCH -N 1


######################
# Begin work section #
######################



# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

~/runpsipred -i $LOCATION/peptide_$SLURM_ARRAY_TASK_ID.fasta -d ~/db/uniprot_swissprot_current -r $LOCATION



# Do some work based on the SLURM_ARRAY_TASK_ID
# For example:
# ./my_process $SLURM_ARRAY_TASK_ID
#
# where my_process is you executable