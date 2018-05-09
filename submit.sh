#! /bin/bash
#-----------------------------------------------------
#$ -S /bin/bash                       # Working Shell
#$ -o $JOB_NAME_$JOB_ID.out              # Log file
#$ -j y                               # Write error in log file
#$ -l h_rt=1:00:00                  # Request resource: hard run time hours:minutes:seconds
#$ -l h_vmem=1G                      # Request resource: memory requirements/per slot
#$ -N test                      # Set job name
#$ -cwd                               # Change into directory where you wrote qsub
#$ -binding linear:1                  # activate extractSequential(process_gauges)
#----------------------------------------------------

# a single core instruction job
# time python extract.py

# extract 29 basins from an given input list using 12 cores off one single node
#time seq 0 29 | parallel -j 12 python extract.py -n {}

