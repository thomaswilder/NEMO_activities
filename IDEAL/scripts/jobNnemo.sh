#!/bin/bash
#
# Script Name: jobNnemo.pbs
#
# Author: T.Wilder
# Date: 20/03/2023
#
# Description: This script submits the job scripts to run nemo and tidy up the file space. 
#
# Run Information: Execute from the Monsoon2 command line using bash jobNnemo.sh

# ----------------------------------------------------------------------------------
# ----------------------------------- Main -----------------------------------------
# ----------------------------------------------------------------------------------

# shell variables
LAST=1 # should one less than the total jobs submitted.
OUT_FREQ=($(seq 1 1 $LAST)) # ($(seq FIRST STEP LAST))

# submit jobs
FIRST=$(qsub jobnemo.pbs)
echo $FIRST
	
# tidy up file and update namelist_cfg
SECOND=$(qsub -W depend=afterok:$FIRST jobnamup.pbs)
echo $SECOND

JOB1=$SECOND

# loop n times to achieve desired spin up length
for val in ${OUT_FREQ[@]}; do
	echo $val

	JOB2=$(qsub -W depend=afterok:$JOB1 jobnemo.pbs)
	echo $JOB2
	
	JOB3=$(qsub -W depend=afterok:$JOB2 jobnamup.pbs)
	echo $JOB3
		
	JOB1=$JOB3
	
done

