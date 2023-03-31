#!/bin/bash
#! comout.sh
#! Script to combine the model output into one file and also tidy up the file space

export BASE_DIR=$HOME/NEMO/NEMO_4.0.4_mirror
export MODEL=IDEAL
export EXPNAM=EXPH10
export NUM_XIOS=2
export NUM_CPU=60

export OUT_FREQ=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $2 }')
export OUT_START=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $3 }')
export OUT_END=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $4 }')

for val in ${OUT_FREQ[@]}; do
	# use rebuild_nemo to combine output files
	$BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_T $NUM_XIOS
	$BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_U $NUM_XIOS
	$BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_V $NUM_XIOS
	$BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_W $NUM_XIOS

	# move files into OUTPUTS folder
	mv ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_T.nc $BASE_DIR/cfgs/$MODEL/$EXPNAM/OUTPUTS
	mv ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_U.nc $BASE_DIR/cfgs/$MODEL/$EXPNAM/OUTPUTS
	mv ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_V.nc $BASE_DIR/cfgs/$MODEL/$EXPNAM/OUTPUTS
	mv ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_W.nc $BASE_DIR/cfgs/$MODEL/$EXPNAM/OUTPUTS

	# delete remaining output files in exp folder
	rm -v ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_T_*.nc
	rm -v ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_U_*.nc 
	rm -v ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_V_*.nc 
	rm -v ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_W_*.nc 
done

# move restart files in to RESTARTS. These files may be combined later so `ln_xios_read=.true.` for XIOS reading by NEMO.
# mv ${MODEL}_*_restart_*.nc $BASE_DIR/cfgs/$MODEL/$EXPNAM/RESTARTS

# tidy up the file space, remove uneeded files
rm -v nam_rebuild*
rm -v submit*.pbs.*
