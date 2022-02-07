#------------------------------
# model test with modeltest-NG
#------------------------------

#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### local dependencies ####
  TOOLS=/N/u/danschw/Carbonate/my_tools
 
  module load singularity
    # singularity version 3.6.4 loaded

##### Define paths #####
PARENT=~/GitHub/spore_amg/metaG/data/coat
cd $PARENT

ODIR=${PARENT}/align-trim-tree/model_test
mkdir -p $ODIR

##### Get modeltest-NG #####
# using Docker image
  
  # # download image
MSNG="$TOOLS/modeltest"
  # mkdir $MSNG
  # cd MSNG
  # singularity pull docker://nanozoo/modeltest-ng:0.1.6--06cdfc1
  
  # to run from container
MSCMD="singularity exec $MSNG/modeltest-ng_0.1.6--06cdfc1.sif modeltest-ng"
  
##### Run modeltest-NG #####

# got this error on first run:
# modeltest-ng: /opt/conda/conda-bld/modeltest-ng_1572470208140/work/src/model/parameter_branches.cpp:138: virtual double modeltest::ParameterBranches::optimize(modeltest::mt_opt_params_t*, double, double, bool): Assertion `fabs(cur_loglh - loglh) < 1e-6' failed.
# adding "- h ug"
# following https://github.com/ddarriba/modeltest/issues/13 

$MSCMD \
-d aa \
-i $PARENT/align-trim-tree/cotJB_MafftEinsi.trim \
-o $ODIR/model_test \
# -p 8 \
-r 123 \
-t ml -h ug

#                          Model         Score        Weight
# ----------------------------------------------------------
#        BIC             LG+G4+F   129759.5012        0.8799
#        AIC             LG+G4+F   127358.9619        0.6799
#       AICc               LG+G4  1934376.9079        1.0000


