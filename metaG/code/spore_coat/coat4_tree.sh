#--------------------------
# phylogeny with  RAxML-NG
#--------------------------


#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

module load raxmlng
  # Flex 2.6.4 loaded.
  # CMake version 3.8.0 loaded.
  # raxmlng version 0.9.0-pthreads loaded.

##### Define paths #####
PARENT=~/GitHub/spore_amg/metaG/data/coat
ODIR=${PARENT}/align-trim-tree/check_msa
mkdir -p $ODIR
ALN=$PARENT/align-trim-tree/cotJB_MafftEinsi.trim
CODEDIR=~/GitHub/spore_amg/metaG/code/spore_coat
#modeltest-ng selected model : 
MODEL="LG+G4+F"
cd $ODIR

raxml-ng --check --msa $ALN  \
--model $MODEL --data-type AA \
--prefix check-msa 

# [00:00:00] Loaded alignment with 476 taxa and 88 sites
# 
# WARNING: Sequences phage_gvd_6 and phage_blastp_12 are exactly identical!
# WARNING: Sequences phage_gvd_8 and phage_gvd_14 are exactly identical!
# WARNING: Sequences phage_gvd_8 and phage_gvd_24 are exactly identical!
# WARNING: Sequences phage_gvd_17 and phage_gvd_21 are exactly identical!
# WARNING: Sequences phage_gvd_20 and bacteria_uniref_365 are exactly identical!
# WARNING: Sequences phage_gvd_26 and phage_blastp_4 are exactly identical!
# WARNING: Sequences phage_gvd_33 and phage_gvd_35 are exactly identical!
# WARNING: Sequences phage_gvd_37 and bacteria_uniref_180 are exactly identical!
# WARNING: Duplicate sequences found: 8
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/spore_amg/metaG/data/coat/align-trim-tree/check_msa/check-msa.raxml.reduced.phy


ALN=$ODIR/check-msa.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model $MODEL --prefix parse-msa 

# Binary MSA file created: parse-msa.raxml.rba
ALN=$ODIR/parse-msa.raxml.rba

# * Estimated memory requirements                : 51 MB
# * Recommended number of threads / MPI processes: 2

# get tree
ODIR=$PARENT/align-trim-tree/tree
mkdir -p $ODIR
cd $ODIR
# raxml-ng --msa $ALN --model $MODEL --threads 4 --seed 123 --prefix $ODIR/MafftEinsi.trim

# parallelized ML search
mkdir -p $ODIR/ml_search
cd $ODIR/ml_search

#-- Vars in batch_raxML-ng.sh
# SEED=$1
# THREADS=$2
cpus=2
# TREES=$3 # type{number}
# MODEL=$4
# ALN=$5
# ODIR=$6 

# 100  ML searches starting with parsimony trees (10 runs x 10 per run)
seeds=($(seq 11 10 101))
for i in ${seeds[@]}; do
sbatch --job-name=MLprs$i --time=5:30:00 --cpus-per-task=$cpus \
$CODEDIR/batch_raxML-ng.sh $i $cpus "pars{10}" $MODEL $ALN "$ODIR/ml_search"
done

# 100  ML searches starting with random trees (10 runs x 10 per run)
seeds=($(seq 13 10 103))
for i in ${seeds[@]}; do
sbatch --job-name=MLrnd$i --time=5:30:00 --cpus-per-task=$cpus \
$CODEDIR/batch_raxML-ng.sh $i $cpus "rand{10}" $MODEL $ALN "$ODIR/ml_search"
done


grep "Final LogLikelihood:" *.raxml.log | sort -k 3 
cat *.raxml.mlTrees > mltrees
raxml-ng --rfdist --tree mltrees --redo --prefix RF
# Loaded 200 trees with 661 taxa.
# Average absolute RF distance in this tree set: 289.487739
# Average relative RF distance in this tree set: 0.219975
# Number of unique topologies in this tree set: 200
#does not seem to converge$CODEDIR/batch_raxML-ng.sh $i $cpus "rand{10}" $MODEL $ALN "$ODIR/ml_search"
done



# Adding more trees to reach 500
# 150  ML searches starting with parsimony trees (10 runs x 10 per run)
seeds=($(seq 311 10 451))
for i in ${seeds[@]}; do
sbatch --job-name=MLpars$i --time=5:59:00 --cpus-per-task=4 \
$CODEDIR/batch_raxML-ng.sh  $i 4 "pars{10}" $MODEL $ALN "$ODIR/ml_search"
done

# 150  ML searches starting with random trees (10 runs x 10 per run)
seeds=($(seq 313 10 453))
for i in ${seeds[@]}; do
sbatch --job-name=MLrand$i --time=5:59:00 --cpus-per-task=4 \
$CODEDIR/batch_raxML-ng.sh  $i 4 "rand{10}" $MODEL $ALN "$ODIR/ml_search"
done

grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head
cat *.raxml.mlTrees > mltrees
raxml-ng --rfdist --tree mltrees --redo --prefix RF
  # Reading input trees from file: mltrees
  # Loaded 500 trees with 468 taxa.
  # 
  # Average absolute RF distance in this tree set: 414.395206
  # Average relative RF distance in this tree set: 0.445586
  # Number of unique topologies in this tree set: 500
#does not seem to converge

best_trees=($( grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head -n 5 | cut -c-8| tr -d . ))

best_trees=( "${best_trees[@]/%/.raxml.bestTree}" )
cat ${best_trees[@]} > mltrees_top
raxml-ng --rfdist --tree mltrees_top --redo --prefix RF_top
  # Reading input trees from file: mltrees_top
  # Loaded 5 trees with 468 taxa.
  # 
  # Average absolute RF distance in this tree set: 343.800000
  # Average relative RF distance in this tree set: 0.369677
  # Number of unique topologies in this tree set: 5
    
# using Best scoring ML tree
BEST=( $(grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head -n 1 | grep -o ".*.raxml")) 
cat "$BEST.bestTree" > ../bestScoreML.tree

# # #reruns to complete time out
# i=103
# sbatch --job-name=MLpars$i --time=2:59:00 --cpus-per-task=2 \
# $CODEDIR/batch_raxML-ng.sh $i 2 "pars{10}" $MODEL $ALN "$ODIR/ml_search"
# # 
# i=413 # 443
# sbatch --job-name=MLrand$i --time=1:59:00 --cpus-per-task=2 \
# $CODEDIR/batch_raxML-ng.sh $i 2 "rand{10}" $MODEL $ALN "$ODIR/ml_search"
# 
# cat *.raxml.mlTrees | wc -l
