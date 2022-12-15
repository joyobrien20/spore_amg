#!/bin/bash

#This was executed on Carbonate interactive job

#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

# local dependencies #
TOOLS=/N/u/danschw/Carbonate/my_tools
MUSCLE=${TOOLS}/muscle5.1.linux_intel64
## muscle 5.1.linux64 [12f0e2]
#TRIMAL=/geode2/home/u020/danschw/Carbonate/my_tools/trimal-trimAl/source
## trimAl v1.4.rev22 build[2015-05-21]. 2009-2015. Salvador Capella-Gutierrez and Toni Gabaldón.
## trimAl webpage: http://trimal.cgenomics.org
## (used in Consuelo Gazitúa wt (2020) paper (Sullivan lab))


##### Define paths #####
PARENT=~/GitHub/spore_amg/genes_of_interest/spo0A_cog
ODIR=${PARENT}/data/align
mkdir -p ${ODIR}


# alignment
$MUSCLE -align ${PARENT}/data/seq2align.faa -stratified -output ${ODIR}/ensemble.efa -threads 8 -amino
$MUSCLE -disperse $ODIR/ensemble.efa
#D_LP=0.02256 D_Cols=0.4727

# If the dispersion is zero, then all the MSAs in the ensemble
# are the same and your alignment is robust. Quite likely,
# it has no errors. If you see a large dispersion, 
# say bigger than 0.05, then there is significant variation 
# between the alignments, and this is necessarily explained 
# by alignment errors.

$MUSCLE -efastats $ODIR/ensemble.efa -log $ODIR/efastats.log

$MUSCLE -maxcc $ODIR/ensemble.efa -output $ODIR/maxcc.afa

## trim alignment
#cd $TRIMAL
#
#./trimal -in ${ODIR}/maxcc.afa  -colnumbering -out ${ODIR}/maxcc.afa.trim -gappyout#
#
## cat ${ODIR}/maxcc.afa.trim
