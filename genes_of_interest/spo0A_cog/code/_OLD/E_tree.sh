#--------------------------
# phylogeny with  IQtree
#--------------------------


#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

# local dependencies #
TOOLS=/N/u/danschw/Carbonate/my_tools
IQTREE=${TOOLS}/iqtree-2.1.3-Linux/bin/iqtree2


##### Define paths #####
PARENT=~/GitHub/spore_amg/genes_of_interest/spo0A_cog
ALN=$PARENT/data/align-trim-tree/maxcc.afa.trim

$IQTREE -s $ALN --alrt 1000 -B 1000 -T AUTO  --threads-max 8 
