#!/bin/bash
#
# makerootvisuals.sh - execute as root user to rebuild rootvisuals_C.so and dependencies
#
# author: richard.t.jones at uconn.edu
# version: july 21, 2026
#
# usage: (as root on nod30)
#   # cd /nfs/direct/packages/gpfs/gpfs1/osgusers/spotfinder
#   # bash makerootvisuals.sh

# make sure that . is writable
touch x || exit 1
/bin/rm x

sudo -u gluex bash -c '
source setup.sh
rm -f CobremsGeneration_cc.so CobremsGeneration_cc.d CobremsGeneration_cc_ACLiC_dict_rdict.pcm
rm -f Couples_C.so Couples_C.d Couples_C_ACLiC_dict_rdict.pcm
rm -f Map2D_cc.so Map2D_cc.d Map2D_cc_ACLiC_dict_rdict.pcm
rm -f rootvisuals_C.so rootvisuals_C.d rootvisuals_C_ACLiC_dict_rdict.pcm
$ROOTSYS/bin/root -l -b -q makerootvisuals.C
'

if [[ $? = 0 ]]; then
    echo "Fresh libraries built, now copy to deployment area as:"
    echo "cp CobremsGeneration_cc* Couples_C* Map2D_cc* rootvisuals_C* /var/www/html/tools/spotfinder"
fi
