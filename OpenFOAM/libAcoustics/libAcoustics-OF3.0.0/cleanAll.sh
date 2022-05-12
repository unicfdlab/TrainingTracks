#!/bin/bash

#
# clean cases
#
cd cases
CASES_DIRS=`ls -d ?*`
for CASE in $CASES_DIRS
do
    cd $CASE
    ./Allclean
    cd ../
done
cd ../

#
# clean sources
#

cd lib; wclean; cd ../
cd "test/SimpleFoamFunctionObject"; wclean; cd ../../
cd "test/SimpleFoamLibrary"; wclean; cd ../../


#
#END-OF-FILE
#


