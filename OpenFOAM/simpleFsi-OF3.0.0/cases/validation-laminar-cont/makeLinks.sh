#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Invalid number of arguments"
    echo "Please, specify which time"
    echo "will be linked to processor"
    echo "directories"
    exit -1
fi

PROCESSORS=`ls -d proce*`

for PROC in $PROCESSORS
do
    cd $PROC/$1
    ln -s ../../$1/weaklyCoupledFsiDict.gz
    cd ../../
done

#
#
#


