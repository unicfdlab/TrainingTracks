#!/bin/bash
#
#

C0=100
FREQ=186.047
INFILE="complexPressureData/medium_surfacePressure/$FREQ/dataWithMesh.msh"
python3 bem/pulsosphere_imported_data.py $C0 $FREQ $INFILE

#
#END-OF-FILE
#

