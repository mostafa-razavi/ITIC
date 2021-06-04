#!/bin/bash
WORK_DIR='/home/mostafa/Git/ITIC'
ZUres_input_data='C8_MiPPE.zures'
ITIC_conditions='C8.itic'
Outname='C8_MiPPE.sat'

$WORK_DIR/src/ITIC $ZUres_input_data $ITIC_conditions | tee "ITIC.log"
cat "ITIC.log" | tail -n7  > $Outname
