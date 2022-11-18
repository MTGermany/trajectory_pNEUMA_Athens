#!/bin/bash

if(ls 20181024* > tmp 2> tmp); then
    echo "removing prefix 20181024_";
    mmv "20181024_*" "#1";
fi

IDM=0
ACC=2
GIP=3
FVDM=5
LCM=10

SSE_s=1
SSE_lns=2
SSE_a=0
noLandscape=0
landscape=1

###################################################
usedData="d1_0900_0930"
model=$IDM;  str_model="IDM";  str_calib="IDM_s"
#model=$IDM;  str_model="IDM";  str_calib="IDM_lns"
#model=$ACC;  str_model="ACC";  str_calib="ACC_s"
#model=$GIP;  str_model="GIP";  str_calib="GIP_s"
#model=$FVDM; str_model="FVDM"; str_calib="FVDM_s"
#model=$LCM;  str_model="LCM";  str_calib="LCM_s"
###################################################


for f in "d8_0900_0930_veh848.FCdata"; do   #test
#for f in "d1_0900_0930_veh952.FCdata"; do   #test
#for f in `ls ${usedData}_*.FCdata`; do
    proj=`basename $f .FCdata`
    fullProj="${proj}_${str_calib}"
    echo "proj=$proj, fullProj=$fullProj"
    cpcalib templateData $proj $str_calib
    $fullProj.run
done

