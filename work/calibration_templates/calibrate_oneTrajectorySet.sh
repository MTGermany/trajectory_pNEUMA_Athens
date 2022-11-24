#!/bin/bash

if (($#!=2)); then
    echo "usage: calibrate_oneTrajectorySet.sh drone time"
    echo "ex: calibrate_oneTrajectorySet.sh d1 0900_0930"
    echo "GOF and model fixed variables within"
    exit -1;
fi
###################################################

drone=$1
time=$2
usedData="${drone}_${time}"

####################################################

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
model=$IDM; str_model="IDM"; GOF=$SSE_s; str_GOF="s"
#model=$IDM; str_model="IDM"; GOF=$SSE_lns; str_GOF="lns"
#model=$ACC; str_model="ACC"; GOF=$SSE_s; str_GOF="s"
#model=$GIP; str_model="GIP"; GOF=$SSE_s; str_GOF="s"
#model=$FVDM; str_model="FVDM"; GOF=$SSE_s; str_GOF="s"
#model=$LCM; str_model="LCM"; GOF=$SSE_s; str_GOF="s"
###################################################


for f in `ls ${usedData}_*.FCdata`; do
    proj=`basename $f .FCdata`
    echo "proj=$proj"
    calibTraj $GOF $model $proj $noLandscape
done

resultBasename="${usedData}_results_${str_model}_SSE_${str_GOF}"
params=""
if( (($model<2)) || (($model==3)) ); then params="v0 T s0 a b"; fi
if(($model==2)); then params="v0 T s0 a b cool"; fi
if(($model==5)); then params="v0 T s0 a gamma"; fi
if(($model==10)); then params="beta0 beta1 beta2 beta3"; fi

# option -P is perl-like regexp. Otherwise, I cannot grep a tab!

for param in $params; do
    echo "summarizing results for $param in $resultBasename.$param"
    #grep -P "\t${param}=" ${usedData}_*_${str_model}_${str_GOF}.out > $resultBasename.$param;
    grep -P "\t${param}=" ${usedData}_*.out > $resultBasename.$param;
done


echo "summarizing results for GOF in $resultBasename.GOF"
grep "resulting_SSE" *.out > $resultBasename.GOF

echo "summarizing combined results in $resultBasename"

#with param errors
#paste -d '\t' $resultBasename.[^G]* $resultBasename.GOF |
#    cut -d $'\t' -f 1,2,3,5,6,8,9,11,12,14,15,18 > $resultBasename

#without param errors, only GOF
if( (($model<2)) || ((model==3)) || ((model==5)) );then
    paste -d '\t' $resultBasename.[^G]* $resultBasename.GOF |
	cut -d $'\t' -f 1,2,5,8,11,14,18 > $resultBasename;
fi
if (($model==2));then
    paste -d '\t' $resultBasename.[^G]* $resultBasename.GOF |
	cut -d $'\t' -f 1,2,5,8,11,14,17,21 > $resultBasename;
fi
if (($model==10));then
    paste -d '\t' $resultBasename.[^G]* $resultBasename.GOF |
	cut -d $'\t' -f 1,2,5,8,11,15 > $resultBasename;
fi

for f in `ls ${usedData}_*.FCdata`; do num=`grep -c -P $'\t6\t' $f`; if(($num>0)); then echo "$f has $num TL lines"; fi; done



echo ""; echo "making parameter lists as input for makeHistograms.sh"

for param in $params "GOF"; do
    echo "manipulating file $resultBasename.$param"
    perl -i -p -e 's/#//g' $resultBasename.$param
    perl -i -p -e 's/://g' $resultBasename.$param
    perl -i -p -e 's/://g' $resultBasename.$param
    perl -i -p -e 's/=/\t/g' $resultBasename.$param
done

perl -i -p -e 's/resulting_SSE.*avg_error/avg_error/g' $resultBasename.GOF
perl -i -p -e 's/m//g' $resultBasename.GOF
