#!/bin/bash

#extraction template. Only the corresp csv files are needed

if (($#!=3)); then
    echo "usage: cp_PNEUMAproject.sh drone date time"
    echo "example: cp_PNEUMAproject.sh d8 20181024 0900_0930"
    exit -1;
fi

drone=$1
date=$2
time=$3

workdir="$HOME/trafficSim/data/traj_pNEUMA_Athens/work"


cd $workdir
project="${date}_${drone}_${time}"
echo "drone=${drone}"
echo "time=${time}"
if test -f $project.csv; then
    cp template_${drone}.param $project.param
    cp template_${drone}.run $project.run
    cp template_${drone}.plot $project.plot
    for f in `ls template_${drone}.*`; do
	base=${f%.*}
	ext=${f#$base.}
	#cp -v $f ${project}.$ext
	cp $f ${project}.$ext
    done
    
    if test -f template_${drone}_logical_xy.gnu; then
	cp template_${drone}_logical_xy.gnu ${project}_logical_xy.gnu
    fi
    if test -f template_${drone}_logical_yt.gnu; then
	cp template_${drone}_logical_yt.gnu ${project}_logical_yt.gnu
    fi
    cp template_${drone}_heatmap.gnu ${project}_heatmap.gnu
    cp template_${drone}_lanes.gnu ${project}_lanes.gnu
    cp template_${drone}_logical_xt.gnu ${project}_logical_xt.gnu
    
    for f in `ls ${project}_*.gnu $project.run $project.plot`; do
	perl -i -p -e "s/template/${date}/g" $f
	perl -i -p -e "s/${drone}/${drone}_${time}/g" $f
    done
    echo "copied project template_${drone} to $project"
    
else
    echo "file $project.csv does not exist"
    
fi

