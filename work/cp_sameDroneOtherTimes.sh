#!/bin/bash

#extraction template. Only the corresp csv files are needed

date=20181024
drone="d1"
workdir="$HOME/trafficSim/data/traj_pNEUMA_Athens/work"


cd $workdir

#for time in 0800_0830 0830_0900 0900_0930 0930_1000 1000_1030 1030_1100; do
for time in 0900_0930; do
    project="${date}_${drone}_${time}"
    if test -f $project.csv; then
        cp template_${drone}.param $project.param
        cp template_${drone}_heatmap.gnu ${project}_heatmap.gnu
        cp template_${drone}_lanes.gnu ${project}_lanes.gnu
        cp template_${drone}_logical_xt.gnu ${project}_logical_xt.gnu
        cp template_${drone}.run $project.run
        cp template_${drone}.plot $project.plot
        perl -i -p -e "s/template_${drone}/${project}/g" ${project}_*.gnu $project.run $project.plot
        echo "copied template project template_${drone} to $project"
    else
	echo "file $project.csv does not exist"
    fi
done 
