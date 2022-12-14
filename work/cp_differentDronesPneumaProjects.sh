#!/bin/bash

#extraction template. Only the corresp csv files are needed

date="20181024"
time="0900_0930"
refdrone="d8"
workdir="$HOME/trafficSim/data/traj_pNEUMA_Athens/work"


cd $workdir

for drone in d2 d3 d4 d5 d6 d7 d9 d10; do
    project="${date}_${drone}_${time}"
    if test -f $project.csv; then
     # cp template_${refdrone}.param $project.param # only at beginning!
     # cp template_${refdrone}_heatmap.gnu ${project}_heatmap.gnu
      cp template_${refdrone}_lanes.gnu ${project}_lanes.gnu
      cp template_${refdrone}_logical_xt.gnu ${project}_logical_xt.gnu
      perl -i -p -e "s/template_${refdrone}/${project}/g" ${project}_*.gnu
      #extractTraj_pNEUMA $project.csv 2
      gnuplot ${project}_heatmap.gnu
     # extractTraj_pNEUMA $project.csv 3
     # extractTraj_pNEUMA $project.csv 4 2 0
     # extractTraj_pNEUMA $project.csv 4 4 1
     # extractTraj_pNEUMA $project.csv 6 2
     # extractTraj_pNEUMA $project.csv 6 4
    else
	echo "file $project.csv does not exist"
    fi
done 
