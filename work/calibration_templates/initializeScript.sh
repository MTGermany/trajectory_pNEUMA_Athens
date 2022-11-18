

#!/bin/bash

#extraction and calibration template. Only the corresp csv files are needed


#prepare calibration by simulating in work:

date="20181029"
drone="d8"
workdir="$HOME/trafficSim/data/traj_pNEUMA_Athens/work"


cd $workdir

for time in 0800_0830 0830_0900 0900_0930 0930_1000 1000_1030 1030_1100; do
    project="${date}_${drone}_${time}"
    if test -f $project.csv; then
      cp template_${drone}.param $project.param
      cp template_${drone}_heatmap.gnu ${project}_heatmap.gnu
      cp template_${drone}_lanes.gnu ${project}_lanes.gnu
      cp template_${drone}_logical_xt.gnu ${project}_logical_xt.gnu
      perl -i -p -e "s/20181024_d8_0900_0930/${project}/g" ${project}_*.gnu
      extractTraj_pNEUMA $project.csv 3
      extractTraj_pNEUMA $project.csv 4 2 0
      extractTraj_pNEUMA $project.csv 4 4 1
      extractTraj_pNEUMA $project.csv 6 2
      extractTraj_pNEUMA $project.csv 6 4
    else
	echo "file $project.csv does not exist"
    fi
done

cp ../calibration_saveGnuScripts/* .
perl -i -p -e 's/20181024/20181029/g' *sh

