#!/bin/bash

#extraction template. Only the corresp csv files are needed
if (($#!=2)); then
    echo "calling syntax: extractAllFCdata.sh <csvname> <drone>"
    echo "example: extractAllFCdata.sh 0181024_d1_0900_0930.csv 1"
    exit -1;
fi

csvname=$1
drone=$2
project=`basename $csvname .csv`
templatedir="$HOME/trafficSim/data/traj_pNEUMA_Athens/work"

if test -f "$csvname";
then echo "$csvname is ok";
else echo "$csvname does not exist"; exit -1;
fi

# copy standard control and info files; renaming F... difficult
# use tmp dir to be able to use mmv

cp $templatedir/vehicleProperties.param .

if test -d tmpdir; then rm -r tmpdir; fi
mkdir tmpdir
cp $templatedir/template_d$drone.* $templatedir/template_d${drone}_* tmpdir
cd tmpdir
mmv "template_d${drone}*" "${project}#1"
perl -i -p -e "s/template_d${drone}/${project}/g" $project.run $project.plot ${project}_*.gnu
cd ..
mv tmpdir/* .
rmdir tmpdir

#deactivate traffic lights and save them by renaming them to *_hidden

for f in `ls $project.trafficLights*[0-9]`; do
    mv $f ${f}_hidden;
done



#rm $project.trafficLights*_hidden
#mmv "$project.trafficLights*" "$project.trafficLights#1_hidden"

#==============================================
echo "creating .FCdata files for project $project by running $project.run"
#comment out making of global hetamaps
perl -i -p -e 's/^(.+)\.csv 2$/#\1\.csv 2/g' $project.run 
$project.run
#==============================================

for f in `ls ${project}_lanes.gnu ${project}_logical_xt.gnu`; do
    gnuplot $f;
done

exit 0

############################################################################
