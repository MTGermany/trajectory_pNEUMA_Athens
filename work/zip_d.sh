#!/bin/bash
#make standardized zip files from contents of d1, ..., d8
if (($#!=1)); then
    echo "calling syntax:zip_d.sh <droneDir>"
    echo "example: zip_d.sh d1"
    exit -1;
fi
drone=$1
zip -r $drone.zip $drone/*.FCdata $drone/*.traj $drone/$xt.eps $drone/*lanes $drone/*lanes.eps

