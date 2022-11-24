#!/bin/bash

time="0900_0930"
for drone in d1 d2 d3 d4 d6 d7 d8; do
    calibrate_oneTrajectorySet.sh $drone $time
done
