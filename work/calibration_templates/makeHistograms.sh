#!/bin/bash

time="0900_0930"
for drone in d1 d2 d3 d4 d6 d7 d8; do
    for proj in ${drone}_${time}_results_IDM_SSE_s; do
        analyzeCalibr "$proj.T" -0.1 0.2 20
        analyzeCalibr "$proj.v0" 6 1 20
        analyzeCalibr "$proj.s0" -1 0.5 30
        analyzeCalibr "$proj.a" -0.1 0.2 20
        analyzeCalibr "$proj.b" -0.1 0.2 30
        analyzeCalibr "$proj.GOF" 0 0.2 30
        analyzeCalibr "$proj.GOF" -2.15 0.1 60 1  #log scaling
    done
done
