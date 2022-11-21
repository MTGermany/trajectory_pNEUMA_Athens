#!/bin/bash

#for proj in d1_0900_0930_results_IDM_SSE_s d8_0900_0930_results_IDM_SSE_s; do
for proj in d8_results_IDM_SSE_s; do
    analyzeCalibr "$proj.T" -0.1 0.2 20
    analyzeCalibr "$proj.v0" 6 1 20
    analyzeCalibr "$proj.s0" -1 0.5 30
    analyzeCalibr "$proj.a" -0.1 0.2 20
    analyzeCalibr "$proj.b" -0.1 0.2 30
    analyzeCalibr "$proj.GOF" 0 0.2 30
    analyzeCalibr "$proj.GOF" -2.15 0.1 60 1  #log scaling
done
