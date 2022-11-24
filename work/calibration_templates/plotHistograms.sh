#!/bin/bash

time="0900_0930"
for drone in d1 d2 d3 d4 d6 d7 d8; do
    cp template_histograms.gnu plotHistogram.gnu
    perl -i -p -e "s/drone=.*$/drone=\"${drone}\"/g" plotHistogram.gnu
    perl -i -p -e "s/time=.*$/time=\"${time}\"/g" plotHistogram.gnu
    gnuplot plotHistogram.gnu
done
