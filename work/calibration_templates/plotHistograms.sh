#!/bin/bash

time="0900_0930"
model="IDM"
GOF="s"
#for drone in d6; do
for drone in d1 d2 d3 d4 d6 d7 d8; do

    histogramProject="${drone}_${time}_${model}_${GOF}"
    cp template_histograms.gnu plotHistogram.gnu
    perl -i -p -e "s/drone=.*$/drone=\"${drone}\"/g" plotHistogram.gnu
    perl -i -p -e "s/time=.*$/time=\"${time}\"/g" plotHistogram.gnu
    gnuplot plotHistogram.gnu

    cp templateDataset_${model}_${GOF}_hist.fig ${histogramProject}_hist.fig
    perl -i -p -e "s/template_IDM_s/${histogramProject}/g" \
         ${histogramProject}_hist.fig
    fig2eps ${histogramProject}_hist.fig
    echo "made ${histogramProject}_hist.eps"
    
done
