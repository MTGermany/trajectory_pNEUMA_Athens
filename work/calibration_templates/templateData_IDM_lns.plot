#!/bin/bash

if gnuplot templateData_IDM_lns.gnu; 
  then echo "plotting successful"
  else echo "plotting failed"; exit -1;
fi

fig2dev -L pdf templateData_IDM_lns.fig templateData_IDM_lns.pdf

echo "produced templateData_IDM_lns.pdf"
echo "note:  when copying projects use cpcalib, not cpinp!\":"


