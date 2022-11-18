#!/bin/bash

if gnuplot templateData_IDM_s.gnu; 
  then echo "plotting successful"
  else echo "plotting failed"; exit -1;
fi

fig2dev -L pdf templateData_IDM_s.fig templateData_IDM_s.pdf
#fig2dev -L eps templateData_IDM_s.fig templateData_IDM_s.eps
#convert templateData_IDM_s.eps templateData_IDM_s.pdf
#evince  templateData_IDM_s.eps &

echo "produced templateData_IDM_s.pdf"
echo "note:  when copying projects use cpcalib, not cpinp!\":"

