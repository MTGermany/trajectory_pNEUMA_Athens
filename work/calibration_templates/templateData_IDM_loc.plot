if gnuplot templateData_IDM_loc.gnu; 
  then echo "plotting successful"
  else echo "plotting failed"; exit -1;
fi

fig2dev -L pdf templateData_IDM_loc.fig templateData_IDM_loc.pdf
#fig2dev -L eps templateData_IDM_loc.fig templateData_IDM_loc.eps
#convert templateData_IDM_loc.eps templateData_IDM_loc.pdf
#evince  templateData_IDM_loc.eps &

echo "produced templateData_IDM_loc.pdf"
echo "note:  when copying projects use cpcalib, not cpinp!\":"

