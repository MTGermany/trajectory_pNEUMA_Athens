
set style line 1 lt 3 lw 4 pt 7 ps 1.9  lc rgb "#000000" 
set style line 2 lt 1 lw 4 pt 2 ps 1.5  lc rgb "#CC0022" 
set style line 3 lt 8 lw 4 pt 4 ps 1.2 
set style line 4 lt 6 lw 1 pt 4 ps 1.5  lc rgb "#FFBB00" 
set style line 5 lt 1 lw 4 pt 5 ps 1.5  lc rgb "#00DD22" 
set style line 6 lt 5 lw 4 pt 4 ps 1.5  lc rgb "#00AAAA" 
set style line 7 lt 3 lw 4 pt 4 ps 2.0  lc rgb "#1100AA" 
set style line 8 lt 4 lw 4 pt 8 ps 1.5 
set style line 9 lt 7 lw 4 pt 9 ps 1.5  lc rgb "#999999"


#######################################################
str_data="templateData"
str_data_v=sprintf("%s v", str_data)
str_data_vl=sprintf("%s v_l", str_data)
str_proj="templateData_IDM_s"
str_tseries_s=sprintf("figElem/%s.tseries_s.eps",str_proj)
str_tseries_v=sprintf("figElem/%s.tseries_v.eps",str_proj)
print str_tseries_s
print str_tseries_v
# format $s, %i, %1.1f etc
#######################################################

# plot time series of deviations data/mdoel for final beta values


set size 1,0.55
set nogrid

#columns  tdata, sdata,vdata,vldata,ydata,haty

set term post eps enhanced color dash "Helvetica" 20

#######################################################
set out str_tseries_s
print "plotting ".str_tseries_s
#######################################################

set key top right box
set xlabel ""
set auto x
set ylabel "Gap [m]" offset 1,0
set auto y
#set yrange [0:50]  # necessary because of huge gaps if no leader

s_plot(s)=(s<150) ? s : NaN

plot\
  str_proj.".out" u 1:(s_plot($2)) t str_data w l ls 2,\
  str_proj.".out" u 1:(s_plot($5)) t "Calibrated IDM" w l ls 1



#######################################################
set out str_tseries_v
print "plotting ". str_tseries_v
#######################################################

set size 1,0.62
#set nokey
set xlabel "Time [s]"
set ylabel "Speed [m/s]" offset 1,0
set auto y
set ytics 1
plot\
  str_proj.".out"  u 1:($3) t str_data_v w l ls 2,\
  str_proj.".out"  u 1:($4) t str_data_vl w l ls 4,\
  str_proj.".out"  u 1:($6) t "Calibrated IDM" w l ls 1





#######################################################
# objective function landscape
# fixed scaling from program; min(SSE, 2*SSEopt) plotted to be able to
# use autoscale
#######################################################

str_beta0="Desired speed V_0 [m/s]"
str_beta1="Time gap T [s]"
str_beta2="Minimum gap s_0 [m]"
str_beta3="Max acceleration a [m/s^2]"
str_beta4="Comfortable decel. b [m/s^2]"

str_beta0beta1="figElem/".str_proj.".beta0_beta1.eps"
str_beta0beta2="figElem/".str_proj.".beta0_beta2.eps"
str_beta0beta3="figElem/".str_proj.".beta0_beta3.eps"
str_beta0beta4="figElem/".str_proj.".beta0_beta4.eps"
str_beta1beta2="figElem/".str_proj.".beta1_beta2.eps"
str_beta1beta3="figElem/".str_proj.".beta1_beta3.eps"
str_beta1beta4="figElem/".str_proj.".beta1_beta4.eps"
str_beta2beta3="figElem/".str_proj.".beta2_beta3.eps"
str_beta2beta4="figElem/".str_proj.".beta2_beta4.eps"
str_beta3beta4="figElem/".str_proj.".beta3_beta4.eps"

set term post eps enhanced color solid "Helvetica" 22

unset contour               #no contour lines
set contour surface       #Aktiviert Kontourlinien auf 3D-Flaeche
unset clabel  # dann lauter gleiche Kontourlinien; 
               #       Farbe und Typ mit "w l ls" beim splot-Kommando


#set cntrparam bspline 
set cntrparam levels 15                
#set cntrparam levels incr 0,2,1000
#set cntrparam levels incr zmin,5,zmaxCont    
#set cntrparam levels discrete -8,-1.6,0,1

set nokey
set size 0.7,1
set lmargin at screen 0.18   # bugfix wrong clippings
set rmargin at screen 0.70


set palette defined ( 0 "#dd00ff", 3 "blue", 7 "#00cccc",\
      11 "green", 18 "yellow", 30 "orange", 40 "#ff4422", 55 "#ffaaaa", 100 "white") 


set style line 99 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
unset pm3d                            
set pm3d                              
set pm3d  hidden3d 99          
set view 20,330
set pm3d; set pm3d  map        


set surface  
unset surface 


set nogrid 

#set zrange [:1.1*zmax]
set auto z
#set cbrange [:zmax]
set cbtics 


############################################################
set out str_beta0beta1
print "plotting ".str_beta0beta1
############################################################

set colorbox
set label 1 "SSE" at screen 0.89,1.20

set auto x
set auto y
set autoscale fix
set ytics auto

set xlabel "" offset 0, 0.3;  #unset xtics
#set ylabel str_beta1 offset 0.7, 0; # because bug at offsets
set ylabel ""
set label 2 str_beta1 at screen 0.05,0.5 center rotate by 90

splot  str_proj.".beta0_beta1" u 1:2:($3)   w l ls 99

unset colorbox
unset label 1


############################################################
set out str_beta0beta2
print "plotting ".str_beta0beta2
############################################################

#set ylabel str_beta2 offset -1, 0;   # because bug at offsets
set label 2 str_beta2 at screen 0.05,0.5 center rotate by 90
splot  str_proj.".beta0_beta2" u 1:2:($3) w l ls 99


############################################################
set out str_beta0beta3
print "plotting ".str_beta0beta3
############################################################

#set ylabel str_beta3  offset 0.9, 0;  # because bug at offsets
set label 2 str_beta3 at screen 0.05,0.5 center rotate by 90
splot  str_proj.".beta0_beta3" u 1:2:($3) w l ls 99


############################################################
set out str_beta0beta4
print "plotting ".str_beta0beta4
############################################################

set xlabel str_beta0 offset 0, 0.3; set xtics
#set ylabel str_beta4 offset -1, 0; 
set label 2 str_beta4 at screen 0.05,0.5 center rotate by 90
splot  str_proj.".beta0_beta4" u 1:2:($3) w l ls 99

unset label 2

############################################################
set out str_beta1beta2
print "plotting ".str_beta1beta2
############################################################

#set xlabel str_beta1; 
set xlabel ""; 
set ylabel "";
unset ytics
splot  str_proj.".beta1_beta2" u 1:2:($3) w l ls 99




############################################################
set out str_beta1beta3
print "plotting ".str_beta1beta3
############################################################

set ylabel ""; #unset ytics;
splot  str_proj.".beta1_beta3" u 1:2:($3) w l ls 99


############################################################
set out str_beta1beta4
print "plotting ".str_beta1beta4
############################################################

set xlabel str_beta1; set xtics
set ylabel ""; #unset ytics
splot  str_proj.".beta1_beta4" u 1:2:($3) w l ls 99


############################################################
set out str_beta2beta3
print "plotting ".str_beta2beta3
############################################################

set xlabel "";  #unset xtics
set ylabel "";  #unset ytics
splot  str_proj.".beta2_beta3" u 1:2:($3) w l ls 99


############################################################
set out str_beta2beta4
print "plotting ".str_beta2beta4
############################################################


set xlabel str_beta2; set xtics
set ylabel ""; #unset ytics
splot  str_proj.".beta2_beta4" u 1:2:($3) w l ls 99


############################################################
set out str_beta3beta4
print "plotting ".str_beta3beta4
############################################################

set xlabel str_beta3; set xtics
#set ylabel str_beta4; set ytics
set ylabel "";  #unset ytics

splot  str_proj.".beta3_beta4" u 1:2:($3)   w l ls 99




quit

