
set style line 1 lt 1 lw 1 pt 7 ps 1.9 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 1 pt 9 ps 1.5 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 1 pt 4 ps 1.2 dt 1 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 1 pt 5 ps 1.5 dt 1 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#9999FF"
set style line 7 lt 1 lw 1 pt 4 ps 2.0 dt 2 lc rgb "#1100FF" #blau,
set style line 8 lt 1 lw 3 pt 8 ps 1.5 dt 1 lc rgb "#AA00AA"
set style line 9 lt 1 lw 4 pt 9 ps 1.5 dt 1 lc rgb "#999999" #grau
set style line 12 lt 1 lw 4 pt 9 ps 1.5 dt 1 lc rgb "#880000"

set term post eps enhanced color solid "Helvetica" 12

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN
proj="20181024_d8_0900_0930"


lane(y,laneRef)=round(laneRef+y/3.1)

set param
set key opaque box
#set size 0.77,1
#set size square
#set size ratio -1

set xlabel "t [s]"
#set xrange [150:400]
set auto x


set ylabel "x_{logical} [m]"
#set yrange [0:400]
set auto y



#############################################
laneRef=4   #!! only lanePlot<=2 goes to the right
lanePlot=4
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt.eps", proj, laneRef, lanePlot)
str_lanePlotAll=sprintf("lane=%i, all",lanePlot)
str_lanePlotMoto=sprintf("lane=%i, motorcycles",lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile
plot\
 infile u (filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotAll w l ls 1,\
 infile u (filterData($2,0)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotMoto w l ls 2,\
 infile u \
   (filterData($2,6)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Red Traffic Lights", vehSelect31) w l ls 12



# infile u \
#   (filterData($1,vehSelect11)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
#   t sprintf("Veh %i", vehSelect11) w l ls 1,\
# infile u \
#   (filterData($1,vehSelect21)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
#   t sprintf("Veh %i", vehSelect21) w l ls 8,\

quit

#############################################
laneRef=4   #!! only lanePlot<=2 goes to the right
lanePlot=3
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt.eps", proj, laneRef, lanePlot)
str_lanePlotAll=sprintf("lane=%i, all",lanePlot)
str_lanePlotMoto=sprintf("lane=%i, motorcycles",lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile
replot
