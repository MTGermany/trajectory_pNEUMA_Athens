
set style line 1 lt 1 lw 2 pt 7 ps 1.9 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 2 pt 9 ps 1.5 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 2 pt 4 ps 1.2 dt 1 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 2 pt 4 ps 1.5 dt 1 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 2 pt 5 ps 1.5 dt 1 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 2 pt 4 ps 1.5 dt 1 lc rgb "#9999FF"
set style line 7 lt 1 lw 2 pt 4 ps 2.0 dt 1 lc rgb "#4477FF" #blau,
set style line 8 lt 1 lw 3 pt 8 ps 1.5 dt 1 lc rgb "#AA00AA"
set style line 9 lt 1 lw 2 pt 9 ps 1.5 dt 1 lc rgb "#999999" #grau
set style line 12 lt 1 lw 4 pt 9 ps 1.5 dt 1 lc rgb "#DD0000"

set term post eps enhanced color solid "Helvetica" 14

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN

proj="20181024_d8_0900_0930"

w=3.2 # as on the real road
lane(y,laneRef)=round(laneRef+y/(0.85*w)) # excluding motorcycles
laneMotoLeft(y,laneRef)=round(laneRef+(y-0.5*w)/(0.5*w))
laneMotoRight(y,laneRef)=round(laneRef+(y+0.5*w)/(0.5*w))


set param
set key opaque box bottom right


set xlabel "Time [s]"



set ylabel "Distance x [m]"
set yrange [0:375]


#############################################
# laneRef=2; lanePlot=2; set xrange [650:900]  # left w/o motorcycles
# laneRef=1; lanePlot=1; set xrange [190:350]   # right w/ motorcycles
# laneRef=3; lanePlot=3; set xrange [650:900] # reverse left w/o motorcycles
 laneRef=4; lanePlot=4; set xrange [650:900] # reverse right w/ motorcycles

infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt_TGF.eps", proj, laneRef, lanePlot)
str_lanePlotAll=sprintf("lane=%i, all",lanePlot)
str_lanePlotMoto=sprintf("lane=%i, motos center",lanePlot)
str_lanePlotMotoLeft=sprintf("lane=%i, motos left",lanePlot)
str_lanePlotMotoRight=sprintf("lane=%i, motos right",lanePlot)
#############################################

#print "infile= ",infile; quit
set out epsfile
print "plotting ",epsfile
plot\
 infile u (filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotAll w l ls 1,\
 infile u \
   (filterData($2,0)*filterData(laneMotoLeft($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotMotoLeft w l ls 5,\
 infile u \
   (filterData($2,0)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotMoto w l ls 9,\
 infile u \
   (filterData($2,0)*filterData(laneMotoRight($5,laneRef),lanePlot)*$3):($4)\
   t str_lanePlotMotoRight w l ls 7,\
 infile u \
   (filterData($2,6)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Red Traffic Lights") w l ls 12

quit

#############################################
laneRef=2
lanePlot24
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt_TGF.eps", proj, laneRef, lanePlot)
str_lanePlotAll=sprintf("lane=%i, all",lanePlot)
str_lanePlotMoto=sprintf("lane=%i, motorcycles",lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile
replot


quit

