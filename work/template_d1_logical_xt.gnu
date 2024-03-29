######################################################################
##!!!!! geometrical drone error at t=540 s y shift down by 0.7m!!
######################################################################


set style line 1 lt 1 lw 1 pt 7 ps 1.9 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 1 pt 9 ps 1.5 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 1 pt 4 ps 1.2 dt 1 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 1 pt 5 ps 1.5 dt 1 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#9999FF"
set style line 7 lt 1 lw 1 pt 4 ps 2.0 dt 1 lc rgb "#4477FF" #blau,
set style line 8 lt 1 lw 3 pt 8 ps 1.5 dt 2 lc rgb "#AA00AA"
set style line 9 lt 1 lw 4 pt 9 ps 1.5 dt 1 lc rgb "#999999" #grau
set style line 12 lt 1 lw 3 pt 9 ps 1.5 dt 1 lc rgb "#DD0000"
 

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN

set term post eps enhanced color solid "Helvetica" 12

proj="template_d1"

lane(y,laneRef)=round(laneRef+y/2.8)

set param
set key opaque box


set xlabel "t [s]"
set auto x

set ylabel "x_{logical} [m]"
set auto y



#############################################
laneRef=1 
lanePlot=1
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt.eps", proj, laneRef, lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile


plot\
 infile u (filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("lane=%i, all",lanePlot) w l ls 1,\
 infile u (filterData($2,0)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("lane=%i, motorcycles",lanePlot) w l ls 2,\
 infile u \
   (filterData($2,6)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Red Traffic Lights") w l ls 12


quit

#############################################
lanePlot=1
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt.eps", proj, laneRef, lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile
replot

#############################################
lanePlot=3
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt.eps", proj, laneRef, lanePlot)
#############################################

set out epsfile
print "plotting ",epsfile
replot

#############################################
lanePlot=1
tmin=400
tmax=700
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_lane%i_xt_t%i-%i.eps",\
  proj, laneRef, lanePlot, tmin, tmax)
#############################################

tmin=400
tmax=700
set xrange [tmin:tmax]

replot

quit

# single trajs

vehSelect11=1010  # must be on relevant lane
vehSelect21=1008  # must be on relevant lane
vehSelect31=1004  # must be on relevant lane

 infile u \
   (filterData($1,vehSelect11)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Veh %i", vehSelect11) w l ls 1,\
 infile u \
   (filterData($1,vehSelect21)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Veh %i", vehSelect21) w l ls 8,\
 infile u \
   (filterData($1,vehSelect31)*filterData(lane($5,laneRef),lanePlot)*$3):($4)\
   t sprintf("Veh %i", vehSelect31) w l ls 9,\

plot\
 infile u (filterData($2,1)*$4):($5) t "Cars" w l ls 1,\
 infile u (filterData($2,0)*$4):($5) t "Motorcycles" w l ls 2,\
 infile u (filterData($2,2)*$4):($5) t "MedVehs" w l ls 5,\
 infile u (filterData($2,3)*$4):($5) t "Trucks" w l ls 8,\
 infile u (filterData($2,4)*$4):($5) t "Taxis" w l ls 9,\
 infile u (filterData($2,5)*$4):($5) t "Buses" w l ls 4



