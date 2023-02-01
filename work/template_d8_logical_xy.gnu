
set style line 1 lt 1 lw 1 pt 7 ps 1.9 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 1 pt 9 ps 1.5 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 1 pt 4 ps 1.2 dt 1 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 1 pt 5 ps 1.5 dt 1 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 1 pt 4 ps 1.5 dt 1 lc rgb "#00AAAA"
set style line 7 lt 1 lw 1 pt 4 ps 2.0 dt 1 lc rgb "#1100FF" #blau,
set style line 8 lt 1 lw 1 pt 8 ps 1.5 dt 1 lc rgb "#9999FF"
set style line 9 lt 1 lw 1 pt 9 ps 1.5 dt 9 lc rgb "#999999" #grau


# types=[Motorcycle, Car, Medium Vehicle, Heavy Vehicle, Taxi, Bus, redTL]


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)

filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,minVal)=(data>=minVal) ? 1 : NaN
filterInterval(data,minVal,maxVal)=((data>=minVal)&&(data<=maxVal)) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN

proj="template_d8"

set term post eps enhanced color solid "Helvetica" 12
set param
set key opaque box

set xlabel "x_{logical} [m]"
set auto x

set ylabel "y_{logical} to the left [m]"
set auto y


#############################################
laneRef=2
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_xy.eps", proj, laneRef)
#############################################

set out epsfile
print "plotting ",epsfile

plot\
 infile u (filterCarsTrucks($2)*$4):($5) t "Cars and trucks" w l ls 1,\
 infile u (filterData($2,0)*$4):($5) t "Motorcycles" w l ls 2


#############################################
laneRef=2
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_xy_TaxisBuses.eps", proj, laneRef)
set out epsfile
print "plotting ",epsfile
#############################################


plot\
 infile u (filterCarsTrucks($2)*$4):($5) t "Cars/Trucks" w l ls 1,\
 infile u (filterData($2,0)*$4):($5) t "Motorcycles" w l ls 2,\
 infile u (filterData($2,4)*$4):($5) t "Taxis" w l ls 4,\
 infile u (filterData($2,5)*$4):($5) t "Buses" w l ls 6


quit

#############################################
laneRef=1
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i.eps", proj, laneRef)
#############################################

set out epsfile
print "plotting ",epsfile
replot


quit









plot\
 infile u (filterData($2,1)*$4):($5) t "Cars" w l ls 1,\
 infile u (filterData($2,0)*$4):($5) t "Motorcycles" w l ls 2,\
 infile u (filterData($2,2)*$4):($5) t "MedVehs" w l ls 5,\
 infile u (filterData($2,3)*$4):($5) t "Trucks" w l ls 8,\
 infile u (filterData($2,4)*$4):($5) t "Taxis" w l ls 9,\
 infile u (filterData($2,5)*$4):($5) t "Buses" w l ls 4



