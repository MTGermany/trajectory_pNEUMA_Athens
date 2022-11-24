
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
set style line 8 lt 1 lw 1 pt 8 ps 1.5 dt 1 lc rgb "#AA00AA"
set style line 9 lt 1 lw 1 pt 9 ps 1.5 dt 2 lc rgb "#999999" #grau
set style line 12 lt 1 lw 3 pt 9 ps 1.5 dt 1 lc rgb "#880000"


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)

# types=[Motorcycle, Car, Medium Vehicle, Heavy Vehicle, Taxi, Bus, redTL]

filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN
filterInterval(data,minVal,maxVal)=((data>=minVal)&&(data<=maxVal)) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN

proj="20181024_d1_0900_0930"

set term post eps enhanced color solid "Helvetica" 12
set param
set key opaque box

set xlabel "t [s]"
set xrange [0:920]

set ylabel "y_{logical} [m]"
set yrange [-5:5]
#set auto y

####################################################
# !!! Drone d1_0900_0930 has y shift by about 0.8 m between t1=525 and t2=575
####################################################

t1=525.
t2=575.


#############################################
laneRef=2
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_yt_demo_yshift.eps", proj, laneRef)
#############################################

set out epsfile
print "plotting ",epsfile


plot\
 infile u (filterInterval($3,0,t1)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, t<%i",t1) w l ls 1,\
 infile u (filterInterval($3,t1,t2)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, %i<=t<%i",t1,t2) w l ls 3,\
 infile u (filterGe($3,t2)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, t>=%i",t2) w l ls 8


#############################################
laneRef=2
infile=sprintf("%s.road%i.traj", proj, laneRef)
epsfile=sprintf("%s_road%i_yt_demo_yshift_detail.eps", proj, laneRef)
#############################################

set out epsfile
print "plotting ",epsfile
set xrange [500:600]

plot\
 infile u (filterInterval($3,0,t1)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, t<%i",t1) w l ls 1,\
 infile u (filterInterval($3,t1,t2)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, %i<=t<%i",t1,t2) w l ls 3,\
 infile u (filterGe($3,t2)*filterCarsTrucks($2)*$3):($5)\
   t sprintf("Cars and trucks, t>=%i",t2) w l ls 8





