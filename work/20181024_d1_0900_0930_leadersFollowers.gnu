######################################################################
##!!!!! geometrical drone error at t=540 s y shift down by 0.7m!!
######################################################################


set style line 1 lt 1 lw 6 pt 7 ps 2.0 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#FF5500" #closedSquare
set style line 4 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#FFBB00" #gelb,
set style line 5 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#22FF00" #gruen,
set style line 6 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#009999"
set style line 7 lt 1 lw 2 pt 7 ps 0.6 dt 1 lc rgb "#0000FF" #blau,
set style line 8 lt 1 lw 3 pt 7 ps 0.6 dt 2 lc rgb "#AA00AA"
set style line 9 lt 1 lw 4 pt 9 ps 0.6 dt 1 lc rgb "#999999" #grau
set style line 12 lt 1 lw 3 pt 9 ps 0.6 dt 1 lc rgb "#DD0000"
 

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN
filterCarsTrucks(data)=( (data==1)||(data==2)||(data==3)) ? 1 : NaN

set term post eps enhanced color solid "Helvetica" 12

proj="20181024_d1_0900_0930"

lane(y,laneRef)=round(laneRef+y/2.8)

set param
set key opaque box
#set size 0.77,1
#set size square
#set size ratio -1





###########################################################
#vehID=986
vehID=952
#vehID=975
#vehID=835
#vehID=839
#vehID=842
infileL=sprintf("%s_veh%i.leaders", proj, vehID)
infileF=sprintf("%s_veh%i.followers", proj, vehID)
###########################################################



###########################################################
#epsfile=sprintf("%s_veh%i_surroundingVehs_xt.eps", proj, vehID)
epsfile=sprintf("%s_surroundingVehs_xt.eps", proj)
set out epsfile
print "plotting ",epsfile
###########################################################

set xlabel "t [s]"
set auto x

set ylabel "distance headway [m]"
set auto y


plot\
 infileL u ($1):(0.1*$3) t sprintf("0.1*x subject veh %i",vehID) w l ls 1,\
 infileL u ($1):($7)  t "gap leader 1"    w p ls 2,\
 infileL u ($1):($11) t "gap leader 2"    w p ls 3,\
 infileL u ($1):($15) t "gap leader 3"    w p ls 4,\
 infileF u ($1):($7)  t "gap follower 1"  w p ls 5,\
 infileF u ($1):($11) t "gap follower 2"  w p ls 6,\
 infileF u ($1):($15) t "gap follower 3"  w p ls 7


###########################################################
#epsfile=sprintf("%s_veh%i_surroundingVehs_yt.eps", proj, vehID)
epsfile=sprintf("%s_surroundingVehs_yt.eps", proj)
set out epsfile
print "plotting ",epsfile
###########################################################

set ylabel "lateral offset [m]"
set auto y


plot\
 infileL u ($1):($4) t sprintf("y subject veh %i",vehID) w l ls 1,\
 infileL u ($1):($8)  t "gap leader 1"    w p ls 2,\
 infileL u ($1):($12) t "gap leader 2"    w p ls 3,\
 infileL u ($1):($16) t "gap leader 3"    w p ls 4,\
 infileF u ($1):($8)  t "gap follower 1"  w p ls 5,\
 infileF u ($1):($12) t "gap follower 2"  w p ls 6,\
 infileF u ($1):($16) t "gap follower 3"  w p ls 7



###########################################################
epsfile=sprintf("%s_surroundingVehs_xy.eps", proj)
set out epsfile
print "plotting ",epsfile
###########################################################

set xlabel "distance headway [m]"
set auto x
set ylabel "lateral offset [m]"
set auto y


plot\
 infileL u (0.0001*$1):(0.00001*$3) t sprintf("Subject veh %i",vehID) w p ls 1,\
 infileL u ($7):($8)  t "gap leader 1"    w p ls 2,\
 infileL u ($11):($12) t "gap leader 2"    w p ls 3,\
 infileL u ($15):($16) t "gap leader 3"    w p ls 4,\
 infileF u ($7):($8)  t "gap follower 1"  w p ls 5,\
 infileF u ($11):($12) t "gap follower 2"  w p ls 6,\
 infileF u ($15):($16) t "gap follower 3"  w p ls 7
