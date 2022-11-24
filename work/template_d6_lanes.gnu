max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN

set term post eps enhanced color solid "Helvetica" 12
proj="template_d6"

#############################################
infile=sprintf("%s.peaks",proj)
epsfile=sprintf("%s_peaks.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

set param
set key opaque box

set xlabel "x [m]"
set auto x

set ylabel "y [m]"
set auto y

set surface; unset pm3d; set pm3d map
set auto cb
set cbrange [0:50]  

zVal(z)=(z>0.1) ? z : NaN;
toNorth(heading)=(heading>=0) ? 1 : NaN;
toSouth(heading)=(heading<0) ? 1 : NaN;


set palette defined ( 0 "white", 2 "yellow", 10 "orange", \
      30 "red", 100 "#550000")
splot infile u ($1):($2):(zVal($4))  t "" w p pt 7 ps 0.5 palette




#############################################
infile=sprintf("%s.lanes",proj)
epsfile=sprintf("%s_lanes.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

unset label 1
unset label 2

set palette defined ( 0 "white", 5 "yellow", 20 "orange", \
      50 "red", 100 "#550000")
splot infile u ($2):($3):(zVal($5))  t "" w l  palette lw 3





