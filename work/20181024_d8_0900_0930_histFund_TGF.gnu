

#####################################################################
# dashed now explicit with dt (implicit over "dashed" + ls no longer works)
# BUG (2021, Ubuntu20): dt does not work with png (there always solid)
######################################################################

# post eps dashed no longer works but dashtype (dt) in ls
# specs! dt 1 = solid
 
set style line 1 lt 1 lw 2 pt 7 ps 1.9 dt 1 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 2 pt 9 ps 1.5 dt 1 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 2 pt 4 ps 1.2 dt 1 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 2 pt 4 ps 1.5 dt 1 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 2 pt 5 ps 1.5 dt 1 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 2 pt 4 ps 1.5 dt 1 lc rgb "#00AAAA"
set style line 7 lt 1 lw 2 pt 4 ps 2.0 dt 7 lc rgb "#4477FF" #blau,
set style line 8 lt 1 lw 2 pt 8 ps 1.5 dt 1 lc rgb "#220088"
set style line 9 lt 1 lw 2 pt 9 ps 1.5 dt 9 lc rgb "#999999" #grau

set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" 
set style line 12 lt 1 lw 6 pt 2 ps 1.5 dt 2 lc rgb "#CC0022" 
set style line 13 lt 8 lw 6 pt 4 ps 1.2  lc rgb "#FF3300"
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA"
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100FF"
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#220088"
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"



##############################################################
# Beispiele fuer Funktionen 
##############################################################

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(round(data)==number) ? 1 : NaN
filterDataInverse(data,number)=(round(data)==number) ? NaN : 1

selectRange(x,xmin,xmax)=((x>=xmin) && (x<=xmax)) ? 1 : NaN

##############################################################
set term post eps enhanced color solid "Helvetica" 20
#set term png notransparent truecolor medium font "Helvetica,12"
#set term pngcairo enhanced color notransparent crop font "Helvetica,12" #better

proj="20181024_d8_0900_0930"
xcenter=300
laneIndex=1

##############################################################
infile=sprintf("%s.road%i_x%i.hist",proj,laneIndex,xcenter)
epsfile=sprintf("../figs/%s_road%i_x%i_hist_TGF.eps",proj,laneIndex,xcenter)
##############################################################


set out epsfile
print "plotting ", epsfile
set key box top right

set xlabel "Distance y to the right [m]"
set xrange [-4:6]
set ylabel "\#Vehicles"
set boxwidth 0.9 relative
plot\
  infile u (-$1):2 t "Motorcycles"\
     w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50,\
  infile u (-$1):($3+$4+$5+$6+$7) t "Other vehicles"\
     w boxes lc rgb "#00000000" lw 3 fs transparent


##############################################################
xdet=300
infile=sprintf("%s.road%i_x%i.fund",proj,laneIndex,xdet)
epsfile=sprintf("../figs/%s_road%i_x%i_fund.eps",proj,laneIndex,xdet)
##############################################################

set out epsfile
print "plotting ", epsfile," from ", infile
set key box top right

set xlabel "Density Q/V [veh./km]"
set xrange[0:800]
#set auto x
set ylabel "Flow [veh./h]"
 
plot\
  infile u (($3==0) ? 0 : $2/$3) : ($2) t sprintf("x=%i",xdet)\
     w p ls 1

xdet=250
infile=sprintf("%s.road%i_x%i.fund",proj,laneIndex,xdet)
epsfile=sprintf("../figs/%s_road%i_x%i_fund.eps",proj,laneIndex,xdet)
set out epsfile
print "plotting ", epsfile
replot

xdet=200
infile=sprintf("%s.road%i_x%i.fund",proj,laneIndex,xdet)
epsfile=sprintf("../figs/%s_road%i_x%i_fund.eps",proj,laneIndex,xdet)
set out epsfile
print "plotting ", epsfile
replot


##############################################################
xdet=300
infile=sprintf("%s.road%i_x%i.fund",proj,laneIndex,xdet)
epsfile=sprintf("../figs/%s_road%i_x%i_fundEdie.eps",proj,laneIndex,xdet)
##############################################################
set out epsfile
print "plotting ", epsfile
plot infile u ($4) : ($5) t sprintf("x=%i",xdet)\
     w p ls 1


quit 

# histograms with separate car and other vehs
##############################################################
infile=sprintf("%s.road%i_x%i.hist",proj,laneIndex,xcenter)
epsfile=sprintf("%s_road%i_x%i_hist_TGF_3classes.eps",proj,laneIndex,xcenter)
##############################################################

set out epsfile
print "plotting ", epsfile
set key top left

set xlabel "Distance y to the right [m]"
set xrange [-4:6]
set ylabel "\#Vehicles"
set boxwidth 0.9 relative
plot\
  infile u (-$1):2 t "Motorcycles"\
     w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50,\
  infile u (-$1):3 t "Cars"\
     w boxes lc rgb "#00000000" lw 3 fs solid 0.40,\
  infile u (-$1):($4+$5+$6+$7) t "Other vehicles"\
     w boxes lc rgb "#99ff0000" lw 3 fs transparent

