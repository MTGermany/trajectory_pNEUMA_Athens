

##############################################################
drone="d8"
##############################################################



# Beispiele fuer Funktionen 


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(round(data)==number) ? 1 : NaN
filterDataInverse(data,number)=(round(data)==number) ? NaN : 1

selectRange(x,xmin,xmax)=((x>=xmin) && (x<=xmax)) ? 1 : NaN

##############################################################
set term post eps enhanced color solid "Helvetica" 20
#set term pngcairo enhanced color notransparent crop font "Helvetica,12" #better

##############################################################
infile=sprintf("%s.accDistr",drone)
epsfile=sprintf("%s_accDistr.eps",drone)
##############################################################


set out epsfile
print "plotting ", epsfile
set key top left

set xlabel "acceleration [m/s^2]"
set xrange [-3:3]
set ylabel "count"
set boxwidth 0.9 relative
plot\
  infile u 1:2 t "" w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50
