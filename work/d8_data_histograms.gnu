

##############################################################
proj="d8"
#proj2="20181024_d8_0930_1000.road2"
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
infile=sprintf("%s.accDistr",proj)
epsfile=sprintf("%s_accDistr.eps",proj)
##############################################################


set out epsfile
print "plotting ", epsfile
set key opaque box top right

set xlabel "acceleration [m/s^2]"
set xrange [-3:3]
set ylabel "count"
set boxwidth 0.9 relative
plot\
  infile u 1:2 t "" w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50

##############################################################

infile=sprintf("%s.accDistr",proj)
epsfile=sprintf("%s_invaccDistr.eps",proj)
set out epsfile
print "plotting ", epsfile
set xlabel "inverse acceleration [s^2/m]"
set xrange [-3:3]
set ylabel "count*acc^2"
set boxwidth 0.9 relative
plot\
  infile u (1/$1):($2*$1**2) t "" w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50


#infile=sprintf("%s.accDistr",proj2)
#epsfile=sprintf("%s_accDistr.eps",proj2)
#set out epsfile
#print "plotting ", epsfile
#replot

###############################################################

infile=sprintf("%s.TTCdistr",proj)
epsfile=sprintf("%s_TTCdistr.eps",proj)
set out epsfile
print "plotting ", epsfile
set xlabel "TTC [s]"
offset=0.03
set xrange [0.4:8]
set ylabel "count"
set boxwidth 0.5 relative
plot\
  infile u ($1-offset):2 t "v<=v_c" w boxes lc rgb "#0000cc" lw 2 fs solid 0.50,\
  infile u ($1+offset):3 t "v>v_c" w boxes lc rgb "#cc0033" lw 2 fs solid 0.50


###############################################################

infile=sprintf("%s.TTCdistr",proj)
epsfile=sprintf("%s_bkindistr.eps",proj)
set out epsfile
print "plotting ", epsfile
set xlabel "bkin [m/s^2]"
offset=0.01
set xrange [0.2:4]
set ylabel "count"
set boxwidth 0.5 relative
plot\
  infile u ($4-offset):5 t "v<=v_c" w boxes lc rgb "#0000cc" lw 2 fs solid 0.50,\
  infile u ($4+offset):6 t "v>v_c" w boxes lc rgb "#cc0033" lw 2 fs solid 0.50

###############################################################

infile=sprintf("%s.TTCdistr",proj)
epsfile=sprintf("%s_invbkindistr.eps",proj)
set out epsfile
print "plotting ", epsfile
set xlabel "bkin^{-1} [s^2/m]"
offset=0.01
set xrange [0.2:4.5]
set ylabel "count"
set boxwidth 0.5 relative
plot\
  infile u ($7-offset):8 t "v<=v_c" w boxes lc rgb "#0000cc" lw 2 fs solid 0.50,\
  infile u ($7+offset):9 t "v>v_c" w boxes lc rgb "#cc0033" lw 2 fs solid 0.50

###############################################################

infile=sprintf("%s.TTCdistr",proj)
epsfile=sprintf("%s_critIDMdistr.eps",proj)
set out epsfile
print "plotting ", epsfile
set xlabel "critIDM [s^2/m]"
offset=0.01
set xrange [0.2:4.5]
set ylabel "count"
set boxwidth 0.5 relative
plot\
  infile u ($10-offset):11 t "v<=v_c" w boxes lc rgb "#0000cc" lw 2 fs solid 0.50,\
  infile u ($10+offset):12 t "v>v_c" w boxes lc rgb "#cc0033" lw 2 fs solid 0.50

