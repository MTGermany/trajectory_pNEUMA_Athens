

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
set term post eps enhanced color solid "Helvetica" 16

projName="d8_IDM_s"
histBaseName="d8_results_IDM_SSE_s"

##############################################################
dataName=sprintf("%s.v0.hist",histBaseName)
epsName=sprintf("%s_hist_v0.eps",projName)
set out epsName
print "plotting ", epsName #" from ",dataName
##############################################################

set key top left

set xlabel "v0 [m/s]"
set auto x
set ylabel "\#CF pairs"
set boxwidth 0.9 relative
plot dataName u 1:2 t ""\
     w boxes lc rgb "#aa4477ff" lw 2 fs solid 0.50

#,\
#  "mixedTraffOutput3.hist" u 1:3 t "Other Vehicles"\
#     w boxes lc rgb "#44000000" lw 3 fs transparent

param="T"
dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot

param="s0"
dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot

param="a"
dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot

param="b"
dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot

param="GOF"
dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot

param="GOF_log"

dataName=sprintf("%s.%s.hist",histBaseName,param)
epsName=sprintf("%s_hist_%s.eps",projName,param)
set out epsName
print "plotting ", epsName
set xlabel sprintf("%s [SI]",param)
replot
