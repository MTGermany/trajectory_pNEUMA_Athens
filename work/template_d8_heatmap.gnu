max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN


proj="template_d8"


set param
set key opaque box


# auto x has sometimes
# gnuplot bugs with multiplot but good for general scale findung

set xlabel "x [m]"
#set ylabel "y to the left [m]"
set ylabel "y [m]"

set surface; unset pm3d; set pm3d map


zVal(z)=(z>0.1) ? z : NaN;
toNorth2(heading)=(heading>=0) ? 1 : NaN;
toSouth2(heading)=(heading<0) ? 1 : NaN;
toNorth3(heading)=(abs(heading)>=0.5*pi) ? 1 : NaN;
toSouth3(heading)=(abs(heading)<0.5*pi) ? 1 : NaN;


################ WhatToDo=3 #################
infile=sprintf("%s.heatmap3",proj)
epsfile=sprintf("%s_heatmap3.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

set cbrange [0:50]  

# replot geht nicht mit multiplot!
# bei WhatToDo=3 verwende ich aber anyway nur die toSouth-Schattierung

set term post eps enhanced color solid "Helvetica" 16
set label 1 "Northern directions" at screen 0.59,0.815 front textcolor rgbcolor "#00aa33"
set label 2 "Southern directions" at screen 0.59,0.850 front textcolor rgbcolor "#ff2200"

set size noratio
set size 0.90,1

set xrange [-550:-150]
set yrange [1155:1185]

set multiplot
set palette defined ( 0 "white", 15 "green", 40 "blue", \
      100 "#000088")
splot infile u ($1):($2):(toNorth3($3)*zVal($4))  t "" w p palette ps 0.30
set palette defined ( 0 "white", 10 "yellow", 30 "orange", \
      50 "red", 100 "#880000")
splot infile u ($1):($2):(toSouth3($3)*zVal($4))  t "" w p palette ps 0.30
unset multiplot

#quit

################ WhatToDo=2 #################
infile=sprintf("%s.heatmap2",proj)   # unrotated, with WhatToDo=2
epsfile=sprintf("%s_heatmap2.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

#set size square
set size ratio -1
set term post eps enhanced color solid "Helvetica" 12


set xrange [-450:100]    # for unrotated setting
set yrange [1000:1380]   # for unrotated setting

set multiplot
set palette defined ( 0 "white", 5 "green", 30 "blue", \
      100 "#000088")
splot infile u ($1):($2):(toNorth2($3)*zVal($4))  t "" w p palette ps 0.10

set palette defined ( 0 "white", 2 "yellow", 10 "orange", \
      30 "red", 100 "#880000")
splot infile u ($1):($2):(toSouth2($3)*zVal($4))  t "" w p palette ps 0.10
unset multiplot




