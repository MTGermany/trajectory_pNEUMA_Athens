max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterData(data,number)=(data==number) ? 1 : NaN
filterGe(data,number)=(data>=number) ? 1 : NaN

set term post eps enhanced color solid "Helvetica" 12
proj="20181024_d10_0900_0930"


set param
set key opaque box


# auto x has sometimes
# gnuplot bugs with multiplot but good for general scale findung

set xlabel "x [m]"
set xrange [-400:0] # for unrotated setting


set ylabel "y to the left [m]"
set yrange [1000:1400]   # for unrotated setting

set surface; unset pm3d; set pm3d map


zVal(z)=(z>0.1) ? z : NaN;
toNorth(heading)=(heading>=0) ? 1 : NaN;
toSouth(heading)=(heading<0) ? 1 : NaN;

set label 1 "Northern directions" at screen 0.59,0.825 front textcolor rgbcolor "#00aa33"
set label 2 "Southern directions" at screen 0.59,0.860 front textcolor rgbcolor "#ff2200"


################ WhatToDo=3 #################
infile=sprintf("%s.heatmap3",proj)
epsfile=sprintf("%s_heatmap3.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

set cbrange [0:50]  

# replot geht nicht mit multiplot!
# bei WhatToDo=3 verwende ich aber anyway nur die toSouth-Schattierung

set size noratio
set size 1,1

set auto x
set auto y

set palette defined ( 0 "white", 15 "yellow", 40 "orange", \
      70 "red", 100 "#880000")
 plot infile u ($1):($2):(zVal($4))  t "" w p palette ps 0.3


#quit

################ WhatToDo=2 #################
infile=sprintf("%s.heatmap2",proj)   # unrotated, with WhatToDo=2
epsfile=sprintf("%s_heatmap2.eps",proj)
set out epsfile
print "plotting ",epsfile
#############################################

#set size square
set size ratio -1


set xrange [-950:-450]    # for unrotated setting d10
set yrange [250:600]   # for unrotated setting d10

set multiplot
set palette defined ( 0 "white", 5 "green", 30 "blue", \
      100 "#000088")
splot infile u ($1):($2):(toNorth($3)*zVal($4))  t "" w p palette ps 0.10

set palette defined ( 0 "white", 2 "yellow", 10 "orange", \
      30 "red", 100 "#880000")
splot infile u ($1):($2):(toSouth($3)*zVal($4))  t "" w p palette ps 0.10
unset multiplot






#quit

