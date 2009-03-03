set multiplot layout 1,2
set title "Elevation map"
plot '../data/map_2.dat' using 1:2:3   with image not,\
     'state.dat' index i ps 10 not, \
     'cloud.dat' index i pt 3 ps 0.2 not

set title "Elevation measured by the aircraft"
plot 'output.dat' u 1 w l not,\
     "tmp.dat" index i pt 4 not
unset multiplot
     i=i+1
if(i<150) reread

