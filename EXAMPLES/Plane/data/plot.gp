set term gif animate  opt delay 10 size 600,250  x000000 xffffff

set output 'plane.gif'
set noborder
set notics
set palette gray
set key horiz
i = 0;
unset colorbox
! awk 'BEGIN{i=0;}{print i,$1,"\n\n"; i++;}' output.dat>tmp.dat

load 'loop.gp'
! rm tmp.dat
set output
