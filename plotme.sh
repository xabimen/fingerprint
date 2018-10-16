set term pngcairo transparent size 1200,1000 font "Helvetica,20"
set output 'output.png'

set autoscale xfix
set autoscale yfix
set cbrange [0:0.5]
set palette defined (0 "white", 1 "green", 2 "yellow", 3 "orange", 4 "red", 5 "blue", 6 "purple", 7 "black" )

p 'output' matrix w image
