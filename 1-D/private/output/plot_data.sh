#!/bin/bash
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot1_1.eps'
  set xlabel "x"
  set title 'IC 1 h = 1/64'
  set ylabel "u,v"
  plot "upwind1_u1.txt" using 1:2 with lines title 'IC1 upwind1', \
  "lf1_u1.txt" using 1:2 with lines title 'IC1 lf1', \
  "lw1_u1.txt" using 1:2 with lines title 'IC1 lw1', \
  "lf1_u1.txt" using 1:3 with lines title 'IC1 exact'
EOF
epstopdf plots/plot1_1.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot1_3.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 1 h = 1/256'
  plot "upwind1_u3.txt"  using 1:2 with lines title 'IC1 upwind3', \
  "lf1_u3.txt" using 1:2 with lines title 'IC1 lf3', \
  "lw1_u3.txt" using 1:2 with lines title 'IC1 lw3', \
  "lf1_u3.txt" using 1:3 with lines title 'IC1 exact'
EOF
epstopdf plots/plot1_3.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot1_5.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 1 h = 1/1024'
  plot "upwind1_u5.txt" using 1:2 with lines title 'IC1 upwind5', \
  "lf1_u5.txt" using 1:2 with lines title 'IC1 lf5', \
  "lw1_u5.txt" using 1:2 with lines title 'IC1 lw5', \
  "lf1_u5.txt" using 1:3 with lines title 'IC1 exact'
EOF
epstopdf plots/plot1_5.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot1_7.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 1 h = 1/4096'
  plot "upwind1_u7.txt" using 1:2 with lines title 'IC1 upwind7', \
  "lf1_u7.txt" using 1:2 with lines title 'IC1 lf7', \
  "lw1_u7.txt"  using 1:2 with lines title 'IC1 lw7', \
  "upwind1_u7.txt" using 1:3 with lines title 'IC1 exact'
EOF
epstopdf plots/plot1_7.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot1_9.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 1 h = 1/16384'
  plot "upwind1_u9.txt" using 1:2 with lines title 'IC1 upwind9', \
  "lf1_u9.txt" using 1:2 with lines title 'IC1 lf9', \
  "lw1_u9.txt"  using 1:2 with lines title 'IC1 lw9', \
  "upwind1_u9.txt" using 1:3 with lines title 'IC1 exact'
EOF
epstopdf plots/plot1_9.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot2_1.eps'
  set xlabel "x"
  set title 'IC 2 h = 1/64'
  set ylabel "u,v"
  plot "upwind2_u1.txt" using 1:2 with lines title 'IC2 upwind1', \
  "lf2_u1.txt" using 1:2 with lines title 'IC2 lf1', \
  "lw2_u1.txt" using 1:2 with lines title 'IC2 lw1', \
  "lf2_u1.txt" using 1:3 with lines title 'IC2 exact'
EOF
epstopdf plots/plot2_1.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot2_3.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 2 h = 1/256'
  plot "upwind2_u3.txt"  using 1:2 with lines title 'IC2 upwind3', \
  "lf2_u3.txt" using 1:2 with lines title 'IC2 lf3', \
  "lw2_u3.txt" using 1:2 with lines title 'IC2 lw3', \
  "lf2_u3.txt" using 1:3 with lines title 'IC2 exact'
EOF
epstopdf plots/plot2_3.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot2_5.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 2 h = 1/1024'
  plot "upwind2_u5.txt" using 1:2 with lines title 'IC2 upwind5', \
  "lf2_u5.txt" using 1:2 with lines title 'IC2 lf5', \
  "lw2_u5.txt" using 1:2 with lines title 'IC2 lw5', \
  "lf2_u5.txt" using 1:3 with lines title 'IC2 exact'
EOF
epstopdf plots/plot2_5.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot2_7.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 2 h = 1/4096'
  plot "upwind2_u7.txt" using 1:2 with lines title 'IC2 upwind7', \
  "lf2_u7.txt" using 1:2 with lines title 'IC2 lf7', \
  "lw2_u7.txt"  using 1:2 with lines title 'IC2 lw7', \
  "upwind2_u7.txt" using 1:3 with lines title 'IC2 exact'
EOF
epstopdf plots/plot2_7.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot2_9.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 2 h = 1/16384'
  plot "upwind2_u9.txt" using 1:2 with lines title 'IC2 upwind9', \
  "lf2_u9.txt" using 1:2 with lines title 'IC2 lf9', \
  "lw2_u9.txt"  using 1:2 with lines title 'IC2 lw9', \
  "lf2_u9.txt" using 1:3 with lines title 'IC2 exact'
EOF
epstopdf plots/plot2_9.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot3_1.eps'
  set xlabel "x"
  set title 'IC 3 h = 1/64'
  set ylabel "u,v"
  plot "upwind3_u1.txt" using 1:2 with lines title 'IC3 upwind1', \
  "lf3_u1.txt" using 1:2 with lines title 'IC3 lf1', \
  "lw3_u1.txt" using 1:2 with lines title 'IC3 lw1', \
  "lf3_u1.txt" using 1:3 with lines title 'IC3 exact'
EOF
epstopdf plots/plot3_1.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot3_3.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 3 h = 1/256'
  plot "upwind3_u3.txt"  using 1:2 with lines title 'IC3 upwind3', \
  "lf3_u3.txt" using 1:2 with lines title 'IC3 lf3', \
  "lw3_u3.txt" using 1:2 with lines title 'IC3 lw3', \
  "lf3_u3.txt" using 1:3 with lines title 'IC3 exact'
EOF
epstopdf plots/plot3_3.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot3_5.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 3 h = 1/1024'
  plot "upwind3_u5.txt" using 1:2 with lines title 'IC3 upwind5', \
  "lf3_u5.txt" using 1:2 with lines title 'IC3 lf5', \
  "lw3_u5.txt" using 1:2 with lines title 'IC3 lw5', \
  "lf3_u5.txt" using 1:3 with lines title 'IC3 exact'
EOF
epstopdf plots/plot3_5.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot3_7.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 3 h = 1/4096'
  plot "upwind3_u7.txt" using 1:2 with lines title 'IC3 upwind7', \
  "lf3_u7.txt" using 1:2 with lines title 'IC3 lf7', \
  "lw3_u7.txt"  using 1:2 with lines title 'IC3 lw7', \
  "lf3_u7.txt" using 1:3 with lines title 'IC3 exact'
EOF
epstopdf plots/plot3_7.eps
gnuplot << EOF
  set term postscript eps color 
  set output 'plots/plot3_9.eps'
  set xlabel "x"
  set ylabel "u,v"
  set title 'IC 3 h = 1/16384'
  plot "upwind3_u9.txt" using 1:2 with lines title 'IC3 upwind9', \
  "lf3_u9.txt" using 1:2 with lines title 'IC3 lf9', \
  "lw3_u9.txt"  using 1:2 with lines title 'IC3 lw9', \
  "lf3_u9.txt" using 1:3 with lines title 'IC3 exact'
EOF
epstopdf plots/plot3_9.eps
