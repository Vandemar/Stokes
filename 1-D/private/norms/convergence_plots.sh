#/bin/bash
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind1_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 1 norm 1 convergence plot'
  plot for [upwind1_1 in system("ls *upwind1*.txt")] upwind1_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind1_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 1 norm 2 convergence plot'
  plot for [upwind1_2 in system("ls *upwind1*.txt")] upwind1_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind1_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 1 norm 3 convergence plot'
  plot for [upwind1_3 in system("ls *upwind1*.txt")] upwind1_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind2_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 2 norm 1 convergence plot'
  plot for [upwind2_1 in system("ls *upwind2*.txt")] upwind2_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind2_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 2 norm 2 convergence plot'
  plot for [upwind2_2 in system("ls *upwind2*.txt")] upwind2_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind2_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 2 norm 3 convergence plot'
  plot for [upwind2_3 in system("ls *upwind2*.txt")] upwind2_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind3_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 3 norm 1 convergence plot'
  plot for [upwind3_1 in system("ls *upwind3*.txt")] upwind3_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind3_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 3 norm 2 convergence plot'
  plot for [upwind3_2 in system("ls *upwind3*.txt")] upwind3_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/upwind3_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'upwind 3 norm 3 convergence plot'
  plot for [upwind3_3 in system("ls *upwind3*.txt")] upwind3_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf1_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 1 norm 1 convergence plot'
  plot for [lf1_1 in system("ls *lf1*.txt")] lf1_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf1_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 1 norm 2 convergence plot'
  plot for [lf1_2 in system("ls *lf1*.txt")] lf1_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf1_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 1 norm 3 convergence plot'
  plot for [lf1_3 in system("ls *lf1*.txt")] lf1_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf2_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 2 norm 1 convergence plot'
  plot for [lf2_1 in system("ls *lf2*.txt")] lf2_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf2_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 2 norm 2 convergence plot'
  plot for [lf2_2 in system("ls *lf2*.txt")] lf2_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf2_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 2 norm 3 convergence plot'
  plot for [lf2_3 in system("ls *lf2*.txt")] lf2_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf3_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 3 norm 1 convergence plot'
  plot for [lf3_1 in system("ls *lf3*.txt")] lf3_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf3_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 3 norm 2 convergence plot'
  plot for [lf3_2 in system("ls *lf3*.txt")] lf3_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lf3_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lf 3 norm 3 convergence plot'
  plot for [lf3_3 in system("ls *lf3*.txt")] lf3_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw1_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 1 norm 1 convergence plot'
  plot for [lw1_1 in system("ls *lw1*.txt")] lw1_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw1_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 1 norm 2 convergence plot'
  plot for [lw1_2 in system("ls *lw1*.txt")] lw1_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw1_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 1 norm 3 convergence plot'
  plot for [lw1_3 in system("ls *lw1*.txt")] lw1_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw2_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 2 norm 1 convergence plot'
  plot for [lw2_1 in system("ls *lw2*.txt")] lw2_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw2_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 2 norm 2 convergence plot'
  plot for [lw2_2 in system("ls *lw2*.txt")] lw2_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw2_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 2 norm 3 convergence plot'
  plot for [lw2_3 in system("ls *lw2*.txt")] lw2_3 using 1:4 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw3_conv1.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 3 norm 1 convergence plot'
  plot for [lw3_1 in system("ls *lw3*.txt")] lw3_1 using 1:2 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw3_conv2.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 3 norm 2 convergence plot'
  plot for [lw3_2 in system("ls *lw3*.txt")] lw3_2 using 1:3 
EOF
gnuplot << EOF
  set term postscript eps color
  set output 'plots/lw3_conv3.eps'
  set logscale xy
  set xlabel 'log(h)'
  set ylabel 'log(norm)'
  set title 'lw 3 norm 3 convergence plot'
  plot for [lw3_3 in system("ls *lw3*.txt")] lw3_3 using 1:4 
EOF
