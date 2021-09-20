 set nokey
 set xlabel "Time (years)"
set ylabel                                                                                "e"
 set title "1"
 set terminal postscript monochrome
 set output "giffvbw.eps"
 plot 'giffv.tmp'  with lines     
