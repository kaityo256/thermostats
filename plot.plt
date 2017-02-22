set term png
set out "data.png"
T=1.0
set style data line
set xla "exp(-E)"
set yla "P(E)"
set xtics 0.2
set ytics 0.2
p "nose_hoover.dat" u (exp(-$1/T)):2 t "NH" lw 3\
, "langevin.dat" u (exp(-$1/T)):2 t "Langevin" lw 3 \
, "kinetic_moments.dat" u (exp(-$1/T)):2 t "KM" lw 3\
, "nose_hoover_chain.dat" u (exp(-$1/T)):2 t "NHC" lw 3\

set xla "p"
set yla "q"
set xtics 1.0
set ytics 1.0

set style data points
set size square

set out "nose_hoover_ps.png"
p "nose_hoover_ps.dat" u 2:3 t "NH"

set out "langevin_ps.png"
p "langevin_ps.dat" u 2:3 t "Langevin"

set out "kinetic_moments_ps.png"
p "kinetic_moments_ps.dat" u 2:3 t "KM"

set out "nose_hoover_chain_ps.png"
p "nose_hoover_chain_ps.dat" u 2:3 t "NHC"
