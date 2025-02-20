set terminal pngcairo size 800,600

# Plot angular velocities
set output 'angular_velocity_omega_x.png'
set title "Angular Velocity (omega_x)"
set xlabel "Time (s)"
set ylabel "Angular Velocity (rad/s)"
set grid
plot 'simulation_data.dat' using 1:6 with lines title 'True omega_x', \
     'simulation_data.dat' using 1:13 with lines title 'Estimated omega_x'

set output 'angular_velocity_omega_y.png'
set title "Angular Velocity (omega_y)"
set xlabel "Time (s)"
set ylabel "Angular Velocity (rad/s)"
set grid
plot 'simulation_data.dat' using 1:7 with lines title 'True omega_y', \
     'simulation_data.dat' using 1:14 with lines title 'Estimated omega_y'

set output 'angular_velocity_omega_z.png'
set title "Angular Velocity (omega_z)"
set xlabel "Time (s)"
set ylabel "Angular Velocity (rad/s)"
set grid
plot 'simulation_data.dat' using 1:8 with lines title 'True omega_z', \
     'simulation_data.dat' using 1:15 with lines title 'Estimated omega_z'

# Plot quaternion components
set output 'quaternion_q_x.png'
set title "Quaternion (q_x)"
set xlabel "Time (s)"
set ylabel "Quaternion Components"
set grid
plot 'simulation_data.dat' using 1:2 with lines title 'True q_x', \
     'simulation_data.dat' using 1:9 with lines title 'Estimated q_x'

set output 'quaternion_q_y.png'
set title "Quaternion (q_y)"
set xlabel "Time (s)"
set ylabel "Quaternion Components"
set grid
plot 'simulation_data.dat' using 1:3 with lines title 'True q_y', \
     'simulation_data.dat' using 1:10 with lines title 'Estimated q_y'

set output 'quaternion_q_z.png'
set title "Quaternion (q_z)"
set xlabel "Time (s)"
set ylabel "Quaternion Components"
set grid
plot 'simulation_data.dat' using 1:4 with lines title 'True q_z', \
     'simulation_data.dat' using 1:11 with lines title 'Estimated q_z'

set output 'quaternion_q_w.png'
set title "Quaternion (q_w)"
set xlabel "Time (s)"
set ylabel "Quaternion Components"
set grid
plot 'simulation_data.dat' using 1:5 with lines title 'True q_w', \
     'simulation_data.dat' using 1:12 with lines title 'Estimated q_w'

# Plot estimation errors and square root of diagonal elements of P
set output 'error_and_sqrt_p_q_x.png'
set title "Estimation Error and sqrt(P) (q_x)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:16 with lines title 'Error q_x', \
     'simulation_data.dat' using 1:23 with lines title 'sqrt_p_0'

set output 'error_and_sqrt_p_q_y.png'
set title "Estimation Error and sqrt(P) (q_y)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:17 with lines title 'Error q_y', \
     'simulation_data.dat' using 1:24 with lines title 'sqrt_p_1'

set output 'error_and_sqrt_p_q_z.png'
set title "Estimation Error and sqrt(P) (q_z)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:18 with lines title 'Error q_z', \
     'simulation_data.dat' using 1:25 with lines title 'sqrt_p_2'

set output 'error_and_sqrt_p_q_w.png'
set title "Estimation Error and sqrt(P) (q_w)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:19 with lines title 'Error q_w', \
     'simulation_data.dat' using 1:26 with lines title 'sqrt_p_3'

set output 'error_and_sqrt_p_omega_x.png'
set title "Estimation Error and sqrt(P) (omega_x)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:20 with lines title 'Error omega_x', \
     'simulation_data.dat' using 1:27 with lines title 'sqrt_p_4'

set output 'error_and_sqrt_p_omega_y.png'
set title "Estimation Error and sqrt(P) (omega_y)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:21 with lines title 'Error omega_y', \
     'simulation_data.dat' using 1:28 with lines title 'sqrt_p_5'

set output 'error_and_sqrt_p_omega_z.png'
set title "Estimation Error and sqrt(P) (omega_z)"
set xlabel "Time (s)"
set ylabel "Value"
set grid
plot 'simulation_data.dat' using 1:22 with lines title 'Error omega_z', \
     'simulation_data.dat' using 1:29 with lines title 'sqrt_p_6'