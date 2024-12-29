# Descriptions for the codes !!!!!!!!!!!


# For the detail of this Brownian particle simulation code please read the manual file named brownian_cluster_dynamics_simulation_code_manual.pdf and also Vikki_Ph.D_thesis.pdf



code: bcd_hard_spheroids.c
input file: runner_bhs.sh

Description: This code simulate the system of hard brownian spheroids and calculate the dynamics, structure following the formulation of brownian cluster dynamics. 




code: bcd_hard_spheroids_widoms.c
input file: runner_bhsw.sh

Description: This code simulate the system of hard Brownian spheroids and calculate the dynamics, structure as well as widom chemical potential of the system, following the formulation of Brownian cluster dynamics. 




code: bcd_hard_spheroids_obstacle.c
input file: runner_bhso.sh

Description: This code simulate the system of hard Brownian spheroids and calculate the dynamics in the presence of immobile obstacles arranged in a periodic manner. The code follows the formulation of Brownian cluster dynamics. 




code: bcd_hard_spheroids_obstacle_widoms.c
input file: runner_bhsow.sh

Description: This code simulate the system of hard Brownian spheroids and calculate the dynamics and structure in the presence of immobile obstacles arranged in a periodic manner. It also calculates the Widom's chemical potential. The code follows the formulation of Brownian cluster dynamics. 




code: bcd_hard_patchy_spheroids.c
input file: runner_bhps.sh

Description: This code simulate the system of hard brownian spheroids decorated with patches and calculate the dynamics, structure following the formulation of brownian cluster dynamics. The purpose of this code is to study the structure and the kinetic of the self assembly in simple colloidal system.




code: nvt_hard_patchy_spheroids.c
input file: runner_nvhps.sh

Description: This code simulate the system of hard brownian spheroids decorated with patches and calculate the dynamics, structure following the formulation of Monte Carlo simulation in NVT ensemble. The purpose of this code is to study the structure, phase diagram and to equillibriate faster.




code: npt_hard_patchy_spheroids.c
input file: runner_nphps.sh

Description: This code simulate the system of hard brownian spheroids decorated with patches and calculate the dynamics, structure following the formulation of Monte Carlo simulation in NpT ensemble. The purpose of this code is to study the structure, phase diagram at a constant pressure.




code: hybrid_monte_carlo_integral_solver.c 
input file: runner_integral_solver.sh

Description" This code integrate the classical density functional free energy expression by calculating the mayer f function around a cylinder as discussed in the published work in the given manual.





code: gibbs_ensemble_spherical.c
input file: runner_gibbs_ensemble_spherical.sh

Description: This code simulate the system over a spherical surface following the Gibbs ensemble simulation framework. We use it predict the structure of the closed system like spherical shells, as referred in the manual.





code: obstacle_analyser.c
input file: runner_analyser.sh

Description: This code is an example for the analysis of output data obtained from the code for the simulation of spheroids in the presence of cylindrical obstacles.

