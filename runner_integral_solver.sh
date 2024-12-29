#!/bin/bash

echo "5" > Base_data.txt #number of particles


echo "0.01" >> Base_data.txt #step length



echo "200000" >> Base_data.txt #simulation steps



echo "0.1" >> Base_data.txt #volume fraction



echo "0.95" >> Base_data.txt #theta angle



echo "2" >> Base_data.txt   #kind of particles in the system e.g. 2 


echo "50.0" >> Base_data.txt 

echo "0" >> Base_data.txt   #sco     0: constant volume  1: constant axis
echo "0" >> Base_data.txt
echo "0" >> Base_data.txt

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt 

echo "2.0" >> Base_data.txt   #aspect ratio defined for each kind of particles default is 1.0
echo "1.00" >> Base_data.txt
echo "1.5" >> Base_data.txt

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt 

echo "2.00" >> Base_data.txt   #size enhancement parameter 1.0 in ordinary case  !!!!!!!
echo "1.0" >> Base_data.txt
echo "1.0" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt

echo "0.6" >> Base_data.txt   #fraction of each kind of particles.
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt






gcc hybrid_monte_carlo_integral_solver.c -o hybrid_monte_carlo_integral_solver -lm
./hybrid_monte_carlo_integral_solver
		
