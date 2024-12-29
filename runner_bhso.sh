#!/bin/bash



#simulation parameters



echo "0"  >  Base_data.txt  #0: to start from beginning or 1: to start from a relaxed_system.txt file


echo "1" >> Base_data.txt  #for finite volume fraction :1    single particle :0


echo "4.5" >> Base_data.txt #radius of the obstacles


echo "0.5" >> Base_data.txt #area fraction of the obstacles


echo "0.005" >> Base_data.txt  #step length


echo "500000000" >> Base_data.txt  #simulation steps


echo "20" >> Base_data.txt  #probation time is the time the system is left to randomize after that the simulation will start


echo "0.1" >> Base_data.txt  #volume fraction


echo "1" >> Base_data.txt  #kind of particles in the system e.g. 1 or 2 







#particle's shape and size



echo "50.0" >> Base_data.txt  #sco 0: constant volume  1: constant axis

echo "0" >> Base_data.txt  

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt  #aspect ratio define for each kind of particles default is 1.0

echo "1.0" >> Base_data.txt 

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt  #size enhancement parameter default is 1 

echo "1.0" >> Base_data.txt   

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt  #fraction of each kind of particles

echo "1" >> Base_data.txt 

echo "50.0" >> Base_data.txt






gcc bcd_hard_spheroids_obstacle.c -o bcd_hard_spheroids_obstacle -lm
./bcd_hard_spheroids_obstacle
