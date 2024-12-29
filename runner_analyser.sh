#!/bin/bash

echo "0" > Base_data.txt #To start the Simulation from the zero: 0, to start from a relaxed file: 1



echo "1" >> Base_data.txt #For single particle diffsusion 0: for finite volume fraction: 1



echo "9.657" >> Base_data.txt #radius of the obstacles



echo "0.785" >> Base_data.txt #area fraction of the obstacles



echo "0.005" >> Base_data.txt #step length



echo "100000000" >> Base_data.txt #simulation steps



echo "20" >> Base_data.txt #probation periods



echo "0.2" >> Base_data.txt #volume fraction



echo "1" >> Base_data.txt #Type of particles in the system









echo "50.0" >> Base_data.txt 

echo "0" >> Base_data.txt #sco     0: constant volume  1: constant axis

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt 

echo "3.000" >> Base_data.txt #aspect ratio define for each kind of particles  default is 1.0

echo "50.0" >> Base_data.txt




echo "50.0" >> Base_data.txt 

echo "1.0" >> Base_data.txt   #size enhancement parameter 1.0 in ordinary case  !!!!!!!

echo "50.0" >> Base_data.txt



echo "50.0" >> Base_data.txt 

echo "1.0" >> Base_data.txt   #fraction of the particle enhancement parameter 1.0 in ordinary case  !!!!!!!

echo "50.0" >> Base_data.txt







echo "1" >> Base_data.txt   # t-phy-1
echo "500" >> Base_data.txt  # t-phy-2
echo "1" >> Base_data.txt  # Put the number of time over which the system will time average the data
echo "40" >> Base_data.txt  # Time interval between the data for the time average




gcc obstacle_analyser.c -o obstacle_analyser -lm
./obstacle_analyser
