#!/bin/bash

echo "100" > Base_data.txt #number of particles





echo "0.01" >> Base_data.txt #step length





echo "500000000" >> Base_data.txt #simulation steps





echo "0.15" >> Base_data.txt # Initial area fraction of both the spheres





echo "1" >> Base_data.txt #kind of particles in the system e.g. 2 







echo "50.0" >> Base_data.txt 

echo "0" >> Base_data.txt   #sco =   0: constant volume  1: constant axis
echo "0" >> Base_data.txt
echo "0" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt 

echo "1.0" >> Base_data.txt   #aspect ratio defined for each kind of particles  default is 1.0  !!!!!!!
echo "1.0" >> Base_data.txt
echo "1.5" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt 

echo "1.0" >> Base_data.txt   #size enhancement parameter 1.0 in ordinary case  !!!!!!!
echo "1.0" >> Base_data.txt
echo "1.0" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt

echo "0.5" >> Base_data.txt #fraction of each kind of particles.
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt
	
echo "2" >> Base_data.txt  #kind of interaction	0: hard core; 1: HC+Jenus; 2: HC+Patchy; 3: HC+Isotropic; 4: HC+jenus+Isotropic; 5: HC+patchy+Isotropic; 
echo "3" >> Base_data.txt
echo "3" >> Base_data.txt
echo "0" >> Base_data.txt
echo "0" >> Base_data.txt
echo "2" >> Base_data.txt
echo "2" >> Base_data.txt
echo "2" >> Base_data.txt
echo "2" >> Base_data.txt

echo "50.0" >> Base_data.txt





#Jenus potential
echo "50.0" >> Base_data.txt #epsi n+(n/2)(n-1)

echo "0.0" >> Base_data.txt 
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt





echo "50.0" >> Base_data.txt #beta

echo "0.0" >> Base_data.txt
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt 





#Patchy potential
echo "50.0" >> Base_data.txt #epsi

echo "0.25" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt
echo "0.1" >> Base_data.txt

echo "50.0" >> Base_data.txt 




echo "50.0" >> Base_data.txt #beta

echo "0.0008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt
echo "0.00008" >> Base_data.txt

echo "0.0" >> Base_data.txt


echo "50.0" >> Base_data.txt





echo "50.0" >> Base_data.txt #omega

echo "0.22" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt

echo "50.0" >> Base_data.txt 




echo "50.0" >> Base_data.txt #wom

echo "0.075" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt
echo "0.9" >> Base_data.txt

echo "50.0" >> Base_data.txt 






#Isotropic potential
echo "50.0" >> Base_data.txt #epsi

echo "0.1" >> Base_data.txt 
echo "0.1" >> Base_data.txt 
echo "0.1" >> Base_data.txt 
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt






echo "50.0" >> Base_data.txt #beta

echo "0.0106" >> Base_data.txt
echo "0.0106" >> Base_data.txt 
echo "0.0106" >> Base_data.txt 
echo "0.0" >> Base_data.txt

echo "50.0" >> Base_data.txt


gcc gibbs_ensemble_spherical.c -o gibbs_ensemble_spherical -lm
./gibbs_ensemble_spherical

		
