#!/bin/bash


#simulation parameters


echo "0"  >  Base_data.txt  #0: to start from beginning or 1: to start from a relaxed_system.txt file



echo "0"  >>  Base_data.txt  #0: cubic boxes  1: floppy box   



echo "1000" >> Base_data.txt  #number of particles



echo "0.01" >> Base_data.txt  #step length



echo "500000000" >> Base_data.txt  #simulation steps



echo "0.1" >> Base_data.txt  #pressure: it should be higher than 0.004 otherwise you would need to go and change the fvf value in the main script



echo "1" >> Base_data.txt  #kind of particles in the system e.g. 2 







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

echo "1.0" >> Base_data.txt 

echo "50.0" >> Base_data.txt







#patch properties




echo "1" >> Base_data.txt  #kind of patches in the system e.g. 2 



echo "50.0" >> Base_data.txt  #0: switch off the interaction between patches or 1: to switch on the interaction
	
echo "1" >> Base_data.txt  # between all patches with type id 0 and 0

echo "50.0" >> Base_data.txt







#patchy potential



echo "50.0" >> Base_data.txt  #epsi

echo "0.1" >> Base_data.txt  # between all patches with type id 0 and 0

echo "50.0" >> Base_data.txt 



echo "50.0" >> Base_data.txt  #beta

echo "0.00002">> Base_data.txt  # between all patches with type id 0 and 0

echo "50.0" >> Base_data.txt



echo "50.0" >> Base_data.txt  #omega

echo "0.95" >> Base_data.txt  # for patches with type id 0 

echo "50.0" >> Base_data.txt 



echo "50.0" >> Base_data.txt  #del_omega

echo "0.05" >> Base_data.txt  # for patches with type id 0 

echo "50.0" >> Base_data.txt 






#patch decoration over the particles



echo "50.0" >> Base_data.txt  #number of patches on each particle types

echo "2" >> Base_data.txt  #2 patches on type 0 particle

echo "50.0" >> Base_data.txt



echo "50.0" >> Base_data.txt  #patches on each particle type

echo "0 0" >> Base_data.txt  #decorated patches on each particles type where two patches with patch type id 0 and 0 is decorated over 0 type of particle # the patch vector to each patch type decorated on the particles must be supplied from the Patch_vector.txt file, in the same sequence, where patch vector is written in the file row by row

echo "50.0" >> Base_data.txt 




gcc npt_hard_patchy_spheroids.c -o npt_hard_patchy_spheroids -lm
./npt_hard_patchy_spheroids
