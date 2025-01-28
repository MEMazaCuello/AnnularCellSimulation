# AnnularCellSimulation
C++ simulation of a 2D system of rods inside a vibrating annular cell via a Monte Carlo algorithm.

Rods are represented as rectangles defined by a characteristic WIDTH and LENGTH and their positions are given in cartesian coordinates, plus the angle between their long axis and the OX axis (in radians).

The simulation parameters defining the lengths of the Rods and the AnnularCell are set in 'GlobalParameters.hpp'. 
The 'main.cpp' file includes some use cases.
Compile using C++20 standard. 

A previous version of this code was used for the simulations of

> A. Díaz-De Armas _et al._, _Domain walls in vertically vibrated monolayers of cylinders confined in annuli_, Phys. Rev. Research **2**, 033436 (2020)
