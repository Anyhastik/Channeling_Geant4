# Channeling_Geant4
This is an implementation of a new C++ module for channeling radiation (CR) simulation in Geant4 [1]. The module is built into Geant4 as a discrete set of physical processes [2]. Such a form of integration makes it possible to combine new physical processes in single crystals with those already built into Geant4 (ordinary bremsstrahlung, different types of scattering of both photons and primary particles, etc.), which significantly increases the quality of simulation. The direct simulation of CR in Geant4 allows describing devices such as, e. g. positron sources, without using external programs for CR calculation.
[1] S. Agostinelli, J. Allison, K. Amako et al., Nucl. Instrum. Meth. A 506, 250 (2003). 
[2] A.A. Savchenko and W. Wagner, JINST 16 P12042 (2021).
