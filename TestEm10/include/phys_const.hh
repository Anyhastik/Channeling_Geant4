// #pragma once
#ifndef PHYS_CONST_h
#define  PHYS_CONST_h 1


#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#define ld1 long double
#define d double
#include "globals.hh"
using namespace std;
#define c0 299792458
#define echarge 4.8032 * pow(10, -10) * sqrt(1 / (1.60217733 * pow(10, -6)) * pow(10, -2))
#define nmax 20
#define hbar 6.582 * pow(10, -22)        //(*MeV s*)
#define m0 0.5109996 * pow(10, 6)        //(*ev/c^2*)
#define e0 0.5109996                    //(*MeV*)
#define e 250                          // MeV INPUT!!!!!
#define mom sqrt(pow(e, 2) - pow(e0, 2)) //(*MeV*)
// #define pi 3.14159265358979323846      // GEANT
#define gama  e / e0
#define coeff  1.123443973421022 * pow(10, 28) // coeff = - pow(c0, 2) * pow(10, 20) / (gama * m0);
// #define dp 1.35775 //3.13559
#define tetta 0
#define nbeam 50
#define n_steps 3500
  
class GPHYS_CONST
{
	public:
    GPHYS_CONST();
    ~GPHYS_CONST();
	
static tuple<ld1, array<ld1, 2*nmax+1>> cpotss(G4int sigen, G4int k0, G4int l0, G4int n0, G4String crystaltype);
};
// #include "phys_const.cc"
#endif  /* PHYS_CONST*/