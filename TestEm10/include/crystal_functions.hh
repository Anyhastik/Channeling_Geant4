#pragma once

#ifndef CRYSTAL_FUNCTIONS_HPP
#define CRYSTAL_FUNCTIONS_HPP

// #ifndef CRYSTAL_FUNCTIONS_H
// #define CRYSTAL_FUNCTIONS_H 1

#include "phys_const.hh"
// #include "crystal_functions.cc"

// class CRYSTAL_FUNCTIONS
// {
	// public:
    // CRYSTAL_FUNCTIONS();
    // ~CRYSTAL_FUNCTIONS();

ld1 ducr(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

ld1 ducr_2(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

ld1 ucr(ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

ld1 eve(ld1 i, array<ld1, 2*nmax+1> cpot, ld1 dp, ld1 tetta1, string cr_t);

ld1 functz(ld1 t1z, ld1 z1, ld1 z2, ld1 x2, ld1 k0);

ld1 functiti(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp);

void print(vector<ld1> &vec);

ld1 dWdE(ld1 Nt, ld1 w, ld1 TT, ld1 integral[nbeam]);

ld1 betaz(ifstream &infile);

ld1 integr_f(ld1 x, ld1 eve2, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

ld1 simpson(ld1 a, ld1 b, int n0, ld1 eve2, array<ld1, 2*nmax+1> cpot,ld1 dp, string cr_t);

ld1 find(ld1 x, ld1 eps, ld1 eve1, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

ld1 theta(ld1 x);

ld1 omega(ld1 n, ld1 TT);

ld1 trapezoid(vector<ld1> times, vector<ld1> velocities, int n, ld1 TT);
// };
#endif /*CRYSTAL_FUNCTIONS_HPP*/