#pragma once
#ifndef MDMAKE_OPT
#define MDMAKE_OPT

#include <math.h>
using namespace std;

struct Coordinate{
        double x, y, z;
};

struct optpara{
	string atomcomp;
	string atomratio;
	int seed;
	string stress;
	string strain;
	string potential;
	string cellnum;
	string pbc_x;
	string pbc_y;
	string pbc_z;
	string shape;
        string method; //takizawa
	string periodic; //takizawa
	int omp;
      string stack;
};

#endif

