#pragma once

// UnityDLLExample.h : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>

//model structure
typedef struct {
	double * x; //computational scratch space
	double * x2; //computational scratch space
	double * store; //computational scratch space
	double * storep;
	double * k1; //computational scratch space
	double * k2;  //computational scratch space
	double * k3;  //computational scratch space
	double * k4;  //computational scratch space
	double * XPRIME;		// Computational scratch space
	double * XPRIME1;		// Computational scratch space
	double * XPRIME2;		// Computational scratch space
	double * XPRIME3;		// Computational scratch space
	double * XPRIME4;		// Computational scratch space
	double * mass;  //masses
	double * srad2;
	 
	/*double softening_factor;*/
	double G;  //graviational constant
	double m;  
	int n;  //number of points
	bool LeapFrogFirstIteration;  //Used for LeapFrog method half step on first iteration
	int abmCounter; //restart counter for ABM integration 

	double shield; //shield radius
	double srad_factor;

	double softening_factor; //softening potential factor
	

} Model;

extern "C" __declspec(dllexport) void * createModel(int n);
extern "C" __declspec(dllexport) void stepModelRK2(void * voo, double dt);
extern "C" __declspec(dllexport) void stepModelRK4(void * voo, double dt);
//Added LeapFrog Method attempt 
extern "C" __declspec(dllexport) void stepNbodyModelLeapfrog(void * voo, double dt);

//Added ABM attempt
extern "C" __declspec(dllexport) void stepNbodyModelABM(void * voo, double tStep);

extern "C" __declspec(dllexport) void stepModelEuler(void * voo, double dt);
extern "C" __declspec(dllexport) double getX(void * voo, int i);
extern "C" __declspec(dllexport) double * getXArray(void * voo);
extern "C" __declspec(dllexport) void destroyModel(void * voo);
extern "C" __declspec(dllexport) void ratesOfChange(Model * foo, double* x, double* xdot);
extern "C" __declspec(dllexport) double randRange(double low, double high);
extern "C" __declspec(dllexport) void setX(void * voo, double value, int i);
extern "C" __declspec(dllexport) void setXArray(void * voo, double * value);
extern "C" __declspec(dllexport) void setTStep(void * voo, double tstep);

//softening potential
extern "C" __declspec(dllexport) void setSofteningNBodyModel(void * voo, double softening_factor);
//shield radius
extern "C" __declspec(dllexport) void setSradNbodyModel(void * voo, double srad_factor);
extern "C" __declspec(dllexport) double computeSoftenedRadius(double g_m, double tstep_squared, double srad_factor);
extern "C" __declspec(dllexport) double fcr_guess(double x);
extern "C" __declspec(dllexport) double fcr(double x, int MAX_ITER);

extern "C" __declspec(dllexport) void setG(void * voo, double G);
