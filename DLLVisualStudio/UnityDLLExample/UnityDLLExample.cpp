// UnityDLLExample.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include "UnityDLLExample.h"

extern "C" __declspec(dllexport)  double randRange(double low, double high) {
	double r = (double)rand() / (double)RAND_MAX;
	return low + r*(high - low);
}

//softening potentials
extern "C" __declspec(dllexport) void setSofteningNBodyModel(void * voo, double softening_factor) {
	Model * foo = (Model *)voo;
	foo->softening_factor = softening_factor;
}

//shield radius methods
extern "C" __declspec(dllexport) void setSradNbodyModel(void * voo, double srad_factor)
{
	Model *foo = (Model *)voo;
	foo->srad_factor = srad_factor;
}

extern "C" __declspec(dllexport) double computeSoftenedRadius(double g_m, double tstep_squared, double srad_factor) {
	 //g_m = G*m;
	if (srad_factor>0.0) {
		return srad_factor * fcr(g_m*tstep_squared, 0);
	}
	else {
		return 0.0;
	}
}


/*Shield radius methods... might make into seperate class*/
extern "C" _declspec(dllexport) double fcr(double x, int MAX_ITER) {
	double a, b;
	int done = false;
	int count = 0;
	int sign = 0;

	if (x == 0.0) return 0.0;
	if (x<0.0) { x = -x; sign = 1; }
	if (x == 1.0) { done = true; }
	a = fcr_guess(x);
	b = a;

	while (!done&&count++<MAX_ITER) {
		b = a * (((a*a*a + x) + x) / (a*a*a + (a*a*a + x)));
		if (fabs((b - a) / a)<1.0e-6) done = true;
		else if (fabs(b - a)<1.0e-50) done = true;
		a = b;
	}
	if (sign) return -b;
	else return b;
}

extern "C" _declspec(dllexport) double fcr_guess(double x) {
	if (x > 1.0) {
		if (x < 1.0e1) return 2.0e0;
		else if (x < 1.0e2) return 5.0e0;
		else if (x < 1.0e3) return 1.0e1;
		else if (x < 1.0e4) return 2.0e1;
		else if (x < 1.0e5) return 5.0e1;
		else if (x < 1.0e6) return 1.0e2;
		else if (x < 1.0e7) return 2.0e2;
		else if (x < 1.0e8) return 5.0e2;
		else if (x < 1.0e9) return 1.0e3;
		else if (x < 1.0e10) return 2.0e3;
		else if (x < 1.0e11) return 5.0e3;
		else if (x < 1.0e12) return 1.0e4;
		else if (x < 1.0e15) return 1.0e5;
		else if (x < 1.0e18) return 1.0e6;
		else if (x < 1.0e21) return 1.0e7;
		else if (x < 1.0e24) return 1.0e8;
		else if (x < 1.0e27) return 1.0e9;
		else return 1.0e10;
	}
	else if (x < 1.0) {
		if (x > 1.0e-1) return 5.0e-1;
		else if (x > 1.0e-2) return 2.0e-1;
		else if (x > 1.0e-3) return 1.0e-1;
		else if (x > 1.0e-4) return 5.0e-2;
		else if (x > 1.0e-5) return 2.0e-2;
		else if (x > 1.0e-6) return 1.0e-2;
		else if (x > 1.0e-7) return 5.0e-3;
		else if (x > 1.0e-8) return 2.0e-3;
		else if (x > 1.0e-9) return 1.0e-3;
		else if (x > 1.0e-10) return 5.0e-4;
		else if (x > 1.0e-11) return 2.0e-4;
		else if (x > 1.0e-12) return 1.0e-4;
		else if (x > 1.0e-15) return 1.0e-5;
		else if (x > 1.0e-18) return 1.0e-6;
		else if (x > 1.0e-21) return 1.0e-7;
		else if (x > 1.0e-24) return 1.0e-8;
		else if (x > 1.0e-27) return 1.0e-9;
		else return 1.0e-10;
	}
	else if (x == 1.0) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}

extern "C" __declspec(dllexport) void setMass(void * voo, double m) {
	Model * foo = (Model *)voo;
	foo->m = m;
}

extern "C" __declspec(dllexport) void setMasses(void * voo, double * m) {
	Model * foo = (Model *)voo;
	int n = foo->n;
	foo->m = 0.0;
	for (int i = 0; i < n; i++) {
		foo->mass[i] = m[i];
		foo->m += m[i];
	}
}

extern "C" __declspec(dllexport) void setG(void * voo, double G) {
	Model * foo = (Model *)voo;
	foo->G = G;
}

/*extern "C" _declspec(dllexport) void setTStep(void * voo, double tstep) {
	Model * foo = (Model *)voo;
	foo->t = tstep;
}*/

extern "C" __declspec(dllexport) void * allocModel(int n) {
	Model *foo = (Model *)malloc(sizeof(Model));
	foo->n = n;
	foo->x = (double *)malloc(sizeof(double) * 6 * n);
	foo->x2 = (double *)malloc(sizeof(double) * 6 * n);
	foo->store = (double *)malloc(sizeof(double) * 6 * n);
	foo->storep = (double *)malloc(sizeof(double) * 6 * n);
	foo->k1 = (double *)malloc(sizeof(double) * 6 * n);
	foo->k2 = (double *)malloc(sizeof(double) * 6 * n);
	foo->k3 = (double *)malloc(sizeof(double) * 6 * n);
	foo->k4 = (double *)malloc(sizeof(double) * 6 * n);
	foo->XPRIME = (double *)malloc(sizeof(double)*n * 6);
	foo->XPRIME1 = (double *)malloc(sizeof(double)*n * 6);
	foo->XPRIME2 = (double *)malloc(sizeof(double)*n * 6);
	foo->XPRIME3 = (double *)malloc(sizeof(double)*n * 6);
	foo->XPRIME4 = (double *)malloc(sizeof(double)*n * 6);
	foo->srad2 = (double *)malloc(sizeof(double)*n * 6); //shield radius array
	foo->mass = (double *)malloc(sizeof(double)*n);
	foo->abmCounter = -3; 
	foo->LeapFrogFirstIteration = true;
	foo->shield = 1.0e-4;
	foo->softening_factor = 1.0e-4;
	foo->G = 1.0;
	foo->m = 1.0;

	return foo;
}

extern "C" __declspec(dllexport) void initModel(void * voo) {
	Model *foo = (Model *)voo;
	int n = foo->n;


	for (int i = 0; i < n; i++) {
		foo->x[i * 6 + 0] = randRange(-1, 1);
		foo->x[i * 6 + 1] = randRange(-1, 1);
		foo->x[i * 6 + 2] = randRange(-1, 1);
		foo->x[i * 6 + 3] = 0.0;
		foo->x[i * 6 + 4] = 0.0;
		foo->x[i * 6 + 5] = 0.0;
		foo->mass[i] = foo->m / (double)n;
	}

}


extern "C" __declspec(dllexport) void stepModelEuler(void * voo, double h) {
	Model * foo = (Model *)voo;

	int nEquations = foo->n*6;

	ratesOfChange(foo, foo->x, foo->k1);
	for (int i = 0; i < nEquations; i++) {
		foo->x[i] = foo->x[i] + foo->k1[i] * h;
	}
}


extern "C" __declspec(dllexport) void stepModelRK4(void * voo, double h) {
	Model * foo = (Model *)voo;

	int nEquations = foo->n*6;

	ratesOfChange(foo,foo->x, foo->k1);
	for (int i = 0; i < nEquations; i++) {
		foo->store[i] = foo->x[i] + foo->k1[i] * h / 2.0;
	}
	ratesOfChange(foo,foo->store, foo->k2);
	for (int i = 0; i < nEquations; i++) {
		foo->store[i] = foo->x[i] + foo->k2[i] * h / 2.0;
	}
	ratesOfChange(foo,foo->store, foo->k3);
	for (int i = 0; i < nEquations; i++) {
		foo->store[i] = foo->x[i] + foo->k3[i] * h;
	}
	ratesOfChange(foo,foo->store, foo->k4);
	for (int i = 0; i < nEquations; i++) {
		foo->x[i] = foo->x[i] + (foo->k1[i] + 2.0*foo->k2[i] + 2.0*foo->k3[i] + foo->k4[i]) * h / 6.0;
		//storep is a more accurate guess of actual value 
		//use in predictor corrector for best guess
		foo->storep[i] = (foo->k1[i] + 2.0*foo->k2[i] + 2.0*foo->k3[i] + foo->k4[i]) / 6.0;
	}
}

//Leap frog method 
extern "C" __declspec(dllexport) void stepNbodyModelLeapfrog(void * voo, double tStep) 
{
	Model * foo = (Model *)voo;

	ratesOfChange(foo, foo->x, foo->k1);
	if (foo->LeapFrogFirstIteration) {
		// setup leapfrog on first step, change velocities by a half step
		for (int i = 0; i < foo->n; i++) {
			foo->x[i * 6 + 3] += 0.5*foo->k1[i * 6 + 3] * tStep;
			foo->x[i * 6 + 4] += 0.5*foo->k1[i * 6 + 4] * tStep;
			foo->x[i * 6 + 5] += 0.5*foo->k1[i * 6 + 5] * tStep;
		}
	}
	else {
		// update v,x
		for (int i = 0; i < foo->n * 6; i++) {
			foo->x[i] += foo->k1[i] * tStep;
		}
	}
	foo->LeapFrogFirstIteration=false;
}

/*Predictor Corrector*/
extern "C" __declspec(dllexport) void stepNbodyModelABM(void * voo, double tStep) {
	Model *foo = (Model*)voo;

	double * fk3 = NULL;
	double * fk2 = NULL;
	double * fk1 = NULL;
	double * fk0 = NULL;
	double * fkp = NULL;

	// determine if previous steps exist, if not, populate w/ RK4
	if (foo->abmCounter < 0) {  //counter starts at -3 because ABM needs certain information before it can run
		stepModelRK4(foo, tStep);  //RK4 step first 
		if (foo->abmCounter == -3) {
			for (int i = 0; i < foo->n * 6; i++)
				foo->XPRIME4[i] = foo->storep[i];
		}
		else if (foo->abmCounter == -2) {
			for (int i = 0; i < foo->n * 6; i++)
				foo->XPRIME3[i] = foo->storep[i];
		}
		else {
			for (int i = 0; i < foo->n * 6; i++)
				foo->XPRIME2[i] = foo->storep[i];
		}
	}
	else {
		if (foo->abmCounter % 5 == 0) {
			fk3 = foo->XPRIME4;
			fk2 = foo->XPRIME3;
			fk1 = foo->XPRIME2;
			fk0 = foo->XPRIME1;
			fkp = foo->XPRIME;
		}
		else if (foo->abmCounter % 5 == 1) {
			fk3 = foo->XPRIME3;
			fk2 = foo->XPRIME2;
			fk1 = foo->XPRIME1;
			fk0 = foo->XPRIME;
			fkp = foo->XPRIME4;
		}
		else if (foo->abmCounter % 5 == 2) {
			fk3 = foo->XPRIME2;
			fk2 = foo->XPRIME1;
			fk1 = foo->XPRIME;
			fk0 = foo->XPRIME4;
			fkp = foo->XPRIME3;
		}
		else if (foo->abmCounter % 5 == 3) {
			fk3 = foo->XPRIME1;
			fk2 = foo->XPRIME;
			fk1 = foo->XPRIME4;
			fk0 = foo->XPRIME3;
			fkp = foo->XPRIME2;
		}
		else if (foo->abmCounter % 5 == 4) {
			fk3 = foo->XPRIME;
			fk2 = foo->XPRIME4;
			fk1 = foo->XPRIME3;
			fk0 = foo->XPRIME2;
			fkp = foo->XPRIME1;
		}
		ratesOfChange(foo, foo->x,fk0);
		for (int i = 0; i < foo->n * 6; i++) { //predictor step
			foo->x2[i] = foo->x[i] +
				(tStep / 24.0)*(-9.0*fk3[i] + 37.0*fk2[i]
					- 59.0*fk1[i] + 55.0*fk0[i]);
		}
		ratesOfChange(foo, foo->x2, fkp);
		for (int i = 0; i < foo->n * 6; i++) { //correction step
			foo->x[i] = foo->x[i] +
				(tStep / 24.0)*(fk2[i] - 5.0*fk1[i] +
					19.0*fk0[i] + 9.0*fkp[i]);
		}
	}
	foo->abmCounter++;
}

extern "C" __declspec(dllexport) void stepModelRK2(void * voo, double dt) {
	Model * foo = (Model *)voo;
	for (int i = 0; i < 6*foo->n; i++) {
		foo->store[i] = foo->x[i];
	}
	ratesOfChange(foo, foo->x, foo->k1);
	for (int i = 0; i < 6 * foo->n; i++) {
		foo->x[i] += foo->k1[i] * dt;
	}
	ratesOfChange(foo, foo->x, foo->k2);
	for (int i = 0; i < 6 * foo->n; i++) {
		foo->x[i] = foo->store[i]+0.5*(foo->k1[i]+foo->k2[i]) * dt;
	}
}

extern "C" __declspec(dllexport) void setX(void * voo,double value, int i) {
	Model * foo = (Model *)voo;
	foo->x[i] = value;
}

extern "C" __declspec(dllexport) void setXArray(void * voo, double * value) {
	Model * foo = (Model *)voo;
	int n = foo->n;
	for (int i = 0; i < 6*n; i++) {
		foo->x[i] = value[i];
	}
}

extern "C" __declspec(dllexport) double getX(void * voo, int i) {
	Model * foo = (Model *)voo;
	return foo->x[i];
}

extern "C" __declspec(dllexport) double * getXArray(void * voo) {
	Model * foo = (Model *)voo;
	return foo->x;
}


extern "C" __declspec(dllexport) void destroyModel(void * voo) {
	Model * foo = (Model *)voo;
	free(foo->x);
	free(foo->x2);
	free(foo->mass);
	free(foo->store);
	free(foo->storep);
	free(foo->XPRIME);
	free(foo->XPRIME1);
	free(foo->XPRIME2);
	free(foo->XPRIME3);
	free(foo->XPRIME4);
	free(foo->k1);
	free(foo->k2);
	free(foo->k3);
	free(foo->k4);
	free(foo->srad2);
	// add all the other frees
	free(foo);
}

extern "C" __declspec(dllexport) void ratesOfChange(Model * foo, double* x, double* xdot) {
	int n = foo->n;
	double srad, tstep_squared;
	for (int i = 0; i < n; i++) {
		xdot[i * 6 + 0] = x[i * 6 + 3];
		xdot[i * 6 + 1] = x[i * 6 + 4];
		xdot[i * 6 + 2] = x[i * 6 + 5];
		xdot[i * 6 + 3] = 0.0;
		xdot[i * 6 + 4] = 0.0;
		xdot[i * 6 + 5] = 0.0;
		srad = computeSoftenedRadius(foo->G*foo->mass[i], tstep_squared, foo->srad_factor);
		foo->srad2[i] = srad * srad;
	}
	for (int i = 0; i < n; i++) {
		double xi = x[i * 6 + 0];
		double yi = x[i * 6 + 1];
		double zi = x[i * 6 + 2];
		double axi = 0.0;
		double ayi = 0.0;
		double azi = 0.0;
		{
			for (int j = i + 1; j < n; j++) {
				double xj = x[j * 6 + 0];
				double yj = x[j * 6 + 1];
				double zj = x[j * 6 + 2];
				double dx = xi - xj;
				double dy = yi - yj;
				double dz = zi - zj;

				double dr = sqrt(dx*dx + dy*dy + dz*dz + foo->softening_factor*foo->softening_factor);
				double dr2 = dr * dr;
				double dr3 = dr2*dr;
				double accel = -foo->G / dr3;
				if (dr2 > foo->srad2[j])
				{
					axi += accel * dx * foo->mass[j];
					ayi += accel * dy * foo->mass[j];
					azi += accel * dz * foo->mass[j];
				}
				if (dr2 > foo->srad2[i])
				{
					xdot[j * 6 + 3] -= accel * dx * foo->mass[i];
					xdot[j * 6 + 4] -= accel * dy * foo->mass[i];
					xdot[j * 6 + 5] -= accel * dz * foo->mass[i];
				}
			}
				/*
				double accel = -foo->G / dr2;
				axi += accel*dx / dr*foo->mass[j];
				ayi += accel*dy / dr*foo->mass[j];
				azi += accel*dz / dr*foo->mass[j];
				xdot[j * 6 + 3] -= accel*dx / dr*foo->mass[i];
				xdot[j * 6 + 4] -= accel*dy / dr*foo->mass[i];
				xdot[j * 6 + 5] -= accel*dz / dr*foo->mass[i];
			}
			*/
			{
				xdot[i * 6 + 3] += axi;
				xdot[i * 6 + 4] += ayi;
				xdot[i * 6 + 5] += azi;
			}
		}
	}
}

