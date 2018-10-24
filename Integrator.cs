using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Integrator is an abstract class for integrating a system of ODEs
/// </summary>
abstract public class Integrator  {

	int nEquations;
	double [] store;
	double [] k1;
	double [] k2;
	double [] k3;
	double [] k4;

	public double [] getK3() {
		return k3;
	}

	/// <summary>
	/// Allocate memory for all storage arrays and set number of equations
	/// </summary>
	/// <param name="nEquations">N equations.</param>
	public void Init (int nEquations) {
		// set up temp arrays
		this.nEquations = nEquations;
		store = new double[nEquations];
		k1 = new double[nEquations];
		k2 = new double[nEquations];
		k3 = new double[nEquations];
		k4 = new double[nEquations];
	}

	/// <summary>
	/// Abstract void, override this method to set the ODEs to be
	/// integrated.
	/// </summary>
	/// <param name="x">The values being integrated.</param>
	/// <param name="xdot">The derivatives being calculated.</param>
	abstract public void RatesOfChange (double[] x, double[] xdot);
		
	/// <summary>
	/// Step forward using Euler's method
	/// </summary>
	/// <param name="x">The values being integrated.</param>
	/// <param name="h">The time step.</param>
	public void EulerStep(double [] x, double h) {
		RatesOfChange(x,k1);
		for(int i=0;i<nEquations;i++) {
			x[i] += k1[i]*h;
		}
	}

	/// <summary>
	/// Step forward using 4th order Runge Kutta method
	/// </summary>
	/// <param name="x">The values being integrated.</param>
	/// <param name="h">The time step.</param>
	public void RK4Step(double [] x, double h) {
		RatesOfChange (x,k1);
		for (int i = 0; i < nEquations; i++) {
			store [i] = x [i] + k1 [i] * h / 2.0;
		}
		RatesOfChange (store,k2);
		for (int i = 0; i < nEquations; i++) {
			store [i] = x [i] + k2 [i] * h / 2.0;
		}
		RatesOfChange (store,k3);
		for (int i = 0; i < nEquations; i++) {
			store [i] = x [i] + k3 [i] * h;
		}
		RatesOfChange (store,k4);
		for (int i = 0; i < nEquations; i++) {
			x [i] = x [i] + (k1[i] +2.0*k2[i]+ 2.0*k3 [i]+k4[i]) * h/6.0;
		}
	}
}
