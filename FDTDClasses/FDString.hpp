//
//  StringClass.hpp
//  FDTD_C_String
//
//  Created by admin on 07/11/2017.
//  Copyright © 2017 admin. All rights reserved.
//

#ifndef StringClass_hpp
#define StringClass_hpp

#include <iostream>			// std::cout
#include <cmath>

class FDString {
	
public:
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constructors/Assignments
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//Still Not sure About the Correct Syntax for these
	
	//Contructor
	FDString();
	
	//Destructor
	~FDString(); //{delete EVERYTHING;}
	
	// Copy Assignment
	//	FDString& operator= (const FDString&);
	
	// Move-Constructor
	//	FDString (FDString&&);
	
	// Move-Assignment
	//	FDString& operator= (FDString&&);
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Methods
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void Setup(double, int, bool);	//Semi-Finished
	void setLoss(int);				//DONE
	void setGrid();
	void setCoefs(bool);			//DONE
	void UpdateScheme();			//DONE
	void ProcessIO();
	void printInfo();				//DONE
	void printCoefs();				//DONE
	void setInput(float);			//DONE
	void setOutput(float);			//DONE
	double getOutput(bool);			//DONE
	void setOutType(bool);			//DONE
	void setForce();				//Under Construction
	
private:
    
    const double pi {3.14159265358979323846};
    static const int maxGridSize {1500};		// real-time limit 3000 points approx.
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Conditions
	bool bcFlag;		// boundary condition type: 0: simply supported, 1: clamped
	bool outtype;		// output type: 0: displacement, 1: velocity
	short losstype;		// loss type: 1: independant, 2: dependant
	short itype;		// type of input: 1: struck, 2: plucked
	bool setupFlag;		// Confirm the scheme has been initialised.
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// simulation
	
	double SR;
	// physical parameters
	
	//////// String
	//         E	  nu	Rho
	// Steel : 2e11   0.30	8050
	// Alum  : 7e10   0.35	2700
	// Lead  : 1.6e10 0.44	11340
	// Wood  : 1e10   0.40	480
	
	double gauge;		// string gauge
	double f0 ;			// frequency of note in Hz (see function at EOF)
	
	double E ;			// Young's modulus (Pa) (GPa = 1e9 Pa) of steel
	double nu;			// Poisson Ratios (< .5)
	double rho;			// density (kg/m^3) of steel
	
	double r ;			// string radius (m)
	double L ;          // length (m)
	double loss [4] ;	// loss [freq.(Hz), T60;...]
	
	double T;  // Tension in Newtons
	
	//////// // I/O string
	double xi;		// coordinate of excitation (normalised, 0-1)
	double xo;		// coordinate of readout (normalised, 0-1)

	//////// // Excitation
	double famp;	// peak amplitude of excitation (N)
	double dur;		// duration of excitation (s)
	double exc_st;	// start time of excitation (s)
	double u0, v0;	// initial conditions
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double A, I, hmin, h, K, k, c, lambda, mu;
	int N, li, lo;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Allocate Memory
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Scheme States
	double uDATA[maxGridSize], u1DATA[maxGridSize], u2DATA[maxGridSize];
	double * u = uDATA, * u1 = u1DATA, * u2 = u2DATA;
	double *dummy_ptr;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double sigma0 ,sigma1, z1, z2;
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficient
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// coefficients are named based on position on the x and y axes.
	double A00, B00, B01, B02, BC1, C00, C01;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Excitation Force
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double dist, ind, rc, X, Y;
	
};

#endif /* StringClass_hpp */
