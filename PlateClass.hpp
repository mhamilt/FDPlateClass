//
//  PlateClass.hpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//
#ifndef PlateClass_hpp
#define PlateClass_hpp

#include <iostream>			// std::cout
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#define pi 3.14159265358979323846
#define max_grid_size 1500		// real-time limit 3000 points approx.

// Quick Signum Function
int sgn(double d){
	if(d<=0){
		return 0;
	}
	else{
		return 1;
	}
}

class FD_Plate {
public:
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constructors/Assignments
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//Still Not sure About the Correct Syntax for these
	
	//Contructor
	FD_Plate();
	
	//Destructor
	~FD_Plate(); //{delete EVERYTHING;}
	
	// Copy Assignment
	//	FD_Plate& operator= (const FD_Plate&);
	
	// Move-Constructor
	//	FD_Plate (FD_Plate&&);
	
	// Move-Assignment
	//	FD_Plate& operator= (FD_Plate&&);
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Methods
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void Setup(double, bool);	//Semi-Finished
	void setLoss();					//DONE
	void setCoefs(bool);			//DONE
	void setGrid();
	void setForce();				//Under Construction
	void setOutput(float, float);	//DONE
	void setOutType(bool);			//DONE
	void UpdateScheme();			//DONE
	void ProcessIO();
	void printInfo();				//DONE
	void printCoefs();				//DONE
	double getOutput(bool);			//DONE
	
private:
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int lossFlag;		// record loss conditino type
	bool bcFlag;		// set boundary condition(s);
	bool outFlag;		// set output type 0: displacement, 1: velocity
	bool setupFlag;		// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// // wood
	double E;							// Young's modulus
	double rho;							// density (kg/m^3)
	
	double H;							// thickness (m)
	double Lx;							// x-axis plate length (m)
	double Ly;							// y-axis plate length (m)
	double loss[4];						// loss [freq.(Hz), T60;...]
	double nu;							// Poisson Ratios (< .5)
	
	// I/O Parameters
	double rp [4]; // readout position as percentage.
	
	//Excitation
	double u0 , v0, wid, ctr[2];	// excitation displacement and velocity
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	double D, kappa, hmin, h, mu, k, SR;
	int Nx, Ny, ss, li, lo;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Allocate Memory
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Scheme States
	double uDATA[max_grid_size], u1DATA[max_grid_size], u2DATA[max_grid_size];
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
	double A00, B00, B01, B11, B02, BC1, BC2, C00, C01;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Excitation Force
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double dist, ind, rc, X, Y;
	
	
};


#endif /* PlateClass_hpp */
