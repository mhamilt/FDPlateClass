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
	bool bcFlag_;		// set boundary condition(s);
	bool outFlag_;		// set output type 0: displacement, 1: velocity
	bool setupFlag_;		// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// // wood
	double E_;							// Young's modulus
	double rho_;							// density (kg/m^3)
	
	double H_;							// thickness (m)
	double Lx_;							// x-axis plate length (m)
	double Ly_;							// y-axis plate length (m)
	double loss_[4];						// loss [freq.(Hz), T60;...]
	double nu_;							// Poisson Ratios (< .5)
	
	// I/O Parameters
	double rp_[4]; // readout position as percentage.
	
	//Excitation
	double u0_ , v0_, wid_, ctr_[2];	// excitation displacement and velocity
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	double D_, kappa_, hmin_, h_, mu_, k_, SR_;
	int Nx_, Ny_, ss_, li_, lo_;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Allocate Memory
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Scheme States
	double uDATA_[max_grid_size], u1DATA_[max_grid_size], u2DATA_[max_grid_size];
	double * u_ = uDATA_, * u1_ = u1DATA_, * u2_ = u2DATA_;
	double *dummy_ptr_;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double sigma0_ ,sigma1_, z1_, z2_;
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficient
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// coefficients are named based on position on the x and y axes.
	double A00_, B00_, B01_, B11_, B02_, BC1_, BC2_, C00_, C01_;
	

	
	
};


#endif /* PlateClass_hpp */
