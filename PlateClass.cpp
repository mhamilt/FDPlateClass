//
//  PlateClass.cpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//  Class file for a FDTD Plate

#include "PlateClass.hpp"

FD_Plate::FD_Plate()
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	bcFlag = 0;		// set boundary condition(s);
	outFlag = 1;		// set output type 0: displacement, 1: velocity
	setupFlag = false;	// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//         E	  nu	Rho
	// Steel : 2e11   0.30	8050
	// Alum  : 7e10   0.35	2700
	// Lead  : 1.6e10 0.44	11340
	// Wood  : 1e10   0.40	480
	
	// // wood
	E = 11e9;						// Young's modulus
	rho = 480;						// density (kg/m^3)
	nu = .5;						// Poisson Ratios (< .5)
	
	H = .005;						// thickness (m)
	Lx = 1;							// x-axis plate length (m)
	Ly = 1;							// y-axis plate length (m)
	loss[0] = 100; loss[1] = 2;
	loss[2] = 1000; loss[3] = 1;    // loss [freq.(Hz), T60;...]
	
	// I/O Parameters
	rp [0] = .45; rp[1]=.65; rp[2] = .85; rp[3]= .15; // readout position as percentage.
	
	//Excitation
	ctr[0] = .45; ctr[1] = .45;	// centre point of excitation as percentage
	wid = .5;					// width (m)
	u0 = 0; v0 = 1;				// excitation displacement and velocity
	
}


FD_Plate::~FD_Plate()
{
	
}

// set the type of loss and update the lossFlag
void FD_Plate::setLoss(){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	z1 = 2*kappa*(2*pi*loss[0])/(2*pow(kappa,2));
	z2 = 2*kappa*(2*pi*loss[2])/(2*pow(kappa,2));
	
	sigma0 = 6*log(10)*(-z2/loss[1] + z1/loss[3])/(z1-z2);
	sigma1 = 6*log(10)*(1/loss[1] - 1/loss[3])/(z1-z2);
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setGrid(){
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Spacing
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// stability condition
	hmin = (sqrt(4*k*(sigma1+sqrt(pow(sigma1,2)+pow(kappa,2)))));
	
	Nx = floor(Lx/hmin);		// number of segments x-axis
	Ny = floor(Ly/hmin);		// number of segments y-axis
	h = sqrt(Lx*Ly/(Nx*Ny));;	// adjusted grid spacing x/y
	Nx = Nx+1; Ny = Ny+1;		// grid point number x and y
	mu = (kappa * k)/pow(h,2);	// scheme parameter
	ss = Nx*Ny;					// total grid size.
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setCoefs(bool bcType = 0){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//update flag
	bcFlag = bcType;
	
	// coefficients are named based on position on the x and y axes.
	A00 = 1/(1+k*sigma0); // Central Loss Coeeffient (INVERTED)
	
	//// Current time step (B) coeffients
	// There are six unique coefficients for B coefs
	B00 = (-pow(mu,2)*20 + (2*sigma1*k/pow(h,2))*-4 + 2) * A00;	// center
	B01 = (-pow(mu,2)*-8 + (2*sigma1*k/pow(h,2))) * A00;		// 1-off
	B11 = (-pow(mu,2)*2) * A00;									// diag
	B02 = (-pow(mu,2)*1) * A00;									// 2-off
	
	if(bcType){ // Clamped Boundary Coefficients
		BC1 = (-pow(mu,2)*21 + (2*sigma1*k/pow(h,2))*-4 + 2) * A00; // Side
		BC2 = (-pow(mu,2)*22 + (2*sigma1*k/pow(h,2))*-4 + 2) * A00; // Corner
	}
	else { // Simply Supported Boundary Coefficients
		BC1 = (-pow(mu,2)*19 + (2*sigma1*k/pow(h,2))*-4 + 2) * A00; // Side
		BC2 = (-pow(mu,2)*18 + (2*sigma1*k/pow(h,2))*-4 + 2) * A00; // Corner
	}
	
	//// Previous time step (C) coeffients
	C00 = (-(2*sigma1*k/pow(h,2))*-4 - (1-sigma0*k))  * A00;
	C01 = -(2*sigma1*k/pow(h,2))  * A00;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Setup Scheme
// Takes sample rate as an arguement, 44.1kHz is set by default
void FD_Plate::Setup(double samprate = 44.1e3, bool bctype = 0){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Motion Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	D = (E*(pow(H, 3)))/(12*(1-pow(nu,2)));
	kappa = sqrt(D / (rho*  H) );
	
	SR = samprate;				// internal class sampling rate
	k = 1/SR;					// time step
	
	setLoss();
	setGrid();
	setCoefs(bctype);
	
	//Set Input and Output Indeces
	li = (Ny*( ( (ctr[1]*Nx)-1)) ) + (ctr[0]*Ny)-1;
	lo = (Ny*( ( (rp[1]*Nx)-1)) ) +  (rp[0]*Ny)-1;
	
	//	Update flags
	setupFlag = true;
	bcFlag  = bctype;
	
	//	Clear State Memory
	memset(uDATA, 0, ss * sizeof(double));
	memset(u1DATA, 0, ss * sizeof(double));
	memset(u2DATA, 0, ss * sizeof(double));
	
	setForce();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Print Scheme Info
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::printCoefs(){
	
	printf("--- Coefficient Info --- \n\n");
	printf("Loss A		: %.4fm \n", A00);
	printf("Centre B    : %.4fm \n", B00);
	printf("1-Grid B    : %.4fm \n", B01);
	printf("2-Grid B	: %.4fm \n", B02);
	printf("Diagonal B  : %.4fm \n", B11);
	printf("Centre C	: %.4fm \n", C00);
	printf("1-Grid C    : %.4fm \n", C01);
	printf("Side Bound	: %.4fm \n", BC1);
	printf("Cornr Bound : %.4fm \n\n", BC2);
}

void FD_Plate::printInfo(){
	printf("--- Scheme Info --- \n\n");
	printf("Size			: %.1f m2 \n", Nx*h*Ny*h);
	printf("Thickness(mm)   : %.0f mm \n", H*1000);
	printf("Grid X-Ax		: %d \n", Nx);
	printf("Grid Y-Ax		: %d \n", Ny);
	printf("Total Ps		: %d \n", ss);
	printf("In_cell			: %d\n", li);
	printf("Out_cell		: %d\n", lo);
	printf("TimeStep		: %f\n", k);
	printf("SampRate		: %.0f\n", SR);
	printf("Youngs			: %.2e\n", E);
	printf("Sigma 0			: %f\n", sigma0);
	printf("Sigma 1			: %f\n", sigma1);
	printf("LossType		: %d\n\n", lossFlag);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update Plate State
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// TODO: update should not be run before Setup()
// Do a check to ensure

void FD_Plate::UpdateScheme(){
	
	// Internal Gride Points
	int xi, yi, cp;
	
	for(xi=Nx-4; xi--; ){
		
		for(yi=Ny-4;yi--; ){
			
			cp = (yi+2)+((xi+2) * Ny); // current point
			
			u[cp] = B00*u1[cp] +
			B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
			B02*( u1[cp-2] + u1[cp+2] +u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
			B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
			C00*u2[cp] +
			C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
		}
	}
	
	// Update Side Boundaries
	//X-Axis
	
	for(xi=Nx-4; xi--; ){
		//North
		cp = 1+((xi+2) * Ny); // current point
		u[cp]  = BC1*u1[cp] +
		B01*( u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
		B02*( u1[cp-2] + u1[cp+2] +u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
		B11*( u1[cp+1-Ny] + u1[cp+1+Ny] ) +
		C00*u2[cp] +
		C01*( u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
		
		//South
		cp = Ny-2 +((xi+2) * Ny); // current point
		u[cp]  = BC1*u1[cp] +
		B01*( u1[cp-1] + u1[cp-Ny] + u1[cp+Ny] ) +
		B02*( u1[cp-2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
		B11*( u1[cp-1-Ny] + u1[cp-1+Ny] ) +
		C00*u2[cp] +
		C01*( u2[cp-1] + u2[cp-Ny] + u2[cp+Ny] );
		
		
	}
	
	// Y-Axis
	
	for(yi=Ny-4;yi--; ){
		//West
		cp = yi+Ny+2; // current point
		u[cp]  = BC1*u1[cp] +
		B01*( u1[cp-1] + u1[cp+1] + u1[cp+Ny] ) +
		B02*( u1[cp-2] + u1[cp+2] + u1[cp+(2*Ny)] ) +
		B11*( u1[cp+1+Ny] + u1[cp-1+Ny] ) +
		C00*u2[cp] +
		C01*( u2[cp-1] + u2[cp+1] + u2[cp+Ny] );
		
		//East
		cp = (yi+2) + Ny*(Nx-2); // current point
		u[cp]  = BC1*u1[cp] +
		B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] ) +
		B02*( u1[cp-2] + u1[cp+2] +u1[cp-(2*Ny)] ) +
		B11*( u1[cp-1-Ny] + u1[cp+1-Ny] ) +
		C00*u2[cp] +
		C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] );
		
	}
	
	// Corner Boundaries
	
	cp = Ny+1;
	u[cp] = BC2*u1[cp] +
	B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
	B02*( u1[cp+2] + u1[cp+(2*Ny)] ) +
	B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
	C00*u2[cp] +
	C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
	
	cp = 2*(Ny-1);
	u[cp] = BC2*u1[cp] +
	B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
	B02*( u1[cp-2] + u1[cp+(2*Ny)] ) +
	B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
	C00*u2[cp] +
	C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
	
	cp = Ny*(Nx-2)+1;
	u[cp] = BC2*u1[cp] +
	B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
	B02*( u1[cp+2] + u1[cp-(2*Ny)] ) +
	B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
	C00*u2[cp] +
	C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
	
	cp = Ny*(Nx-1) - 2;
	u[cp] = BC2*u1[cp] +
	B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
	B02*( u1[cp-2] + u1[cp-(2*Ny)] ) +
	B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
	C00*u2[cp] +
	C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
	
	// swap pointers
	dummy_ptr = u2; u2 = u1; u1 = u; u = dummy_ptr;
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// get the output from the plate based on the readout positino and type.
// Will need to implement interpolation, especially f this is meant to be a
// free moving read-out.

double FD_Plate::getOutput(bool outType){
	
	if(outType){
		return (u1[lo]- u2[lo])*SR; // Velocity out
	}
	else{
		return u1[lo]; // Amplitude out
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// Set the readout position on the plate

void FD_Plate::setOutput(float xcoord, float ycoord){
	
	//TODO: CHECK IF ON A VALID GRID POINT
	if((xcoord*Nx)-1 < 1 || (xcoord*Nx)-1 > Nx-2){
		//Do something to ensure it is not on a zero point
	}
	
	if((ycoord*Ny)-1 < 1 || (ycoord*Ny)-1 > Nx-2){
		//Do something to ensure it is not on a zero point
	}
	
	lo = (Ny*( ( (xcoord*Nx)-1)) ) +  (ycoord*Ny)-1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FD_Plate::setOutType(bool outtype){ // set output to velocity amplitude
	outFlag = outtype;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


// Method will set the profile of the input force, for use in JUCE.
// Currently not interpolated, will need to look into that.
// Array values will them selves be multiplied by a rasied cosine/half-cosine

void FD_Plate::setForce(){
	
	int xi, yi, cp;  //loop indeces
	
	// raised cosine in 2D
	for(xi=1;xi<Nx-1;xi++){
		
		X = xi*h;
		
		for(yi=1;yi<Ny-1;yi++){
			cp = yi+(xi * Ny);
			
			Y = yi*h;
			
			dist = sqrt(pow(X-(ctr[0]*Lx),2) + pow(Y-(ctr[1]*Ly),2));
			
			ind = sgn((wid*0.5)-dist);			// displacement (logical)
			rc = .5*ind*(1+cos(2*pi*dist/wid)); // displacement
			
			u2[cp] = u0*rc;
			u1[cp] = v0*k*rc;
			
		}
	}
}
