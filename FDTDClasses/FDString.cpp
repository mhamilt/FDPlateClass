//
//  StringClass.cpp
//  FDTD_C_String
//
//  Created by admin on 07/11/2017.
//  Copyright © 2017 admin. All rights reserved.
//

#include "FDString.hpp"

// Quick Signum Function
int sgn(double d){
	if(d<=0){
		return 0;
	}
	else{
		return 1;
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FDString Constructor
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FDString::FDString()
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// START EDIT HERE
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Conditions
	bcFlag = 0;		// boundary condition type: 0: simply supported, 1: clamped
	outtype = 0;		// output type: 0: displacement, 1: velocity
	losstype = 2;		// loss type: 1: independant, 2: dependant
	itype = 2;		// type of input: 1: struck, 2: plucked
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// simulation
	SR = 44.1e3;
	
	// physical parameters
	
	//////// String
	//         E	  nu	Rho
	// Steel : 2e11   0.30	8050
	// Alum  : 7e10   0.35	2700
	// Lead  : 1.6e10 0.44	11340
	// Wood  : 1e10   0.40	480
	
	gauge = 25;		// string gauge
	f0 = 440;		// frequency of note in Hz (see function at EOF)
	
	E = 2e11;                    // Young's modulus (Pa) (GPa = 1e9 Pa) of steel
	nu = .5;							// Poisson Ratios (< .5)
	rho = 7850;                  // density (kg/m^3) of steel
	
	r = (gauge * 2.54e-5)*.5;	// string radius (m)
	L = .6;                      // length (m)
	loss[0] = 100; loss[1] = 8;
	loss[2] = 1000; loss[3] = 1;    // loss [freq.(Hz), T60;...]
	
	T = pow(((2*f0*r)*L),2)*pi*rho; // Tension in Newtons
	
	// // I/O
	xi = 0.75;		// coordinate of excitation (normalised, 0-1)
	xo = 0.25;		// coordinate of readout (normalised, 0-1)
	
	// // Excitation
	famp = 10;		// peak amplitude of excitation (N)
	dur = 0.001;	// duration of excitation (s)
	exc_st = 0;		// start time of excitation (s)
	u0 = 0; v0 = 1;   // initial conditions
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// END EDIT HERE
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// String Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// derived parameters
	A = pi*pow(r,2);                       // string cross-sectional area
	I = 0.25*pi*pow(r,4);                  // string moment of intertia
	
	c = sqrt(T/(rho*A));        // wave speed
	K = sqrt(E*I/(rho*A));   // stiffness constant
	k = 1/SR;
	
}


FDString::~FDString()
{
	
}

// set the type of loss and update the lossFlag
void FDString::setLoss(int lossType){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	switch (losstype) {
			
  case 0:	// Lossless
			sigma0 = 0;
			sigma1 = 0;
			break;
			
  case 1:	// frequency independant loss
			sigma0 = 6*log(10)/loss[1];
			sigma1 = 0;
			break;
			
  case 2:	// frequency dependant loss
			z1 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*pow(2*pi*loss[0],2)))/(2*pow(K,2));
			z2 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*pow(2*pi*loss[2],2)))/(2*pow(K,2));
			sigma0 = 6*log(10)*(-z2/loss[1] + z1/loss[3])/(z1-z2);
			sigma1 = 6*log(10)*(1/loss[1] - 1/loss[3])/(z1-z2);
			break;
			
  default:	// frequency dependant loss
			z1 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*(2*pi*pow(loss[0],2))))/(2*pow(K,2));
			z2 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*(2*pi*pow(loss[2],2))))/(2*pow(K,2));
			sigma0 = 6*log(10)*(-z2/loss[1] + z1/loss[3])/(z1-z2);
			sigma1 = 6*log(10)*(1/loss[1] - 1/loss[3])/(z1-z2);
			
	}
	
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FDString::setCoefs(bool bcType = 0){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Central Loss Coeffient A (INVERTED)
	A00 = 1/(1+k*sigma0);
	
	// Current time step (B) coeffients
	B00 = (-pow(lambda,2)*2 - pow(mu,2)*6 + 2 - (2*sigma1*k/pow(h,2))*2) * A00; // centre
	B01 = (pow(lambda,2) + pow(mu,2)*4 + (2*sigma1*k/pow(h,2)) ) * A00;	// 1-off
	B02 =  (-pow(mu,2)) * A00;	// 2-off
	
	// Simply Supported Boundary Coefficients
	
	if(bcType){
		BC1 = (pow(lambda,2)*-2 - pow(mu,2)*5 + 2 + (2*sigma1*k/pow(h,2))*-2) * A00;
	}
	else {
		BC1 = (pow(lambda,2)*-2 - pow(mu,2)*7 + 2 + (2*sigma1*k/pow(h,2))*-2) * A00;
	}
	
	//// Previous time step (C) coeffients
	C00 = - ( (1-sigma0*k) + ((2*sigma1*k/pow(h,2))*-2) )  * A00;
	C01 = -((2*sigma1*k/pow(h,2)))  * A00;
	
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


// Set up grid spacing of scheme

void FDString::setGrid(){
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// String Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// derived parameters
	
	A = pi*pow(r,2);			// string cross-sectional area
	I = 0.25*pi*pow(r,4);		// string moment of intertia
	c = sqrt(T/(rho*A));		// wave speed
	K = sqrt(E*I/(rho*A));		// stiffness constant
	k = 1/SR;					// Time Step
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Grid
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//	double hmin = sqrt(0.5* (pow(c,2)*pow(k,2)+sqrt(pow(c,4)*pow(k,4)+16*pow(K,2)*pow(k,2))) );
	
	//// spacing for frequency dependant loss
	hmin = sqrt(pow(c*k,2)+4*sigma1*k + sqrt( pow(pow(c*k,2)+4*sigma1*k,2) +16*pow(K*k,2)));
	
	N = floor(L/hmin);      // number of segments (N+1 is number of grid points)
	h = L/N;                // adjusted grid spacing
	lambda = c*k/h;         // Courant number
	mu = K*k/pow(h,2);      // numerical stiffness constant
	N = N+1;				// change N to be number of grid points
	
}

// Setup Scheme
// Takes sample rate as an arguement, 44.1kHz is set by default
void FDString::Setup(double sampRate = 44.1e3, int losstype = 2, bool bcType = 0){
	
	// configure loss, setCoefs() also run in setLoss()
	setLoss(losstype);
	
	//	Set Scheme Grid Spacing
	setGrid();
	
	// Clear State Memory
	memset(uDATA, 0, N * sizeof(double));
	memset(u1DATA, 0, N * sizeof(double));
	memset(u2DATA, 0, N * sizeof(double));
	
	// Set Scheme Coefficients
	setCoefs(bcType);
	
	// I/O
	setInput(xi);
	setOutput(xo);
	setForce();
	
	//update flags
	setupFlag = true;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Force
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//	double force[Nf];									// force array
	//
	//	// Zero the array
	//	memset(force, 0, Nf * sizeof(double));
	//
	//
	//	double d0 = pow(k,2)/(h*rho*A*(1+k*sigma0))*A00;	// force coeffcient
	//	int durint = floor(dur*SR);			// duration of force signal, in samples
	//	int exc_int = (floor(exc_st*SR));	// start time index for excitation
	//
	//	for(n = exc_int;n<(exc_int+durint-1); n++){
	//
	//		// reevaluate relevant points in force vector
	//		force[n] = famp*0.5*(1-cos((2/itype)*pi*(n/double(durint))))*d0;
	//	}
	
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Print Scheme Info
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FDString::printCoefs(){
	
	printf("--- Coefficient Info --- \n\n");
	printf("Loss A		: %.7fm \n", A00);
	printf("Centre B    : %.6fm \n", B00);
	printf("1-Grid B    : %.6fm \n", B01);
	printf("2-Grid B	: %.6fm \n", B02);
	printf("Centre C	: %.6fm \n", C00);
	printf("1-Grid C    : %.6fm \n", C01);
	printf("BoundCon	: %.6fm \n", BC1);
	printf("Sigma 0		: %f\n", sigma0);
	printf("Sigma 1		: %f\n", sigma1);
	printf("z1		: %f\n", z1);
	printf("z2		: %f\n\n", z2);
}

void FDString::printInfo(){
	
	printf("\n--- Scheme Info --- \n\n");
	printf("Size		: %.1fm \n", N*h);
	printf("Gridmin		: %.7f \n", hmin);
	printf("GridSpace	: %.7f \n", h);
	printf("Grid Num	: %d \n", N);
	printf("In_cell		: %d\n", li);
	printf("Out_cell	: %d\n", lo);
	printf("Youngs		: %.2e\n\n", E);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update Plate State
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// TODO: update should not be run before Setup()
// Do a check to ensure

/*
 In the update Scheme the for loops have been expressed as a test if the
 loop index is equal to 0. This allegedley faster and so far shave a half second
 off over ten seconds. Above each loop is the equivalent not taking into account
 an extra adjustment of the index which is in the caluculation of the
 current point 'cp'.
 
 The update also contains the time loop as well. Though this has had no affect on
 performance
 */

void FDString::UpdateScheme(){
	
	// Internal Gride Points
	int i;

	// Internal Gride Points
	for(i=2; i<N-2; i++){
	//		for(i=Nx-4; i--; ){
		
		u[i] = B00*u1[i] +
		B01*( u1[i-1] + u1[i+1]) +
		B02*( u1[i-2] + u1[i+2]) +
		C00*u2[i] +
		C01*( u2[i-1] + u2[i+1]);
		
	}
	
	// Boundaries
	i = 1;
	u[i] = BC1*u1[i] +
	B01*(u1[i+1]) +
	B02*(u1[i+2]) +
	C00*u2[i] +
	C01*(u2[i+1]);
	
	i = N-2;
	u[i] = BC1*u1[i] +
	B01*( u1[i-1]) +
	B02*( u1[i-2]) +
	C00* u2[i] +
	C01*(u2[i-1]);
	
	// Pointer Swap
	dummy_ptr = u2; u2 = u1; u1 = u; u = dummy_ptr;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// get the output from the plate based on the readout positino and type.
// Will need to implement interpolation, especially f this is meant to be a
// free moving read-out.

double FDString::getOutput(bool outType){
	
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

void FDString::setOutput(float coord = 0.25){
	
	//TODO: CHECK IF ON A VALID GRID POINT
	if((coord*N)-1 < 1 || (coord*N)-1 > N-2){
		//Do something to ensure it is not on a zero point
	}
	
	lo = (coord*N)-1;
}

void FDString::setInput(float coord = 0.75){
	
	//TODO: CHECK IF ON A VALID GRID POINT
	if((coord*N)-1 < 1 || (coord*N)-1 > N-2){
		//Do something to ensure it is not on a zero point
	}
	
	li = (coord*N)-1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FDString::setOutType(bool outtype){ // set output to velocity amplitude
	//	outFlag = outtype;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


// Method will set the profile of the input force, for use in JUCE.
// Currently not interpolated, will need to look into that.
// Array values will them selves be multiplied by a rasied cosine/half-cosine

void FDString::setForce(){
	
	u2[li] = u0;
	u1[li] = u0 + (v0*k);
	
}
