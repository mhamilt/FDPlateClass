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
	
	bcFlag_ = 0;		// set boundary condition(s);
	outFlag_ = 1;		// set output type 0: displacement, 1: velocity
	setupFlag_ = false;	// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//         E	  nu	Rho
	// Steel : 2e11   0.30	8050
	// Alum  : 7e10   0.35	2700
	// Lead  : 1.6e10 0.44	11340
	// Wood  : 1e10   0.40	480
	
	// // wood
	E_ = 11e9;						// Young's modulus
	rho_ = 480;						// density (kg/m^3)
	nu_ = .5;						// Poisson Ratios (< .5)
	
	H_ = .005;						// thickness (m)
	Lx_ = 1;							// x-axis plate length (m)
	Ly_ = 1;							// y-axis plate length (m)
	loss_[0] = 100; loss_[1] = 8;
	loss_[2] = 1000; loss_[3] = 1;    // loss [freq.(Hz), T60;...]
	
	// I/O Parameters
	rp_[0] = .45; rp_[1]=.65; rp_[2] = .85; rp_[3]= .15; // readout position as percentage.
	
	//Excitation
	ctr_[0] = .45; ctr_[1] = .45;	// centre point of excitation as percentage
	wid_ = .5;					// width (m)
	u0_ = 0; v0_ = 1;				// excitation displacement and velocity
	
}


FD_Plate::~FD_Plate()
{
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setLoss(){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	z1_ = 2*kappa_*(2*pi*loss_[0])/(2*pow(kappa_,2));
	z2_ = 2*kappa_*(2*pi*loss_[2])/(2*pow(kappa_,2));
	
	sigma0_ = 6*log(10)*(-z2_/loss_[1] + z1_/loss_[3])/(z1_-z2_);
	sigma1_ = 6*log(10)*(1/loss_[1] - 1/loss_[3])/(z1_-z2_);
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setGrid(){
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Spacing
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// stability condition
	hmin_ = (sqrt(4*k_*(sigma1_+sqrt(pow(sigma1_,2)+pow(kappa_,2)))));
	
	Nx_ = floor(Lx_/hmin_);		// number of segments x-axis
	Ny_ = floor(Ly_/hmin_);		// number of segments y-axis
	h_ = sqrt(Lx_*Ly_/(Nx_*Ny_));;	// adjusted grid spacing x/y
	Nx_ = Nx_+1; Ny_ = Ny_+1;		// grid point number x and y
	mu_ = (kappa_ * k_)/pow(h_,2);	// scheme parameter
	ss_ = Nx_*Ny_;					// total grid size.
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setCoefs(bool bcType = 0){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//update flag
	bcFlag_ = bcType;
	
	// coefficients are named based on position on the x and y axes.
	A00_ = 1/(1+k_*sigma0_); // Central Loss Coeeffient (INVERTED)
	
	//// Current time step (B) coeffients
	// There are six unique coefficients for B coefs
	B00_ = (-pow(mu_,2)*20 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_;	// center
	B01_ = (-pow(mu_,2)*-8 + (2*sigma1_*k_/pow(h_,2))) * A00_;		// 1-off
	B11_ = (-pow(mu_,2)*2) * A00_;									// diag
	B02_ = (-pow(mu_,2)*1) * A00_;									// 2-off
	
	if(bcType){ // Clamped Boundary Coefficients
		BC1_ = (-pow(mu_,2)*21 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Side
		BC2_ = (-pow(mu_,2)*22 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Corner
	}
	else { // Simply Supported Boundary Coefficients
		BC1_ = (-pow(mu_,2)*19 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Side
		BC2_ = (-pow(mu_,2)*18 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Corner
	}
	
	//// Previous time step (C) coeffients
	C00_ = (-(2*sigma1_*k_/pow(h_,2))*-4 - (1-sigma0_*k_))  * A00_;
	C01_ = -(2*sigma1_*k_/pow(h_,2))  * A00_;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Setup Scheme
// Takes sample rate as an arguement, 44.1kHz is set by default
void FD_Plate::Setup(double samprate = 44.1e3, bool bctype = 0){
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Motion Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	D_ = (E_*(pow(H_, 3)))/(12*(1-pow(nu_,2)));
	kappa_ = sqrt(D_ / (rho_*  H_) );
	
	SR_ = samprate;				// internal class sampling rate
	k_ = 1/SR_;					// time step
	
	setLoss();
	setGrid();
	setCoefs(bctype);
	
	//Set Input and Output Indeces
	li_ = (Ny_*( ( (ctr_[1]*Nx_)-1)) ) + (ctr_[0]*Ny_)-1;
	lo_ = (Ny_*( ( (rp_[1]*Nx_)-1)) ) +  (rp_[0]*Ny_)-1;
	
	//	Update flags
	setupFlag_ = true;
	bcFlag_  = bctype;
	
	//	Clear State Memory
	memset(uDATA_, 0, ss_ * sizeof(double));
	memset(u1DATA_, 0, ss_ * sizeof(double));
	memset(u2DATA_, 0, ss_ * sizeof(double));
	
	setForce();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Print Scheme Info
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::printCoefs(){
	
	printf("--- Coefficient Info --- \n\n");
	printf("Loss A		: %.4fm \n", A00_);
	printf("Centre B    : %.4fm \n", B00_);
	printf("1-Grid B    : %.4fm \n", B01_);
	printf("2-Grid B	: %.4fm \n", B02_);
	printf("Diagonal B  : %.4fm \n", B11_);
	printf("Centre C	: %.4fm \n", C00_);
	printf("1-Grid C    : %.4fm \n", C01_);
	printf("Side Bound	: %.4fm \n", BC1_);
	printf("Cornr Bound : %.4fm \n\n", BC2_);
}

void FD_Plate::printInfo(){
	printf("--- Scheme Info --- \n\n");
	printf("Size			: %.1f m2 \n", Nx_*h_*Ny_*h_);
	printf("Thickness(mm)   : %.0f mm \n", H_*1000);
	printf("Grid X-Ax		: %d \n", Nx_);
	printf("Grid Y-Ax		: %d \n", Ny_);
	printf("Total Ps		: %d \n", ss_);
	printf("In_cell			: %d\n", li_);
	printf("Out_cell		: %d\n", lo_);
	printf("TimeStep		: %f\n", k_);
	printf("SampRate		: %.0f\n", SR_);
	printf("Youngs			: %.2e\n", E_);
	printf("Sigma 0			: %f\n", sigma0_);
	printf("Sigma 1			: %f\n", sigma1_);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update Plate State
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// TODO: update should not be run before Setup()
// Do a check to ensure

void FD_Plate::UpdateScheme(){
	
	int xi, yi, cp;

	// Internal Gride Points
	for(xi=Nx_-4; xi--; ){
		
		for(yi=Ny_-4;yi--; ){
			
			cp = (yi+2)+((xi+2) * Ny_); // current point
			
			u_[cp] = B00_*u1_[cp] +
			B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
			B02_*( u1_[cp-2] + u1_[cp+2] +u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
			B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
			C00_*u2_[cp] +
			C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
		}
	}
	
	// Update Side Boundaries
	//X-Axis
	
	for(xi=Nx_-4; xi--; ){
		//North
		cp = 1+((xi+2) * Ny_); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] + u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp+1-Ny_] + u1_[cp+1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
		
		//South
		cp = Ny_-2 +((xi+2) * Ny_); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp-1-Ny_] + u1_[cp-1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
		
		
	}
	
	// Y-Axis
	
	for(yi=Ny_-4;yi--; ){
		//West
		cp = yi+Ny_+2; // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp+Ny_] );
		
		//East
		cp = (yi+2) + Ny_*(Nx_-2); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] +u1_[cp-(2*Ny_)] ) +
		B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] );
		
	}
	
	// Corner Boundaries
	
	cp = Ny_+1;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp+2] + u1_[cp+(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = 2*(Ny_-1);
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp-2] + u1_[cp+(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = Ny_*(Nx_-2)+1;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp+2] + u1_[cp-(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = Ny_*(Nx_-1) - 2;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp-2] + u1_[cp-(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	// swap pointers
	dummy_ptr_ = u2_; u2_ = u1_; u1_ = u_; u_ = dummy_ptr_;
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// get the output from the plate based on the readout positino and type.
// Will need to implement interpolation, especially f this is meant to be a
// free moving read-out.

double FD_Plate::getOutput(bool outType){
	
	if(outType){
		return (u1_[lo_]- u2_[lo_])*SR_; // Velocity out
	}
	else{
		return u1_[lo_]; // Amplitude out
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// Set the readout position on the plate

void FD_Plate::setOutput(float xcoord, float ycoord){
	
	//TODO: CHECK IF ON A VALID GRID POINT
	if((xcoord*Nx_)-1 < 1 || (xcoord*Nx_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	
	if((ycoord*Ny_)-1 < 1 || (ycoord*Ny_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	
	lo_ = (Ny_*( ( (xcoord*Nx_)-1)) ) +  (ycoord*Ny_)-1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FD_Plate::setOutType(bool outtype){ // set output to velocity amplitude
	outFlag_ = outtype;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


// Method will set the profile of the input force, for use in JUCE.
// Currently not interpolated, will need to look into that.
// Array values will them selves be multiplied by a rasied cosine/half-cosine

void FD_Plate::setForce(){
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Excitation Force
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double dist, ind, rc, X, Y;
	int xi, yi, cp;  //loop indeces
	
	// raised cosine in 2D
	for(xi=1;xi<Nx_-1;xi++){
		
		X = xi*h_;
		
		for(yi=1;yi<Ny_-1;yi++){
			cp = yi+(xi * Ny_);
			
			Y = yi*h_;
			
			dist = sqrt(pow(X-(ctr_[0]*Lx_),2) + pow(Y-(ctr_[1]*Ly_),2));
			
			ind = sgn((wid_*0.5)-dist);			// displacement (logical)
			rc = .5*ind*(1+cos(2*pi*dist/wid_)); // displacement
			
			u2_[cp] = u0_*rc;
			u1_[cp] = v0_*k_*rc;
			
		}
	}
}
