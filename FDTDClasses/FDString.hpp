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

class FDString
{
public: // Enum Classes
    enum class BoundaryCondition
    {
        simplySupported,
        clamped
    };
    
    enum class OutputMethod
    {
        amplitude,
        velocity
    };
    
    enum class InputMethod
    {
        pluck,
        strike
    };
    enum class LossModel
    {
        lossless,
        simple,
        frequencyDepenent
    };
public:
	//==========================================================================
	//Contructor
	FDString();
	//==========================================================================
	//Destructor
	~FDString(); //{delete EVERYTHING;}
	//==========================================================================
	// Copy Assignment
	//	FDString& operator= (const FDString&);
	//==========================================================================
	// Move-Constructor
	//	FDString (FDString&&);
	//==========================================================================
	// Move-Assignment
	//	FDString& operator= (FDString&&);
	//==========================================================================
	void setup(double, LossModel, BoundaryCondition);   
    //==========================================================================
    void printInfo();
    void printCoefs();
    //==========================================================================
    void addForce();
    //==========================================================================
    /**
     Set the position of force input

     @param dist a normalised distance between 0 and 1 along the string
     */
    void setInputPosition (float dist);
    /**
     Set the position of reading the output

     @param dist a normalised distance between 0 and 1 along the string
     */
    void setOutputPosition (float dist);
    /**
     set the method of reading output

     @param outType The method of readin output (enum)
     */
    void setOutType (OutputMethod outType);
    //==========================================================================
    /**
     calculate the next time step state
     */
    void updateScheme();
    //==========================================================================
    /**
     read a sample value from the string grid

     @return returns the output of the string
     */
    double getOutput();
private: // functions
    //==========================================================================
    void setLoss (LossModel lossType);
    //==========================================================================
    void setGrid();
    void setMemory();
    //==========================================================================
    void setCoefs (BoundaryCondition bctype);
    //==========================================================================
    void processIO();
    //==========================================================================
    void setInitialCondition();
    //==========================================================================
    /**
     adds force to the `u` time state and updates the force vector and force index `fi`
     */
    void updateForce();
    //==========================================================================
  private:
    //==========================================================================
    const double pi {M_PI};
    /// real-time limit 3000 points approx.
    const int maxGridSize {1500};
    /// max Excitation duration in samples
    const int maxExcDur = 130;
    //==========================================================================
    ///
	OutputMethod outputType;
    ///
	InputMethod inputType;
    ///
	bool setupFlag;
	//==========================================================================
	double sampleRate;
    //==========================================================================
	double gauge;		// string gauge
	double f0 ;			// frequency of note in Hz (see function at EOF)
	double E ;			// Young's modulus (Pa) (GPa = 1e9 Pa) of steel
	double nu;			// Poisson Ratios (< .5)
	double rho;			// density (kg/m^3) of steel
	double r ;			// string radius (m)
	double L ;          // length (m)
	double loss [4] ;	// loss [freq.(Hz), T60;...]
	double T;           // Tension in Newtons
	//==========================================================================
	double xi;		// coordinate of excitation (normalised, 0-1)
	double xo;		// coordinate of readout (normalised, 0-1)
	//==========================================================================
	double famp;	// peak amplitude of excitation (N)
	double dur;		// duration of excitation (s)
	double exc_st;	// start time of excitation (s)
	double u0, v0;	// initial conditions
	//==========================================================================
	double A, I, hmin, h, K, k, c, lambda, mu;
	int N, li, lo;
	//==========================================================================
	/// Scheme States
	double *u, *u1, *u2;
	double *dummy_ptr;
	//==========================================================================
    /// Loss Coefficient
	double sigma0 ,sigma1, z1, z2;
	//==========================================================================
	/// Central Loss Coeffient A (INVERTED)
    double A00;
    /// Current time step (B) coeffients
    double B00, B01, B02;
    /// Boundary Coefficients
    double BC1;
    /// Previous time step (C) coeffients
    double C00, C01;
	//==========================================================================
    /// force vector
    double *force;
    /// force index
    int fi = 0;
};

#endif /* StringClass_hpp */
