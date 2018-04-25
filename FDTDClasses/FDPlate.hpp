//==============================================================================
//  PlateClass.hpp
//  FDTDCPlate
//
//  Created by matthew hamilton on 15/10/2017.
//
//
//  A Finite Difference model for a linear plate based on Chapter 12 of
//  Numerical Sound Synthesis by Stefan Bilbao.
//
//  This class is still under construction
//
//==============================================================================

#ifndef PlateClasshpp
#define PlateClasshpp
#include <iostream>
#include <cmath>

class FDPlate
{
public: // Class Enums
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
public: // Methods
	//==============================================================================
	// Constructors/Assignments
	//==============================================================================
    /** Constructor*/
    FDPlate();
	/** Destructor*/
	~FDPlate(){};
	
	/// Copy Assignment
	//	FDPlate& operator= (const FDPlate&);
	
	/// Move-Constructor
	//	FDPlate (FDPlate&&);
	
	/// Move-Assignment
	//	FDPlate& operator= (FDPlate&&);
	
	//==============================================================================
	// Methods
	//==============================================================================
	/**
	 setup the plate with a given sample rate and boundary condition type: still under construction

	 @param sampRate Sample Rate in Hz
	 @param bctype boundary condition type FDPlate::BoundaryCondition
	 */
	void setup (double sampRate, BoundaryCondition bctype);
    //==============================================================================
	/**
	 <#Description#>

	 @param lowT60 <#lowT60 description#>
	 @param highT60percent <#highT60percent description#>
	 */
	void setLoss (double lowT60, double highT60percent);
	/**
	 set the plate to an initial condition of a raised cosine. This will overwrite
     all current values held on the plate.
	 */
	void setInitialCondition();
    //==============================================================================
	/**
	 <#Description#>

	 @param xcoord <#xcoord description#>
	 @param ycoord <#ycoord description#>
	 */
	void setOutputPosition (double xcoord, double ycoord);
	/**
	 <#Description#>

	 @param lxcoord <#lxcoord description#>
	 @param lycoord <#lycoord description#>
	 */
	void setStereoOutputPosition (double lxcoord, double lycoord);
	/**
	 <#Description#>

	 @param outType <#outtype description#>
	 */
	void setOutType (OutputMethod outType);
	/**
	 <#Description#>

	 @param xCoord <#xCoord description#>
	 @param yCoord <#yCoord description#>
	 */
	void setInterpOut (const double xCoord, const double yCoord);
    //==============================================================================
    /**
     <#Description#>
     
     @param outType <#outType description#>
     @return <#return value description#>
     */
    double getOutput (OutputMethod outType);
    
    /**
     <#Description#>
     
     @return <#return value description#>
     */
    double getInterpOut();
    /**
     <#Description#>
     
     @param outType <#outType description#>
     @param Left <#Left description#>
     @param Right <#Right description#>
     */
    void getStereoOutput (OutputMethod outType, double &Left, double &Right);
    /**
     <#Description#>
     
     @param force <#force description#>
     @return <#return value description#>
     */
    double reverb (double force);
    //==============================================================================
	/**
	 <#Description#>
	 */
	void updateScheme();
	/**
	 <#Description#>

	 @param force force in Newtons
	 */
	void addForce (double force);
	/**
	 Under construction: will add a strike to the plate.
	 */
	void addStrike();
    //==============================================================================
	/**
	 <#Description#>
	 */
	void processIO();
    //==============================================================================
	/**
	 <#Description#>
	 */
	void printInfo();
	/**
	 <#Description#>
	 */
	void printCoefs();
    //==============================================================================
    /**
     Signum Function

     @param double input
     @return returns 0 if d is `<` 0 or 1 if d is `>` 0
     */
	int sgn (double);
	
private: /// Methods
    //==============================================================================
	/**
	 <#Description#>

	 @return <#return value description#>
	 */
	double** getInterpLookTable();
    //==============================================================================
    /**
     <#Description#>
     
     @param bcType <#bcType description#>
     */
    void setCoefs (BoundaryCondition bcType);
    /**
     <#Description#>
     */
    void setGrid();

public: // Variables
    //==============================================================================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Allocate Memory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
    double * u, * u1, * u2;
    /***/
    double *dummyptr;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Derived Parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
    int Nx, Ny, ss, li, lo, lol, lor;


private: // Variables
    //==============================================================================
    /***/
	static constexpr const int interpOrder = 4;
    /***/
	static constexpr const int interpRes = 1000;
    /***/
	int interpZeroIndex;
    /***/
	double interpPointY, interpAlphaY;
    /***/
	double interpPointX, interpAlphaX;
    /***/
	double linearAlphas[4];
    /***/
	int xAlphaIndex, yAlphaIndex;
    /***/
	int  linearInterpInds[4];
    /***/
	int* xInterpIndeces = new int[interpOrder];
    /***/
	int* yInterpIndeces = new int[interpOrder];
    /***/
	double** interpLookTable = getInterpLookTable();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constants
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
	static constexpr double pi {3.14159265358979323846};
    /***/
	static constexpr int maxgridsize {3000};		// real-time limit 3000 points approx.
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
	BoundaryCondition currentBoundCon;		// set boundary condition (s);
    /***/
	OutputMethod outputType;		// set output type 0: displacement, 1: velocity
    /***/
	bool setupFlag;		// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// // wood
    /**Young's modulus*/
	double E;
    /**density (kg/m^3)*/
	double rho;
	
    /**thickness (m)*/
	double H;
    /**x-axis plate length (m)*/
	double Lx;
    /**y-axis plate length (m)*/
	double Ly;
    /**loss [freq.(Hz), T60;...]*/
	double loss[4];
    /**Poisson Ratios (< .5)*/
	double nu;
    
	// I/O Parameters
    /**readout position as percentage.*/
	double rp[4];
	
	//Excitation
    /***/
	double u0 , v0, wid, ctr[2]; // excitation displacement and velocity
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
	double D, kappa, hmin, h, mu, k, SR, readcoordx,readcoordy, readoutpos;
    
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
	double sigma0 ,sigma1, z1, z2;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficient
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// coefficients are named based on position on the x and y axes.
    /***/
	double A00, B00, B01, B11, B02, BC1, BC2, C00, C01, d0;
	
};


#endif /* PlateClasshpp */
