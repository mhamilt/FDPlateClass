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
public:
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
	 @param bctype boundary condition type false: Simply Supported true: Clamped
	 */
	void setup (double sampRate, bool bctype);
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
	void setOutput (double xcoord, double ycoord);
	/**
	 <#Description#>

	 @param lxcoord <#lxcoord description#>
	 @param lycoord <#lycoord description#>
	 */
	void setStereoOutput (double lxcoord, double lycoord);
	/**
	 <#Description#>

	 @param outtype <#outtype description#>
	 */
	void setOutType (bool outtype);
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
    double getOutput (bool outType);
    
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
    void getStereoOutput (bool outType, double &Left, double &Right);
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
	 <#Description#>

	 @param double <#double description#>
	 @return <#return value description#>
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
    void setCoefs (bool bcType);
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
	bool bcFlag;		// set boundary condition (s);
    /***/
	bool outFlag;		// set output type 0: displacement, 1: velocity
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
