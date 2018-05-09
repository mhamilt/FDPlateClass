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
    
    struct PlateParameters
    {
        /// Young's modulus
        float youngs = 11e9;
        /// density (kg/m^3)
        float density = 480;
        /// Poisson Ratios (< .5)
        float poisson = .5;
        /// thickness (m)
        float thickness = .003;
        /// x-axis plate length (m)
        float lengthX = 1;
        /// y-axis plate length (m)
        float lengthY = 1;
        /// T60 decay
        float t60 = 5;
        /// high frequency: percent of T60 (0 < tone < 1)
        float tone = 0.9;
        /// boundary condtions
        BoundaryCondition bcType = BoundaryCondition::simplySupported;
    };
public: // Methods
    //==========================================================================
    // Constructors/Assignments
    //==========================================================================
    /** Constructor: Initialise with sample rate and FDPlate::PlateParameter struct
     or with default settings specifying nothing or just the sample
     rate. The Plate can be re-set using the setup() method*/
    FDPlate(double sampRate, PlateParameters plateParams);
    FDPlate(double sampRate);
    FDPlate();
    /** Destructor*/
    ~FDPlate();
    
    //==========================================================================
    // Methods
    //==========================================================================
    /**
     set the plate to an initial condition of a raised cosine. This will overwrite
     all current values held on the plate.
     */
    void setInitialCondition();
    //==========================================================================
    /**
     Set the read-out point as a normalised position between 0 and 1.
     This position will be rounded to the nearest grid point
     
     @param xcoord x-axis co-ordinate
     @param ycoord y-axis co-ordinate
     */
    void setOutputPosition (double xcoord, double ycoord);
    /**
     set the read-out position for interpolated output
     
     @param xCoord x-axis co-ordinate
     @param yCoord y-axis co-ordinate
     */
    void setInterpOut (const double xCoord, const double yCoord);
    /**
     sets which function the is used when getting output
     */
    void setOutputFunction(OutputMethod outType);
    //==========================================================================
    /**
     Get the output of the plate using the outputFunction pointer.
     
     @return returns value of which ever function outputFunction points to
     */
    double getOutput();
    
    /**
     Process an audio signal using the scheme as a plate reverb unit
     
     @param force the audio signal
     @return output from the plate at the desired point
     */
    double reverb (double force);
    //==========================================================================
    /**
     Update the time state of the scheme
     */
    void updateScheme();
    /**
     add force to the relevant section of the plate
     
     @param force force in Newtons
     */
    void addForce (double force);
    /**
     UNDER CONSTRUCTION: will add a strike force to the plate.
     */
    void addStrike();
    //==========================================================================
    /**
     Print plate parameter information
     */
    void printInfo();
    /**
     Print Internal Coefficients
     */
    void printCoefs();
    //==========================================================================
    /**
     Signum Function
     
     @param input value to operate on
     @return returns 0 if d is `<` 0 or 1 if d is `>` 0
     */
    int sgn (double input);
    
private: /// Methods
    //==========================================================================
    /**
     setup the plate with a given sample rate and boundary condition type: still under construction
     
     @param sampRate Sample Rate in Hz
     @param plateParams the PlateParameters struct which defines the specifications of the plate model
     */
    void setup (double sampRate, PlateParameters plateParams);
    //==========================================================================
    /**
     get the velocity output from the plate. Rounding to the nearest grid point
     
     @return velocity at specified read-out point
     */
    double getVelocityOutput();
    /**
     Gets the amplitude output from the plate. Rounding to the nearest grid point
     
     @return amplitude at specified read-out point
     */
    double getAmplitudeOutput();
    /**
     Use the internal interpolation method to calculate the output from the plate
     
     @return interpolated value at specified read-out point
     */
    double getInterpOut();
    //==========================================================================
    /**
     Populates the internal interpolation lookup table.
     */
    void setInterpLookTable();
    //==========================================================================
    /**
     Sets the loss coefficients
     
     @param t60 the T60 decay in seconds
     @param tone percentage of high frequency decay relative to T60, between .1 and 1.
     Values outside this range will be capped
     */
    void setLoss (double t60, double tone);
    /**
     Sets the coeffiecients for the scheme
     
     @param bcType Boundary Condition
     @param rho plate material density
     @param H plate thickness
     */
    void setCoefs (BoundaryCondition bcType, double rho, double H);
    /**
     Sets the grid spacing for the scheme
     */
    void setGridSpacing();
    
public: // Variables
    //==========================================================================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Allocate Memory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** Time State*/
    double * u, * u1, * u2;
    /***/
    double *dummyptr;
    
private: // Variables
    //==========================================================================
    /***/
    const int interpOrder = 4;
    /***/
    const int interpRes = 1000;
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
    double** interpLookTable;
    /** Pointer to one of the sample output function*/
    double (FDPlate::*outputFunction)();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constants
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /***/
    const double pi {3.14159265358979323846};
    /***/
    const int maxgridsize {3000};		// real-time limit 3000 points approx.
    /***/
    const int maxXgrid = 30.;
    //==========================================================================
    /***/
    BoundaryCondition currentBoundCon;
    /***/
    OutputMethod outputType;
    
    //==========================================================================
    // Physical Parameters
    
    /**x-axis plate length (m)*/
    double Lx;
    /**y-axis plate length (m)*/
    double Ly;
    /**readout position as percentage.*/
    double rp[4];
    /**Excitation*/
    double u0, v0, wid, ctr[2]; // excitation displacement and velocity
    
    //==========================================================================
    /**Loss coefficients*/
    double sigma0 ,sigma1;
    //==========================================================================
    /**Scheme Coefficient*/
    double A00, B00, B01, B11, B02, BC1, BC2, C00, C01, d0;
    //==========================================================================
    /**Derived Parameters*/
    int Nx, Ny, ss, li, lo, lol, lor;
    /**Derived Parameters*/
    double kappa, hmin, h, mu, k, SR, readcoordx,readcoordy, readoutpos;
};


#endif /* PlateClasshpp */
