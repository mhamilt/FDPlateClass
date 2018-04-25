//
//  PlateClass.cpp
//  FDTDCPlate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//  Class file for a FDTD Plate

#include "FDPlate.hpp"

//==============================================================================
FDPlate::FDPlate()
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Flags
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    currentBoundCon = BoundaryCondition::simplySupported;   // set boundary condition (s);
    outputType = OutputMethod::velocity;        // set output type 0: displacement, 1: velocity
    setupFlag = false;    // flag if Setup() has been run
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Physical Parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //         E      nu    Rho
    // Steel : 2e11   0.30    8050
    // Alum  : 7e10   0.35    2700
    // Lead  : 1.6e10 0.44    11340
    // Wood  : 1e10   0.40    480
    
    // // wood
    E = 11e9;                        // Young's modulus
    rho = 480;                       // density (kg/m^3)
    nu = .5;                         // Poisson Ratios (< .5)
    H = .006;                        // thickness (m)
    Lx = 1;                          // x-axis plate length (m)
    Ly = 1;                          // y-axis plate length (m)
    loss[0] = 100; loss[1] = 1;
    loss[2] = 1000; loss[3] = .2;    // loss [freq.(Hz), T60;...]
    
    // I/O Parameters
    rp[0] = .45; rp[1]=.65; rp[2] = .45; rp[3]= .15; // readout position as percentage.
    
    //Excitation
    ctr[0] = .35; ctr[1] = .45;    // centre point of excitation as percentage
    wid = .25;                    // width (m)
    u0 = 0; v0 = 1;                // excitation displacement and velocity
    
};

//==============================================================================

void FDPlate::setup (double sampRate, BoundaryCondition bcType)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Motion Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    D = (E*(pow (H, 3)))/(12*(1-pow (nu,2)));
    kappa = sqrt (D / (rho*  H) );
    
    SR = sampRate;                // internal class sampling rate
    k = 1/SR;                    // time step
    
    setLoss (8,.75);
    setGrid();
    setCoefs (bcType);
    
    //Set Input and Output Indeces
    li = (Ny*(ctr[1]*Nx)) + (ctr[0]*Ny);
    lo = (Ny*(rp[1]*Nx)) +  (rp[0]*Ny);
    
    //    Update flags
    setupFlag = true;
    currentBoundCon  = bcType;
    
    u = new double[ss];
    u1 = new double[ss];
    u2 = new double[ss];
    std::fill (u,  u+ss,  0);
    std::fill (u1, u1+ss, 0);
    std::fill (u2, u2+ss, 0);
}
//==============================================================================

void FDPlate::setLoss (double lowT60, double highT60Percent)
{
    
    loss[1] = lowT60; loss[3] = lowT60*highT60Percent;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Loss coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    z1 = 2*kappa*(2*pi*loss[0])/(2*pow (kappa,2));
    z2 = 2*kappa*(2*pi*loss[2])/(2*pow (kappa,2));
    sigma0 = 6*log (10)*(-z2/loss[1] + z1/loss[3])/(z1-z2);
    sigma1 = 6*log (10)*(1/loss[1] - 1/loss[3])/(z1-z2);
    
    if (setupFlag)
    {
        setGrid();
        setCoefs (currentBoundCon);
    }
}

//==============================================================================

void FDPlate::setGrid()
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Scheme Spacing
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const int maxXgrid = 30.;
    // stability condition
    const double stabilityCondition = (sqrt (4*k*(sigma1+sqrt (pow (sigma1,2)+pow (kappa,2)))));
    
    if (Lx/maxXgrid < stabilityCondition)
    {
        hmin = stabilityCondition;
    }
    else
    {
        hmin =    Lx/maxXgrid;
    }
    Nx = floor (Lx/hmin);        // number of segments x-axis
    Ny = floor (Ly/hmin);        // number of segments y-axis
    
    h = sqrt (Lx*Ly/(Nx*Ny));;    // adjusted grid spacing x/y
    Nx = Nx+1; Ny = Ny+1;        // grid point number x and y
    mu = (kappa * k)/pow (h,2);    // scheme parameter
    ss = Nx*Ny;                    // total grid size.
    
}
//==============================================================================

void FDPlate::setCoefs (BoundaryCondition bcType)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Scheme Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //update flag
    currentBoundCon = bcType;
    
    // coefficients are named based on position on the x and y axes.
    A00 = 1/(1+k*sigma0); // Central Loss Coeffient (INVERTED)
    
    //// Current time step (B) coeffients
    // There are six unique coefficients for B coefs
    B00 = (-pow (mu,2)*20 + (2*sigma1*k/pow (h,2))*-4 + 2) * A00; // center
    B01 = (-pow (mu,2)*-8 + (2*sigma1*k/pow (h,2))) * A00;          // 1-off
    B11 = (-pow (mu,2)*2) * A00;                                      // diag
    B02 = (-pow (mu,2)*1) * A00;                                      // 2-off
    
    switch (bcType)
    {
        case BoundaryCondition::clamped:
        {
            BC1 = (-pow (mu,2)*21 + (2*sigma1*k/pow (h,2))*-4 + 2) * A00; // Side
            BC2 = (-pow (mu,2)*22 + (2*sigma1*k/pow (h,2))*-4 + 2) * A00; // Corner
            break;
        }
        case BoundaryCondition::simplySupported:
        default:
        {
            BC1 = (-pow (mu,2)*19 + (2*sigma1*k/pow (h,2))*-4 + 2) * A00; // Side
            BC2 = (-pow (mu,2)*18 + (2*sigma1*k/pow (h,2))*-4 + 2) * A00; // Corner
            break;
        }
    }
    
    //// Previous time step (C) coeffients
    C00 = (-(2*sigma1*k/pow (h,2))*-4 - (1-sigma0*k))  * A00;
    C01 = -(2*sigma1*k/pow (h,2))  * A00;
    
    //input force coefficient
    d0 = pow (k,2)/(rho*H*pow (h,2))*(1/(1+k*sigma0))*A00;
}

//==============================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Print Scheme Info
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FDPlate::printCoefs()
{
    printf("--- Coefficient Info --- \n\n");
    printf("Loss A      : %.4fm \n", A00);
    printf("Centre B    : %.4fm \n", B00);
    printf("1-Grid B    : %.4fm \n", B01);
    printf("2-Grid B    : %.4fm \n", B02);
    printf("Diagonal B  : %.4fm \n", B11);
    printf("Centre C    : %.4fm \n", C00);
    printf("1-Grid C    : %.4fm \n", C01);
    printf("Side Bound  : %.4fm \n", BC1);
    printf("Cornr Bound : %.4fm \n\n", BC2);
}

void FDPlate::printInfo()
{
    printf("--- Scheme Info --- \n\n");
    printf("Size            : %.1f m2 \n", Nx*h*Ny*h);
    printf("Thickness (mm)  : %.0f mm \n", H*1e3);
    printf("Grid X-Ax       : %d \n", Nx);
    printf("Grid Y-Ax       : %d \n", Ny);
    printf("Total Ps        : %d \n", ss);
    printf("Incell          : %d\n", li);
    printf("Outcell         : %d\n", lo);
    printf("TimeStep        : %.2e\n", k);
    printf("SampRate        : %.2e\n", SR);
    printf("Youngs          : %.2e\n", E);
    printf("Sigma 0         : %f\n", sigma0);
    printf("Sigma 1         : %f\n", sigma1);
}

//==============================================================================
//==============================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update Plate State
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// TODO: update should not be run before Setup()
// Do a check to ensure

void FDPlate::updateScheme()
{
    int cp;
    // Internal Gride Points
    for (int xi = Nx-4; --xi; )
    {
        for (int yi = Ny-4;--yi; )
        {
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
    
    for (int xi = Nx-4; --xi; )
    {
        //North
        cp = 1+((xi+2) * Ny); // current point
        u[cp]  = BC1*u1[cp] +
        B01*( u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp-2] + u1[cp+2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
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
    
    for (int yi = Ny-4;--yi; )
    {
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
    dummyptr = u2; u2 = u1; u1 = u; u = dummyptr;
}
//==============================================================================
//==============================================================================

// get the output from the plate based on the readout positino and type.
// Will need to implement interpolation, especially f this is meant to be a
// free moving read-out.

double FDPlate::getOutput()
{
    double sampleOut;
    switch (outputType)
    {
        case OutputMethod::velocity:
        {
            sampleOut = (u1[lo]- u2[lo])*SR; // Velocity out
            break;
        }
        case OutputMethod::amplitude: // fall through to default
        default:
        {
            sampleOut = u1[lo]; // Amplitude out
            break;
        }
    }
    return sampleOut;
}

//==============================================================================


void FDPlate::getStereoOutput (OutputMethod outType, double &leftOut, double &rightOut)
{
    switch (outType)
    {
        case OutputMethod::velocity:
        {
            leftOut  = (u1[lol]- u2[lol])*SR;
            rightOut = (u1[lor]- u2[lor])*SR;
            break;
        }
        case OutputMethod::amplitude:
        default:
        {
            leftOut  = u1[lol];
            rightOut = u1[lor];
            break;
        }
    }
}

//==============================================================================
//==============================================================================
// Set the readout position on the plate

void FDPlate::setOutputPosition (double xcoord, double ycoord)
{
    int readoutpos = floor((Ny*(xcoord*Nx)) +  (ycoord*Ny));
    
    // Ensure that read out position is a valid grid point
    if (readoutpos % Ny == 0)
    {
        ++readoutpos;
    }
    else if (readoutpos % (Ny) == (Ny-1))
    {
        --readoutpos;
    }
    if (readoutpos - Ny < 0)
    {
        readoutpos += Ny;
    }
    else if (readoutpos - (ss-Ny) > 0)
    {
        readoutpos -= Ny;
    }
    
    lo = readoutpos;
}

void FDPlate::setStereoOutputPosition (double lxcoord, double lycoord)
{
    // For simplicity, this mirrors left and right output
    // down the middle of the y-axis
    double rxcoord = lxcoord; double rycoord = 1-lycoord;
    
    //TODO: CHECK IF ON A VALID GRID POINT
    if((lxcoord*Nx)-1 < 1 || (lxcoord*Nx)-1 > Nx-2){
        //Do something to ensure it is not on a zero point
    }
    if((rxcoord*Nx)-1 < 1 || (rxcoord*Nx)-1 > Nx-2){
        //Do something to ensure it is not on a zero point
    }
    if((lycoord*Ny)-1 < 1 || (lycoord*Ny)-1 > Nx-2){
        //Do something to ensure it is not on a zero point
    }
    
    if((rycoord*Ny)-1 < 1 || (rycoord*Ny)-1 > Nx-2){
        //Do something to ensure it is not on a zero point
    }
    
    lol = (Ny*(lxcoord*Nx)) +  (lycoord*Ny);
    lor = (Ny*(rxcoord*Nx)) +  (rycoord*Ny);
    
}
//==============================================================================
//==============================================================================
// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FDPlate::setOutType (OutputMethod outType)
{
    outputType = outType;  // set output to velocity amplitude
}

//==============================================================================
//==============================================================================
// Method will set the profile of the input force, for use in JUCE.
// Currently not interpolated, will need to look into that.
// Array values will them selves be multiplied by a rasied cosine/half-cosine

void FDPlate::setInitialCondition()
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Excitation Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // raised cosine in 2D
    for (int xi = 1;xi < Nx-1;++xi)
    {
        const double X = xi*h;
        
        for (int yi = 1;yi < Ny-1;++yi)
        {
            const int cp = yi+(xi * Ny);
            const double Y = yi*h;
            const double dist = sqrt (pow (X-(ctr[0]*Lx),2) + pow (Y-(ctr[1]*Ly),2));
            const double ind = sgn((wid*0.5)-dist);            // displacement (logical)
            const double rc = .5*ind*(1+cos (2*pi*dist/wid)); // displacement
            u2[cp] = u0*rc;
            u1[cp] = v0*k*rc;
        }
    }
}
//==============================================================================


void FDPlate::addForce (double force)
{
    u1[li] += d0* force;
}

//==============================================================================

void FDPlate::addStrike()
{
    //TODO: Function to add a 2D rasied cosine multiplied by a gain rc in
    // displacement, dependent on strike velocity
    // know that strike has begun or ended (insert flag)
    // calculate gain rc based on velocity.
    // prepare for multiple strikes, potentially overlapping.
    
    // Need to work out in advance which indeces will actually be affected by a strike.
}

//==============================================================================

double FDPlate::reverb (double force)
{
    updateScheme();
    addForce (force);
    //    return getOutput (true);
    return getInterpOut();
}

//==============================================================================

double FDPlate::getInterpOut()
{
    const int order = interpOrder;
    double interpOut = 0;
    
    for (int xi = 0; xi < order; ++xi)
    {
        for (int yi = 0; yi < order; ++yi)
        {
            // out of bound test to mark point as zero
            if (yInterpIndeces[yi] < 0 || yInterpIndeces[yi] > Ny || xInterpIndeces[xi] < 0 || xInterpIndeces[xi] > Nx)
            {
                //                hiResValue += 0;
            }
            else
            {
                int cp = yInterpIndeces[yi] + ((xInterpIndeces[xi]) * (Ny-1));
                interpOut += (interpLookTable[xi][xAlphaIndex] * interpLookTable[yi][yAlphaIndex]) *
                (u1[cp]-u2[cp]) * SR;
            }
        }
    }
    return interpOut;
}


//======================================================================
//    LINEAR INTERP (COMMENT OUT)
//======================================================================

//double FDPlate::getInterpOut()
//{
//
//        const int cp = interpZeroIndex;
//
//        const double interpOut = (((1-interpAlphaX)* (1-interpAlphaY) * u1[cp])        +
//                           ((1-interpAlphaX)* (interpAlphaY)   * u1[cp+1])    +
//                           ((interpAlphaX)  * (1-interpAlphaY) * u1[cp+Ny]) +
//                           ((interpAlphaX)  * (interpAlphaY)   * u1[cp+1+Ny])
//
//                           -(((1-interpAlphaX)* (1-interpAlphaY) * u2[cp])    +
//                           ((1-interpAlphaX)* (interpAlphaY)   * u2[cp+1])    +
//                           ((interpAlphaX)  * (1-interpAlphaY) * u2[cp+Ny]) +
//                           ((interpAlphaX)  * (interpAlphaY)   * u2[cp+1+Ny])))*SR;
//
//    return interpOut;
//}
//======================================================================
//    LINEAR INTERP (COMMENT OUT)
//======================================================================


//==============================================================================

//Split this so that it is arbitrary of coordinate, so input is xCoord*Nx

void FDPlate::setInterpOut (const double xCoord, const double yCoord)
{
    //    lo = (Ny*(xCoord*Nx)) +  (yCoord*Ny);
    
    const int order = interpOrder;
    for (int i = 0; i < order;++i)
    {
        xInterpIndeces[i] = (i+1 - order/2.) + floor (xCoord*(Nx));
        yInterpIndeces[i] = (i+1 - order/2.) + floor (yCoord*(Ny-1));
    }
    
    interpPointY = (yCoord*(Ny-1));
    interpPointX = (xCoord*Nx);
    
    //        const int order = interpOrder;
    const int res = interpRes;
    xAlphaIndex = floor((interpPointX-floor (interpPointX))*res);
    yAlphaIndex = floor((interpPointY-floor (interpPointY))*res);
}


//======================================================================
//    LINEAR INTERP (COMMENT OUT)
//======================================================================
//void FDPlate::setInterpOut (const double xCoord, const double yCoord)
//{
//    interpAlphaY = (yCoord*(Ny)) - floor (yCoord*(Ny));
//    interpAlphaX = (xCoord*(Nx)) - floor (xCoord*(Nx));
//    interpZeroIndex = int( floor (yCoord*(Ny-1)) + floor (xCoord*(Nx))* (Ny) );
//
//    if (interpZeroIndex > ((Nx*Ny)-1-Ny))
//    {
//        interpZeroIndex = int( floor (yCoord*(Ny-2)) + (Nx-2)* (Ny) );
//        interpAlphaY = 0;
//        interpAlphaX = 0;
//    }
////    printf("Zero Index: %d\n",interpZeroIndex);
//}
//======================================================================
//    LINEAR INTERP (COMMENT OUT)
//======================================================================

//==============================================================================

double **FDPlate::getInterpLookTable()
{
    const int order = interpOrder;
    const int res = interpRes;
    double **alphaTable = new double*[order];
    
    for (int i = 0;i < order;++i)
    {
        alphaTable[i] = new double[res];
        std::fill (alphaTable[i], alphaTable[i]+res-1, 1);
    }
    
    double *polynomialnormaliser = new double [order];
    std::fill (polynomialnormaliser, polynomialnormaliser+order, 1);
    double *alphas = new double [res];
    
    for (int i = 0; i < res;++i)
    {
        alphas[i] = (i/float (res)) - 0.5;
    }
    
    double *anchors = new double [order];
    
    if ((order % 2)== 0)
    {
        for (int i = 0; i < order;++i)
        {
            anchors[i] = -(order - 1)*0.5 + i;
        }
    }
    else
    {
        for (int i = 0; i < order;++i)
        {
            anchors[i] = (-(order)*0.5) + i;
        }
    }
    
    for (int q = 0; q < res;++q) // loop for every value of alpha
    {
        for (int j = 0; j < order;++j) // loop for sub polynomial
        {
            for (int m = 0; m < order; ++m) //loop for each point in subpoly
            {
                if (m != j)
                {
                    if (q == 0)
                    {
                        polynomialnormaliser[j] = polynomialnormaliser[j]*(anchors[j]-anchors[m]);
                    }
                    alphaTable[j][q] *= (alphas[q]-anchors[m]);
                }
            }
            alphaTable[j][q] /= polynomialnormaliser[j];
        }
    }
    delete[] polynomialnormaliser;
    delete[] alphas;
    delete[] anchors;
    return alphaTable;
}

//==============================================================================
int FDPlate::sgn (double d)
{
    if(d<=0)
        return 0;
    else
        return 1;
}
