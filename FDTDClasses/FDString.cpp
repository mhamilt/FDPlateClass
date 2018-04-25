//
//  StringClass.cpp
//  FDTD_C_String
//
//  Created by admin on 07/11/2017.
//  Copyright © 2017 admin. All rights reserved.
//

#include "FDString.hpp"

//==============================================================================
FDString::FDString()
{
    outputType = OutputMethod::velocity;
    inputType = InputMethod::strike;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // physical parameters
    //         E      nu    Rho
    // Steel : 2e11   0.30    8050
    // Alum  : 7e10   0.35    2700
    // Lead  : 1.6e10 0.44    11340
    // Wood  : 1e10   0.40    480
    
    gauge = 10;     // string gauge
    f0 = 100;       // frequency of note in Hz (see function at EOF)
    E = 2e11;       // Young's modulus (Pa) (GPa = 1e9 Pa) of steel
    nu = .3;        // Poisson Ratios (< .5)
    rho = 8050;     // density (kg/m^3) of steel
    r = (gauge * 2.54e-5)*.5;    // string radius (m)
    L = 1.3;                      // length (m)
    loss[0] = 100; loss[1] = 10;
    loss[2] = 1000; loss[3] = 9;    // loss [freq.(Hz), T60;...]
    T = pow(((2*f0*r)*L),2)*pi*rho; // Tension in Newtons
    
    // // I/O
    xi = 0.8;        // coordinate of excitation (normalised, 0-1)
    xo = 0.23;        // coordinate of readout (normalised, 0-1)
    
    // // Excitation
    famp = 1;        // peak amplitude of excitation (N)
    dur = 0.003;    // duration of excitation (s)
    exc_st = 0;        // start time of excitation (s)
    u0 = 0; v0 = 1;   // initial conditions
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // String Derived Parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // derived parameters
    A = pi*pow(r,2);                       // string cross-sectional area
    I = 0.25*pi*pow(r,4);                  // string moment of intertia
    c = sqrt(T/(rho*A));        // wave speed
    K = sqrt(E*I/(rho*A));   // stiffness constant
    sampleRate = 44.1e3;
    k = 1/sampleRate;
    
}

//==============================================================================
FDString::~FDString()
{
    
}

//==============================================================================

void FDString::setup (double sampRate, LossModel lossType, BoundaryCondition bcType)
{
    sampleRate = sampRate;
    setLoss (lossType);
    setGrid();
    setMemory();
    setCoefs (bcType);
    setInputPosition (xi);
    setOutputPosition (xo);
//    setInitialCondition();
    setupFlag = true;
}

//==============================================================================

void FDString::setLoss (LossModel lossType)
{
    switch (lossType)
    {
        case LossModel::lossless:
            sigma0 = 0;
            sigma1 = 0;
            break;
        case LossModel::simple:
            sigma0 = 6*log(10)/loss[1];
            sigma1 = 0;
            break;
        case LossModel::frequencyDepenent:
        default:
            z1 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*(2*pi*pow(loss[0],2))))/(2*pow(K,2));
            z2 = (-pow(c,2) + sqrt(pow(c,4) + 4*pow(K,2)*(2*pi*pow(loss[2],2))))/(2*pow(K,2));
            sigma0 = 6*log(10)*(-z2/loss[1] + z1/loss[3])/(z1-z2);
            sigma1 = 6*log(10)*(1/loss[1] - 1/loss[3])/(z1-z2);
    }
}

//==============================================================================

void FDString::setCoefs (BoundaryCondition bcType)
{
    A00 = 1/(1+k*sigma0);
    B00 = (-pow(lambda,2)*2 - pow(mu,2)*6 + 2 - (2*sigma1*k/pow(h,2))*2) * A00; // centre
    B01 = (pow(lambda,2) + pow(mu,2)*4 + (2*sigma1*k/pow(h,2)) ) * A00;    // 1-off
    B02 =  (-pow(mu,2)) * A00;    // 2-off
    C00 = - ( (1-sigma0*k) + ((2*sigma1*k/pow(h,2))*-2) )  * A00;
    C01 = -((2*sigma1*k/pow(h,2)))  * A00;
    switch (bcType)
    {
        case BoundaryCondition::simplySupported:
            BC1 = (pow(lambda,2)*-2 - pow(mu,2)*7 + 2 + (2*sigma1*k/pow(h,2))*-2) * A00;
            break;
        case BoundaryCondition::clamped:
        default:
            BC1 = (pow(lambda,2)*-2 - pow(mu,2)*5 + 2 + (2*sigma1*k/pow(h,2))*-2) * A00;
            break;
    }
}

//==============================================================================

void FDString::setGrid()
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // String Derived Parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // derived parameters
    
    A = pi*pow(r,2);            // string cross-sectional area
    I = 0.25*pi*pow(r,4);       // string moment of intertia
    c = sqrt(T/(rho*A));        // wave speed
    K = sqrt(E*I/(rho*A));      // stiffness constant
    k = 1/sampleRate;           // Time Step
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Grid
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //    double hmin = sqrt(0.5* (pow(c,2)*pow(k,2)+sqrt(pow(c,4)*pow(k,4)+16*pow(K,2)*pow(k,2))) );
    
    //// spacing for frequency dependant loss
    hmin = sqrt(pow(c*k,2)+4*sigma1*k + sqrt( pow(pow(c*k,2)+4*sigma1*k,2) +16*pow(K*k,2)));
    
    N = floor(L/hmin);      // number of segments (N+1 is number of grid points)
    h = L/N;                // adjusted grid spacing
    lambda = c*k/h;         // Courant number
    mu = K*k/pow(h,2);      // numerical stiffness constant
    N = N+1;                // change N to be number of grid points
}


//==============================================================================

void FDString::printCoefs()
{
    printf("--- Coefficient Info --- \n\n");
    printf("Loss A    : %.7fm \n", A00);
    printf("Centre B  : %.6fm \n", B00);
    printf("1-Grid B  : %.6fm \n", B01);
    printf("2-Grid B  : %.6fm \n", B02);
    printf("Centre C  : %.6fm \n", C00);
    printf("1-Grid C  : %.6fm \n", C01);
    printf("BoundCon  : %.6fm \n", BC1);
    printf("Sigma 0   : %f\n", sigma0);
    printf("Sigma 1   : %f\n", sigma1);
    printf("z1        : %f\n", z1);
    printf("z2        : %f\n\n", z2);
}

void FDString::printInfo()
{
    printf("\n--- Scheme Info --- \n\n");
    printf("Size        : %.1fm \n", N*h);
    printf("Gridmin     : %.7f \n", hmin);
    printf("GridSpace   : %.7f \n", h);
    printf("Grid Num    : %d \n", N);
    printf("In_cell     : %d\n", li);
    printf("Out_cell    : %d\n", lo);
    printf("Youngs      : %.2e\n\n", E);
}

//==============================================================================

void FDString::updateScheme()
{
    // Internal Gride Points
    for(int i = 2; i < N-2; i++)
    {
        u[i] = B00*u1[i] +
        B01*( u1[i-1] + u1[i+1]) +
        B02*( u1[i-2] + u1[i+2]) +
        C00*u2[i] +
        C01*( u2[i-1] + u2[i+1]);
    }
    
    // Boundaries
    int i = 1;
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
    
    updateForce();

    // Pointer Swap
    dummy_ptr = u2; u2 = u1; u1 = u; u = dummy_ptr;
}

//==============================================================================

double FDString::getOutput ()
{
    double sampleOut;
    switch (outputType)
    {
        case OutputMethod::velocity:
        {
            sampleOut = (u1[lo]- u2[lo])*sampleRate; // Velocity out
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

void FDString::setOutputPosition (float dist)
{
    
    //TODO: CHECK IF ON A VALID GRID POINT
    if((dist*N)-1 < 1 || (dist*N)-1 > N-2)
    {
        //Do something to ensure it is not on a zero point
    }
    
    lo = (dist*N)-1;
}

void FDString::setInputPosition (float dist)
{
    
    //TODO: CHECK IF ON A VALID GRID POINT
    if((dist*N)-1 < 1 || (dist*N)-1 > N-2)
    {
        //Do something to ensure it is not on a zero point
    }
    
    li = (dist*N)-1;
}

//==============================================================================

void FDString::setOutType (OutputMethod outType)
{
    outputType = outType;
}

//==============================================================================

void FDString::setInitialCondition()
{
    u2[li] = u0;
    u1[li] = u0 + (v0*k);
}

//==============================================================================

void FDString::setMemory()
{
    u =  new double[maxGridSize];
    std::fill(u, u+maxGridSize, 0);
    u1 = new double[maxGridSize];
    std::fill(u1, u+maxGridSize, 0);
    u2 = new double[maxGridSize];
    std::fill(u2, u+maxGridSize, 0);
    force = new double[maxExcDur];
    std::fill(force, force+maxExcDur, 0);
}

//==============================================================================

void FDString::addForce()
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const double d0 = pow(k,2)/(h*rho*A*(1+k*sigma0))*A00;    // force coeffcient
    switch (inputType)
    {
        case InputMethod::strike:
        {
            const double forcePeriod = 2*pi/double(maxExcDur);
            for (int i = 0; i < maxExcDur; ++i)
            {
                const int cp = (i+fi) % maxExcDur;
                force[cp] += ((1 - cos(i*forcePeriod))*.5)*d0;
            }
            break;
        }
        case InputMethod::pluck:
        default:
        {
            const double forcePeriod = pi/double(maxExcDur);
            for (int i = 0; i < maxExcDur; ++i)
            {
                const int cp = (i+fi) % maxExcDur;
                force[cp] += ((1 - cos(i*forcePeriod))*.5)*d0;
            }
            break;
        }
    }
}

//==============================================================================

void FDString::updateForce()
{
    u[li] += force[fi];
    force[fi] = 0;
//    ++fi;
//    fi %= maxExcDur;
    fi = (fi + 1 == maxExcDur ? 0: fi + 1);
}
