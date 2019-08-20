//
//  PlateClass.cpp
//  FDTDCPlate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//  Class file for a FDTD Plate
//
//

#include "FDPlate.hpp"

//==============================================================================
FDPlate::FDPlate() : FDPlate(44100, PlateParameters())
{
}

FDPlate::FDPlate(float sampleRate) : FDPlate(sampleRate, PlateParameters())
{
}

FDPlate::FDPlate (float sampRate, PlateParameters plateParams)
{
    setup(sampRate, plateParams);
}

//==============================================================================

void FDPlate::setup(float sampRate, FDPlate::PlateParameters plateParams)
{
    // I/O Parameters
    rp[0] = .5; rp[1]=.4; rp[2] = .5; rp[3]= .5; // readout position as percentage.
    
    //Excitation
    ctr[0] = .51; ctr[1] = .49;    // centre point of excitation as percentage
    wid = .25;                     // width (m)
    u0 = 0; v0 = 1;                // excitation displacement and velocity
    
    const float E   = plateParams.youngs;
    const float H   = plateParams.thickness;
    const float nu  = plateParams.poisson;
    const float rho = plateParams.density;
    const float D = (E * pow (H, 3) ) / (12 * (1 - pow (nu,2)));
    kappa = sqrt (D / (rho *  H) );
    
    SR = sampRate;               // internal class sampling rate
    k = 1 / SR;                    // time step
    
    Lx = plateParams.lengthX;
    Ly = plateParams.lengthY;
    
    setLoss (plateParams.t60,plateParams.tone);
    setGridSpacing();
    setCoefs (plateParams.bcType, rho, H);
    
    //Set Input and Output Indeces
    li = (Ny * (ctr[1] * Nx)) + (ctr[0] * Ny);
    lo = (Ny * (rp[1] * Nx)) +  (rp[0] * Ny);
    
    //    Update flags
    currentBoundCon  = plateParams.bcType;
    
    u = new float[ss];
    u1 = new float[ss];
    u2 = new float[ss];
    std::fill (u,  u + ss,  0);
    std::fill (u1, u1 + ss, 0);
    std::fill (u2, u2 + ss, 0);
    
    setInterpLookTable();
    setInterpOut (rp[0], rp[1]);
    setOutputFunction(OutputMethod::velocity);
}
//==============================================================================

FDPlate::~FDPlate()
{
    delete [] u;
    delete [] u1;
    delete [] u2;
    delete [] xInterpIndeces;
    delete [] yInterpIndeces;
    
    for (int i = 0; i < interpOrder; ++i)
        delete [] interpLookTable[i];
    
    delete [] interpLookTable;
}

//==============================================================================

void FDPlate::setLoss (float t60, float tone)
{
    tone = range(tone, 0.01, 0.99);
    // tone = std::clamp(tone, 0.01, 0.99);
    const float lowFrequencyBand  = 100;
    const float highFrequencyBand = 1000;
    const float highT60 = t60 * tone;
    const float z1 = 2 * kappa * (2 * pi * lowFrequencyBand) / (2 * pow (kappa,2));
    const float z2 = 2 * kappa * (2 * pi * highFrequencyBand) / (2 * pow (kappa,2));
    sigma0 = 6 * log (10) * (-z2 / t60 + z1 / highT60) / (z1 - z2);
    sigma1 = 6 * log (10) * (1 / t60 - 1 / highT60) / (z1 - z2);
}

//==============================================================================

void FDPlate::setGridSpacing()
{
    // stability condition
    const float stabilityCondition = (sqrt (4 * k * (sigma1 + sqrt (pow (sigma1,2) + pow (kappa,2)))));
    
    if (Lx / maxXgrid < stabilityCondition)
    {
        hmin = stabilityCondition;
    }
    else
    {
        hmin =    Lx / maxXgrid;
    }
    Nx = floor (Lx / hmin);        // number of segments x-axis
    Ny = floor (Ly / hmin);        // number of segments y-axis
    
    h = sqrt (Lx * Ly / (Nx * Ny));;    // adjusted grid spacing x/y
    Nx = Nx + 1; Ny = Ny + 1;        // grid point number x and y
    mu = (kappa * k) / pow (h,2);    // scheme parameter
    ss = Nx * Ny;                    // total grid size.
    
}
//==============================================================================

void FDPlate::setCoefs (BoundaryCondition bcType, float rho, float H)
{
    //update flag
    currentBoundCon = bcType;
    
    // coefficients are named based on position on the x and y axes.
    A00 = 1 / (1 + k * sigma0); // Central Loss Coeffient (INVERTED)
    
    //// Current time step (B) coeffients
    // There are six unique coefficients for B coefs
    B00 = (-pow (mu,2) * 20 + (2 * sigma1 * k / pow (h,2)) * -4 + 2) * A00; // center
    B01 = (-pow (mu,2) * -8 + (2 * sigma1 * k / pow (h,2))) * A00;        // 1-off
    B11 = (-pow (mu,2) * 2) * A00;                                  // diag
    B02 = (-pow (mu,2) * 1) * A00;                                  // 2-off
    
    switch (bcType)
    {
        case BoundaryCondition::clamped:
        {
            BC1 = (-pow (mu,2) * 21 + (2 * sigma1 * k / pow (h,2)) * -4 + 2) * A00;    // Side
            BC2 = (-pow (mu,2) * 22 + (2 * sigma1 * k / pow (h,2)) * -4 + 2) * A00;    // Corner
            break;
        }
        case BoundaryCondition::simplySupported:
        default:
        {
            BC1 = (-pow (mu,2) * 19 + (2 * sigma1 * k / pow (h,2)) * -4 + 2) * A00;    // Side
            BC2 = (-pow (mu,2) * 18 + (2 * sigma1 * k / pow (h,2)) * -4 + 2) * A00;    // Corner
            break;
        }
    }
    
    // Previous time step (C) coeffients
    C00 = (-(2 * sigma1 * k / pow (h,2)) * -4 - (1 - sigma0 * k))  * A00;
    C01 = -(2 * sigma1 * k / pow (h,2))  * A00;
    
    //input force coefficient
    d0 = pow (k,2) / (rho * H * pow (h,2)) * (1 / (1 + k * sigma0)) * A00;
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
    printf("Actual Size     : %.1f m2 \n", Nx * h * Ny * h);
    //    printf("Thickness (mm)  : %.0f mm \n", H*1e3);
    printf("Grid X-Ax       : %d \n", Nx);
    printf("Grid Y-Ax       : %d \n", Ny);
    printf("Total Ps        : %d \n", ss);
    printf("Incell          : %d\n", li);
    printf("Outcell         : %d\n", lo);
    printf("TimeStep        : %.2e\n", k);
    printf("SampRate        : %.2e\n", SR);
    //    printf("Youngs          : %.2e\n", E);
    printf("Sigma 0         : %f\n", sigma0);
    printf("Sigma 1         : %f\n", sigma1);
}

void FDPlate::printMap()
{
    for (int yi = 0; yi < Ny; ++yi)
    {
        for (int xi = 0; xi < Nx; ++xi)
        {
            const int cp = yi + (xi * Ny);
            
            if ( cp > (2*Ny) &&
                cp < (Ny*(Nx-2)) &&
                (cp % Ny) > 1 &&
                (cp % Ny) < (Ny-2))
            {
                printf("%.2e\t",u[cp]);
            }
            else if ((cp > Ny &&
                      cp < (Ny*(Nx-1)) &&
                      (cp % Ny) > 1 &&
                      (cp % Ny) < (Ny-2)) ||
                     (cp > (2*Ny) &&
                      cp < (Ny*(Nx-2)) &&
                      ((cp % Ny) == 1 ||
                       (cp % Ny) == Nx-2)
                      ) )
            {
                printf("%s\t", "AAAAAAAA");
            }
            else if ((cp == (Ny + 1)) ||
                     (cp == (2*Ny -2)) ||
                     (cp == ((Nx-2)*Ny + 1)) ||
                     (cp == ((Nx-1)*Ny - 2)))
            {
                printf("%s\t", "BBBBBBBB");
            }
            else
            {
                printf("0\t\t\t");
            }
            
            
            
        }
        std::cout << '\n';
    }
}

//==============================================================================

void FDPlate::updateScheme()
{
    //    updateCenter();
    updateRefactor();
//    updateSides();
//    updateCorners();
    dummyptr = u2; u2 = u1; u1 = u; u = dummyptr; // swap pointers
}

void FDPlate::originalUpdate()
{
    for(int xi = 2; xi < (Nx-2); ++xi)
    {
        for(int yi = 2; yi < (Ny-2); ++yi)
        {
            const int cp = yi+(xi * Ny);
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
    for(int xi = 2; xi < Nx-2; ++xi)
    {
        //North
        {
            const int cp = 1+(xi * Ny);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp+1-Ny] + u1[cp+1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
        }
        {
            //South
            const int cp = Ny-2 +(xi * Ny);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp-Ny] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp-1-Ny] + u1[cp-1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp-Ny] + u2[cp+Ny] );
        }
    }
    
    // Y-Axis
    for(int yi = 2; yi < Ny-2; ++yi)
    {
        //West
        {
            const int cp = yi+Ny;
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp+1] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp+1+Ny] + u1[cp-1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp+1] + u2[cp+Ny] );
        }
        
        //East
        {
            const int cp = yi + Ny*(Nx-2);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] +u1[cp-(2*Ny)] ) +
            B11*( u1[cp-1-Ny] + u1[cp+1-Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] );
        }
    }
    
    // Corner Boundaries
    {
        const int cp = Ny+1;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp+2] + u1[cp+(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = 2*(Ny-1);
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp-2] + u1[cp+(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = Ny*(Nx-2)+1;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp+2] + u1[cp-(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = Ny*(Nx-1) - 2;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp-2] + u1[cp-(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    
    // swap pointers
    dummyptr = u2; u2 = u1; u1 = u; u = dummyptr;
}
//==============================================================================

float FDPlate::getOutput()
{
    return (this->*outputFunction)();
}

float FDPlate::getVelocityOutput()
{
    return (u1[lo] - u2[lo]) * SR; // Velocity out
}

float FDPlate::getAmplitudeOutput()
{
    return u[lo];
}

float FDPlate::getInterpOut()
{
    const int order = interpOrder;
    float interpOut = 0;
    
    for (int xi = 0; xi < order; ++xi)
    {
        for (int yi = 0; yi < order; ++yi)
        {
            // out of bound test to mark point as zero
            if (yInterpIndeces[yi] < 0 || yInterpIndeces[yi] > Ny || xInterpIndeces[xi] < 0 || xInterpIndeces[xi] > Nx)
            {
                // hiResValue += 0;
            }
            else
            {
                int cp = yInterpIndeces[yi] + ((xInterpIndeces[xi]) * (Ny - 1));
                interpOut += (interpLookTable[xi][xAlphaIndex] * interpLookTable[yi][yAlphaIndex]) *
                (u1[cp] - u2[cp]) * SR;
            }
        }
    }
    return interpOut;
}

//==============================================================================
// Set the readout position on the plate

void FDPlate::setOutputPosition (float xcoord, float ycoord)
{
    int readoutpos = floor((Ny * (xcoord * Nx)) +  (ycoord * Ny));
    
    // Ensure that read out position is a valid grid point
    if (readoutpos % Ny == 0)
    {
        ++readoutpos;
    }
    else if (readoutpos % (Ny) == (Ny - 1))
    {
        --readoutpos;
    }
    if (readoutpos - Ny < 0)
    {
        readoutpos += Ny;
    }
    else if (readoutpos - (ss - Ny) > 0)
    {
        readoutpos -= Ny;
    }
    
    lo = readoutpos;
}


//==============================================================================
// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FDPlate::setOutputFunction (OutputMethod outType)
{
    switch (outType)
    {
        case OutputMethod::velocity:
        {
            outputFunction = &FDPlate::getVelocityOutput;
            break;
        }
        case OutputMethod::amplitude:     // fall through to default
        default:
        {
            outputFunction = &FDPlate::getAmplitudeOutput;
            break;
        }
    }
}

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
    for (int xi = 1; xi < Nx - 1; ++xi)
    {
        const float X = xi * h;
        
        for (int yi = 1; yi < Ny - 1; ++yi)
        {
            const int cp = yi + (xi * Ny);
            const float Y = yi * h;
            const float dist = sqrt (pow (X - (ctr[0] * Lx),2) + pow (Y - (ctr[1] * Ly),2));
            const float ind = sgn((wid * 0.5) - dist);            // displacement (logical)
            const float rc = .5 * ind * (1 + cos (2 * pi * dist / wid)); // displacement
            u2[cp] = u0 * rc;
            u1[cp] = v0 * k * rc;
        }
    }
}
//==============================================================================


void FDPlate::addForce (float force)
{
    u1[li] += d0 * force;
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

float FDPlate::reverb (float force)
{
    addForce (force);
    updateScheme();
    return getOutput();
}

//==============================================================================
// Split this so that it is arbitrary of coordinate, so input is xCoord*Nx

void FDPlate::setInterpOut (const float xCoord, const float yCoord)
{
    //    lo = (Ny*(xCoord*Nx)) +  (yCoord*Ny);
    
    const int order = interpOrder;
    for (int i = 0; i < order; ++i)
    {
        xInterpIndeces[i] = (i + 1 - order / 2.) + floor (xCoord * (Nx));
        yInterpIndeces[i] = (i + 1 - order / 2.) + floor (yCoord * (Ny - 1));
    }
    
    interpPointY = (yCoord * (Ny - 1));
    interpPointX = (xCoord * Nx);
    
    //        const int order = interpOrder;
    const int res = interpRes;
    xAlphaIndex = floor((interpPointX - floor (interpPointX)) * res);
    yAlphaIndex = floor((interpPointY - floor (interpPointY)) * res);
}

//==============================================================================

void FDPlate::setInterpLookTable()
{
    const int order = interpOrder;
    const int res = interpRes;
    interpLookTable = new float*[order];
    for (int i = 0; i < order; ++i)
    {
        interpLookTable[i] = new float[res];
        std::fill (interpLookTable[i], interpLookTable[i] + res - 1, 1);
    }
    
    float *polynomialnormaliser = new float [order];
    std::fill (polynomialnormaliser, polynomialnormaliser + order, 1);
    float *alphas = new float [res];
    
    for (int i = 0; i < res; ++i)
    {
        alphas[i] = (i / float (res)) - 0.5;
    }
    
    float *anchors = new float [order];
    
    if ((order % 2) == 0)
    {
        for (int i = 0; i < order; ++i)
        {
            anchors[i] = -(order - 1) * 0.5 + i;
        }
    }
    else
    {
        for (int i = 0; i < order; ++i)
        {
            anchors[i] = (-(order) * 0.5) + i;
        }
    }
    
    for (int q = 0; q < res; ++q) // loop for every value of alpha
    {
        for (int j = 0; j < order; ++j) // loop for sub polynomial
        {
            for (int m = 0; m < order; ++m) //loop for each point in subpoly
            {
                if (m != j)
                {
                    if (q == 0)
                    {
                        polynomialnormaliser[j] = polynomialnormaliser[j] * (anchors[j] - anchors[m]);
                    }
                    interpLookTable[j][q] *= (alphas[q] - anchors[m]);
                }
            }
            interpLookTable[j][q] /= polynomialnormaliser[j];
        }
    }
    delete[] polynomialnormaliser;
    delete[] alphas;
    delete[] anchors;
}

//==============================================================================
int FDPlate::sgn (float d)
{
    if(d <= 0)
        return 0;
    else
        return 1;
}

//==============================================================================

float FDPlate::range(float value, float min, float max)
{
    if (min > max) // check values are not switched
    {
        const float temp = min;
        min = max;
        max = temp;
    }
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    
    return value;
}

//float FDPlate::range(float value, float min, float max)
//{
//    if (min > max) // check values are not switched
//    {
//        const float temp = min;
//        min = max;
//        max = temp;
//    }
//    
//    if (value < min)
//        value = min;
//    if (value > max)
//        value = max;
//    
//    return value;
//}

int FDPlate::range(int value, int min, int max)
{
    if (min > max) // check values are not switched
    {
        const int temp = min;
        min = max;
        max = temp;
    }
    
    if (value < min)
        value = min;
    if (value > max)
        value = max;
    
    return value;
}
//==============================================================================
void FDPlate::updateCenter()
{
    for(int x = 2; x < (Nx - 2); ++x)
    {
        for(int y = 2; y < (Ny - 2); ++y)
        {
            const int cp = y + (x * Ny);
            u[cp] = B00 * u1[cp] +
            B01 * ( u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny] ) +
            B02 * ( u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * Ny)] + u1[cp + (2 * Ny)] ) +
            B11 * ( u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny] ) +
            C00 * u2[cp] +
            C01 * ( u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny] );
        }
    }
}

void FDPlate::updateRefactor()
{
    for(int x = 0; x < Nx; ++x)
    {
        for(int y = 0; y < Ny; ++y)
        {
            const int cp = y + (x * Ny);
            const int cpm = (cp % Ny);
            if ( cp > (2*Ny) &&
                cp < (Ny*(Nx-2)) &&
                cpm > 1 &&
                cpm < (Ny-2))
            {
                u[cp] = B00 * u1[cp] +
                B01 * ( u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny] ) +
                B02 * ( u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * Ny)] + u1[cp + (2 * Ny)] ) +
                B11 * ( u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny] ) +
                C00 * u2[cp] +
                C01 * ( u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny] );
            }
            else if ((cp > Ny &&
                      cp < (Ny*(Nx-1)) &&
                      (cp % Ny) > 1 &&
                      (cp % Ny) < (Ny-2)) ||
                     (cp > (2*Ny) &&
                      cp < (Ny*(Nx-2)) &&
                      ((cp % Ny) == 1 ||
                       (cp % Ny) == Nx-2)))
            {    
                u[cp]  = BC1*u1[cp] +
                B01 * ( u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
                B02 * ( u1[cp - 2] + u1[cp + 2] + u1[abs((cp - (2 * Ny)) % ss)] + u1[abs((cp + (2 * Ny))%ss)]) +
                B11 * ( u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny] ) +
                C00 * ( u2[cp] ) +
                C01 * ( u2[cp - 1] + u2[cp + 1] + u2[cp + Ny] + u2[cp - Ny] );
            }
            else if ((cp == (Ny + 1)) ||
                     (cp == (2*Ny - 2)) ||
                     (cp == ((Nx-2)*Ny + 1)) ||
                     (cp == ((Nx-1)*Ny - 2)))
            {
                u[cp] = BC2*u1[cp] +
                B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
                B02*( u1[cp+2] + u1[abs((cp+(2*Ny))%ss)] + u1[cp-2] + u1[abs((cp-(2*Ny))%ss)])   +
                B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
                C00*u2[cp] +
                C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
            }
            
        }
    }
}

void FDPlate::updateSides()
{
    for(int xi = 2; xi < Nx-2; ++xi)
    {
        //North
        {
            const int cp = 1+(xi * Ny);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp+1-Ny] + u1[cp+1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
        }
        {
            //South
            const int cp = Ny-2 +(xi * Ny);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp-Ny] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp-(2*Ny)] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp-1-Ny] + u1[cp-1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp-Ny] + u2[cp+Ny] );
        }
    }
    
    // Y-Axis
    for(int yi = 2; yi < Ny-2; ++yi)
    {
        //West
        {
            const int cp = yi+Ny;
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp+1] + u1[cp+Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] + u1[cp+(2*Ny)] ) +
            B11*( u1[cp+1+Ny] + u1[cp-1+Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp+1] + u2[cp+Ny] );
        }
        
        //East
        {
            const int cp = yi + Ny*(Nx-2);
            u[cp]  = BC1*u1[cp] +
            B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] ) +
            B02*( u1[cp-2] + u1[cp+2] +u1[cp-(2*Ny)] ) +
            B11*( u1[cp-1-Ny] + u1[cp+1-Ny] ) +
            C00*u2[cp] +
            C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] );
        }
    }
    
    //    float i_coef = {1,1,}
    //
    // u[cp]  = BC1 * u1[cp] +
    //
    //
    //          int a[] = {0,1,1,1};
    //
    // B01 * ( u1[cp - 1]        + u1[cp + 1]            + u1[cp - Ny]           + u1[cp + Ny] * a[(i + 3) % 4])
    // B01 * ( u1[cp - 1]        + u1[cp + 1]            + u1[cp - Ny] * a[(i + 2) % 4]   u1[cp + Ny]  )
    // B01 * ( u1[cp - 1]        + u1[cp + 1] * a[(i + 1) % 4] + u1[cp - Ny]     + u1[cp + Ny]  )
    // B01 * ( u1[cp - 1] * a[i % 4] + u1[cp + 1]            + u1[cp - Ny]   + u1[cp + Ny]  )
    // B01 * ( u1[cp - 1] * a[i % 4] + u1[cp + 1] *a[i % 4]           + u1[cp - Ny] *a[i % 4]   + u1[cp + Ny] *a[i % 4] )
    //    int cp;
    //
    //    if ((cp>2+Ny || cp < 2*(Ny - 1)))
    //    {
    //
    //    }
    //    else if ((cp>2+Ny || cp < 2*(Ny - 1)))
    //    {
    //
    //    }
    //    else if ((cp>2+Ny || cp < 2*(Ny - 1)))
    //    {
    //
    //    }
    //    else if ((cp>2+Ny || cp < 2*(Ny - 1)))
    //    {
    //
    //    }
    
    //    const int lim[] = {Nx,Ny};
    //
    //    for(int axis = 0; axis < 2; ++axis)
    //    {
    //        for(int i = 2; i < lim[axis]; ++i)
    //        {
    //            for(int side = 0; side < 2; ++side)
    //            {
    //                const int cp = sp[axis][side];
    //            }
    //        }
    //    }
}

void FDPlate::updateCorners()
{
    // Corner Boundaries
    //    const int corner[] = {Ny + 1,2 * (Ny - 1),Ny*(Nx - 2) + 1,Ny*(Nx - 1) - 2};
    //    for (int i = 0; i < 4; ++i)
    //    {
    //        const int cp = corner[i];
    //        u[cp] = BC2 * u1[cp] +
    //        B01 * ( u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny] ) +
    //        B02 * ( u1[cp + (2 * ((i % 1) ? 1 : -1))] + u1[cp + (2 * Ny)] ) +
    //        B11 * ( u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny] ) +
    //        C00 * u2[cp] +
    //        C01 * ( u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny] );
    //    }
    
    {
        const int cp = Ny+1;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp+2] + u1[cp+(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = 2*(Ny-1);
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp-2] + u1[cp+(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = Ny*(Nx-2)+1;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp+2] + u1[cp-(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
    {
        const int cp = Ny*(Nx-1) - 2;
        u[cp] = BC2*u1[cp] +
        B01*( u1[cp-1] + u1[cp+1] + u1[cp-Ny] + u1[cp+Ny] ) +
        B02*( u1[cp-2] + u1[cp-(2*Ny)] ) +
        B11*( u1[cp-1-Ny] + u1[cp+1-Ny] +u1[cp+1+Ny] + u1[cp-1+Ny] ) +
        C00*u2[cp] +
        C01*( u2[cp-1] + u2[cp+1] + u2[cp-Ny] + u2[cp+Ny] );
    }
}
//==============================================================================
// Fun with overloading operators


void FDPlate::operator++(int)
{
    updateScheme();
    //    originalUpdate();
}
void FDPlate::operator<<(float input)
{
    addForce (input);
}
void FDPlate::operator>>(float& output)
{
    output = getVelocityOutput();
}
