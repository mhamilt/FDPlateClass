//
//  RunPlateClass.cpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 20/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//	Test the Plate Class is working correctly

#include <iostream>
#include <string>
#include "../../FDTDClasses/FDPlate.cpp"
#include "../../FDTDClasses/FDString.cpp"
#include "../../AudioIOClasses/AudioOut.h"
#include "../../CliFileNaming/CliFileName.h"

//==============================================================================

int main (int argc, const char *argv[])
{
    
    //==========================================================================
    // Sampling and Duration
    
    const double sampleRate = 44.1e3;
    const double duration = 4;	 // duration (seconds)
    
    //==========================================================================
    // Set Output File Name
    
    const char *fileName = "/Downloads/Plate.wav";
    char *outputfname = CliSetFilename(argv, fileName);
    
    //==========================================================================
    // Plate Setup
    
    FDPlate::PlateParameters plateParams;
    plateParams.thickness = 0.003;
    plateParams.tone = .5;
    plateParams.lengthX = .6;
    plateParams.lengthY = .4;
    FDPlate plate(sampleRate, plateParams);
    
//    plate.setup(sampleRate, plateParams);
    
    plate.setInitialCondition();
//    plate.setOutputFunction(FDPlate::OutputMethod::amplitude);
    
    plate.printCoefs();
    plate.printInfo();
    
//    FDString string;
//    string.setup(sampleRate, FDString::LossModel::frequencyDepenent, FDString::BoundaryCondition::simplySupported);
//    string.addForce();
//    string.printCoefs();
//    string.printInfo();
    
    //==========================================================================
    // Process Loop

    // This is where the magic happens
    const double Nf = duration * sampleRate; // duration (samples)
    double *out = new double[Nf];
    
    for(int n = 0; n < Nf; ++n)
    {
        plate.updateScheme();
        out[n] = plate.getOutput();
    }
    
    //==========================================================================
    // Output
    
    writeWavMS(out, outputfname, Nf, sampleRate);
    playWavMS(outputfname);
    printf("\nComplete...\n");
    std::cout << "SUCCESS" << '\n';
    return EXIT_SUCCESS;
}
