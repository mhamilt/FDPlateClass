//
//  RunPlateClass.cpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 20/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//    Test the Plate Class is working correctly

#include <iostream>
#include <string>
#include "../../FDTDClasses/FDPlate.cpp"
#include "../../FDTDClasses/FDString.cpp"
#include "../../cpp-cli-audio-tools/src/CliAudioTools.h"
#include "../../CliFileNaming/CliFileName.h"

//==============================================================================

int main (int argc, const char *argv[])
{
    
    //==========================================================================
    // Sampling and Duration
    
    const double sampleRate = 44.1e3;
    const double duration = 4;         // duration (seconds)
    
    //==========================================================================
    // Set Output File Name
    
    const char *fileName = "/Downloads/Plate.wav";
    char *outputfname = CliSetFilename(argv, fileName);
    WavCodec wavCodec;
    
    //==========================================================================
    // Plate Setup
    FDPlate::PlateParameters plateParams;
    plateParams.t60 = 10.3;
    plateParams.thickness = 0.001;
    plateParams.tone = .9;
    plateParams.lengthX = .2;
    plateParams.lengthY = .2;
    plateParams.bcType = FDPlate::BoundaryCondition::simplySupported;
    
    FDPlate plate(sampleRate, plateParams);
    plate.setInitialCondition();
    plate.printCoefs();
    plate.printInfo();
    
    FDString string;
    string.setup(sampleRate, FDString::LossModel::frequencyDepenent, FDString::BoundaryCondition::simplySupported);
    string.addForce();
    string.printCoefs();
    string.printInfo();
    
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
    wavCodec.writeWavMS(out, outputfname, Nf, sampleRate);
    printf("\nComplete...\n");
    std::cout << "SUCCESS" << '\n';
    return EXIT_SUCCESS;
}

