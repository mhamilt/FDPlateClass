//
//  RunPlateClass.cpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 20/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//	Test the Plate Class is working correctly

#include <iostream>
#include <cstring>
#include "../../FDTDClasses/FDPlate.cpp"
#include "../../FDTDClasses/FDString.cpp"
#include "../../AudioIOClasses/AudioOut.h"

int main (int argc, const char *argv[])
{
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sampling and Duration
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const double sampleRate = 44.1e3;
    const double duration = 7;	 // duration (seconds)
    
    //==========================================================================
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set Output File Name
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    char *outputfname;
    int length;
    
    if(!argv[1])
    {
        const char *homedir = getenv("HOME");
        if(!homedir)
        {
            printf("Couldn't find Home Directory. Enter a filename\n");
            return 1;
        }
        printf("NO FILE NAME ENTERED\nUSING DEFAULT FILE DESTINATION\n");
        const char *def_fname = "/Downloads/Plate.wav";
        length = int(strlen(homedir)) + int(strlen(def_fname));
        outputfname = new char[length+1]();
        strncpy(outputfname,homedir, int(strlen(homedir)));
        strcat(outputfname, def_fname);
    }
    else
    {
        length = int(strlen(argv[1]));
        outputfname = new char[length+1]();
        strncpy(outputfname, argv[1], length);
    }
    
    //==========================================================================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Plate Setup
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    FDPlate plateTest;
//    plateTest.setup(sampleRate, FDPlate::BoundaryCondition::clamped);
//    plateTest.setInitialCondition();
//    plateTest.setOutType(FDPlate::OutputMethod::velocity);
//    //    Print info
//    plateTest.printCoefs();
//    plateTest.printInfo();
    
    FDString string;
    string.setup(sampleRate, FDString::LossModel::frequencyDepenent, FDString::BoundaryCondition::simplySupported);
    string.addForce();
    string.printCoefs();
    string.printInfo();
    
    //==========================================================================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Process Loop
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This is where the magic happens
    const double Nf = duration * sampleRate; // duration (samples)
    double *out = new double[Nf];
    
    for(int n = 0; n < Nf; ++n)
    {
        string.updateScheme();
        out[n] = string.getOutput();
    }
    
    //==========================================================================
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    writeWavMS(out, outputfname, Nf, sampleRate);
    playWavMS(outputfname);
    printf("\nComplete...\n");
    std::cout << "SUCCESS" << '\n';
    return 0;
    
}
