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
#include "../../AudioClasses/AudioOut.h"	// For .wav header, write and playback

int main (int argc, const char *argv[]) {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sampling and Duration
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    double SR, Tf;
    int Nf;
    
    SR = 44.1e3; // Sample Rate
    Tf = 3;	 // duration (seconds)
    Nf = Tf *SR; // duration (samples)
    
    double *out;
    out = new double[Nf];
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
    
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
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Plate Setup
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FDPlate plateTest;
    
    //input is SR, and boundary condition type
    plateTest.setup(SR,1);
    plateTest.setInitialCondition();
    
    //	Print info
    plateTest.printCoefs();
    plateTest.printInfo();
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Process Loop
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This is where the magic happens
    
    for(int n = 0; n < Nf; ++n)
    {
        plateTest.updateScheme();
        out[n] = plateTest.getOutput(1);
    }
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    writeWavMS(out, outputfname, Nf, SR);
    delete out;
    playWavMS(outputfname);
    printf("\nComplete...\n");
    
    std::cout << "SUCCESS" << '\n';
    return 0;
    
}
