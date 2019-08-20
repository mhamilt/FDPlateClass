/*
 * Testing out OpenCL on the FD Plate time step update
 *
 */


#include <iostream>
#include "../OclPlate/OclPlate.hpp"
#include "../../cpp-cli-audio-tools/src/MattsAudioTools.h"
#include "../../CliFileNaming/CliFileName.h"

int main(int argc, const char * argv[])
{
    const char *fileName = "/Downloads/Plate.wav";
    char *outputfname = CliSetFilename(argv, fileName);
    const double sampleRate = (44.1e3);
    const double duration = 2.;         // duration (seconds)
    
    FDPlate::PlateParameters plateParams;
    //    plateParams.poisson  = .5;
    //    plateParams.density = 480.;
    //    plateParams.t60 = 3.3;
    plateParams.thickness = 0.001;
    //    plateParams.tone = .8;
    plateParams.lengthX = 1.2;
    plateParams.lengthY = 1.2;
    //    plateParams.youngs = 11e9;
    //    plateParams.bcType = FDPlate::BoundaryCondition::simplySupported;
    
    OclPlate plate;
    plate.setup(sampleRate,plateParams);
    plate.setInitialCondition();
    plate.printCoefs();
    plate.printInfo();
    plate.setupCl();
    // This is where the magic happens
    const int Nf = duration * sampleRate; // duration (samples)
    float *out = new float[Nf];
    plate.clupdate(out, Nf);
    
    //    WavCodec::normaliseBuffer(out, Nf);
    //        plate.printMap();
    //    AudioPlayerOpenAL ap;
    //    ap.playAudioData(out, Nf, 1, sampleRate, 16);
    
    WavCodec wavCodec;
    wavCodec.writeWavMS(out, outputfname, Nf, sampleRate);
    return 0;
}
