/*
 * Testing out OpenCL on the FD Plate time step update
 *
 */


#include <iostream>
#include "../OclPlate/OclPlate.hpp"
#include "../../cpp-cli-audio-tools/src/MattsAudioTools.h"


int main(int argc, const char * argv[])
{
    const double sampleRate = 44.1e3;
    const double duration = 4;         // duration (seconds)
    
    OclPlate plate;
    plate.setInitialCondition();
    plate.printCoefs();
    plate.printInfo();
    
    // This is where the magic happens
    const double Nf = duration * sampleRate; // duration (samples)
    float *out = new float[Nf];
    plate.clupdate();
//    for(int n = 0; n < Nf; ++n)
//    {
//        out[n] = plate.clupdate();
//    }
//    plate.printMap();
//    AudioPlayerOpenAL ap;
//    ap.playAudioData(out, Nf, 1, sampleRate, 16);
    return 0;
}
