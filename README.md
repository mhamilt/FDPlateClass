# FD_Plate_Class

A finite difference Kirchoff Thin-Plate Scheme. Based on the [MATLAB code from my master's thesis](https://github.com/mhamilt/CoupledFDPlateAndString)

## FDTD CLasses

Included are two classes, one for the stiff string and the other for the  Kirchoff Plate.

### FDPlate

Class wrapped method of the [inline code](https://github.com/mhamilt/FDPlate);
The number of grid points is currently limited by the class member

```c
const int maxXgrid = 300;
```

#### Minimal Code

Minimal code to get a sound from the plate:

```cpp
#include <iostream>
#include <string>
#include "../../FDTDClasses/FDPlate.cpp"
#include "../../FDTDClasses/FDString.cpp"
#include "/../../cpp-cli-audio-tools/src/CliAudioTools.h"

int main (int argc, const char *argv[])
{
    const char *fileName = "/absolute/path/to/file.wav";
    const double sampleRate = 44.1e3;
    const double duration = 4;       // duration (seconds)

    FDPlate plate(sampleRate);
    plate.setInitialCondition();

    const double Nf = duration * sampleRate; // duration (samples)
    double *out = new double[Nf];

    for(int n = 0; n < Nf; ++n)
    {
        plate.updateScheme();
        out[n] = plate.getOutput();
    }

    WavCodec wc;
    wc.writeWavMS(out, fileName, Nf, sampleRate);
    return EXIT_SUCCESS;
}
```


### ToDo

- [ ] Methods to couple FDTD systems together.
- [ ] Cleaner method to set system limits on total grid points
- [ ] Move to `std::unique_ptr`
- [ ] complete `addForce (float force)`
- [ ] complete Doxygen documentation

## OpenCL

Out of morbid curiosity, and to force a little code refactoring, included is an OpenCL implementation of the plate. Currently Using the Apple OpenCL API which is now deprecated. There have been the odd machine freeze-ups so it shouldn't be treated as stable.


### ToDo

**Output**
One of the main bottlenecks in the `^kernelBlock` is the `gcl_memcpy` at the end. This was just copying the methodology from Apples own [OpenCL Hello World Example](https://developer.apple.com/library/archive/samplecode/OpenCL_Hello_World_Example/Introduction/Intro.html). To get output efficiently, the next step is to allocate memory on the GPU and read audio data to it with a final read out at the end.

**Work Group Size**

Work Group size appears happiest around `32` rather than the max size of `256`. Defining workgroup size with

```c
gcl_get_kernel_block_workgroup_info(fd_plate_kernel,
                                            CL_KERNEL_WORK_GROUP_SIZE,
                                            sizeof(wgs),
                                            &wgs, NULL);
```

always seems to report back `256`. Stripping out the polling for workgroup size and having a minimal `^kernelBlock` appears to make things run a little quicker. i.e.
```cpp
void (^kernelBlock)() =
    ^{        
        cl_ndrange range =
        {
            1,              // Dimensions
            {0, 0, 0},      // Offset in each dimension.
            {clItems,0, 0}, // global range process.
            {32, 0, 0}      // Workgroup Size
        };

        fd_plate_kernel(&range,
                        (cl_float*)cl_u,
                        (cl_float*)cl_u1,
                        (cl_float*)cl_u2,
                        (cl_float*)cl_B,
                        (cl_float*)cl_C,
                        (cl_int*)cl_Ny,
                        (cl_int*)cl_Nx,
                        (cl_int*)cl_ss);

        void* dummy = cl_u2; cl_u2 = cl_u1; cl_u1 = cl_u; cl_u = dummy; // swap pointers
    };
```
## Cli Audio

Included are a few tools for gerenrating audio from the command line. For the finer details check out the [cpp-cli-audio-tools repo](https://github.com/mhamilt/cpp-cli-audio-tools).

## Notes

Derivations of schemes can be found at the [accompanying GitHub pages site](https://mhamilt.github.io/CoupledFDPlateAndString/)

All code tested with Xcode Version 8.2.1 (8C1002)  
4.2.1 Compatible Apple LLVM 8.0.0 (clang-800.0.42.1)  

[**Remember to add `/usr/local/include` to Header Search Paths**](https://superuser.com/a/898280)
