//
//  OclPlate.hpp
//  OpenCL_Plate
//
//  Created by mhamilt7 on 20/08/2019.
//  Copyright © 2019 mhamilt7. All rights reserved.
//

#ifndef OclPlate_hpp
#define OclPlate_hpp

#include "../../FDTDClasses/FDPlate.hpp"
#include <OpenCL/opencl.h>

class OclPlate : public FDPlate
{
public:
    OclPlate();
    ~OclPlate();
    float clupdate();
    /**
     <#Description#>
     */
    void printClDevice();
private:
    int nextSquare(int N);
    void setupCl();
    
private:
    unsigned int clItems;
    /// cl memory state
    void* cl_u;
    /// cl memory state
    void* cl_u1;
    /// cl memory state
    void* cl_u2;
    /// cl fd coefs
    void* cl_B;
    /// cl fd coefs
    void* cl_C;
    /// cl number of grid point in y axis
    void* cl_Ny;
    void* cl_ss;
    /// queue that can send work to the GPU in our system.
    dispatch_queue_t queue;
    cl_float* u_local;
    cl_float* u1_local;
    cl_float* u2_local;
    //--------------------------------------------------------------------------
    // Coeffs
    cl_float* B_local;
    cl_float* C_local;
    cl_int* Ny_local;
    cl_int* ss_local;
    
};

#endif /* OclPlate_hpp */
