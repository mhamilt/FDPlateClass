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
    void clupdate(float* output, int frames);
    /**
     <#Description#>
     */
    void printClDevice();
    void setupCl();
private:
    int nextSquare(int N);
    
private:
    /// queue that can send work to the GPU in our system.
    dispatch_queue_t queue;
};

#endif /* OclPlate_hpp */
