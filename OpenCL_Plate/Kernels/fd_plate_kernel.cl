// B=  {B00,B01,B02,B11)
__kernel void fd_plate(__global float* u, __global float* u1, __global float* u2,
                          __constant float* B, __constant float* C,
                          __constant int* Ny, __constant int* ss)
{

    size_t i = get_global_id(0);
    if(i < ss)
    {
    u[i] = B[0] * ( u1[i] ) +
           B[1] * ( u1[i - 1] + u1[i + 1] + u1[i - (*Ny)] + u1[i + (*Ny)] ) +
           B[2] * ( u1[i - 2] + u1[i + 2] + u1[i - (2 * (*Ny))] + u1[i + (2 * (*Ny))] ) +
           B[3] * ( u1[i - 1 - (*Ny)] + u1[i + 1 - (*Ny)] + u1[i + 1 + (*Ny)] + u1[i - 1 + (*Ny)] ) +
           C[0] * ( u2[i] ) +
           C[1] * ( u2[i - 1] + u2[i + 1] + u2[i - (*Ny)] + u2[i + (*Ny)] );
    }
    u[i] = 0;
}
