// B=  {B00,B01,B02,B11,BC1,BC2}
__kernel void fd_plate(__global float* u, __global float* u1, __global float* u2,
                       __constant float* B, __constant float* C,
                       __constant int* Ny,__constant int* Nx, __constant int* ss)
{
    
    size_t i = get_global_id(0);
    if(i < (*ss))
    {
        const int cp = i;
        const int cpm = (i % (*Ny));
        if ( cp > (2*(*Ny)) &&
            cp < ((*Ny)*((*Nx)-2)) &&
            cpm > 1 &&
            cpm < ((*Ny)-2))
        {
            u[cp] = B[0] * u1[cp] +
            B[1] * ( u1[cp - 1] + u1[cp + 1] + u1[cp - (*Ny)] + u1[cp + (*Ny)] ) +
            B[2] * ( u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * (*Ny))] + u1[cp + (2 * (*Ny))] ) +
            B[3] * ( u1[cp - 1 - (*Ny)] + u1[cp + 1 - (*Ny)] + u1[cp + 1 + (*Ny)] + u1[cp - 1 + (*Ny)] ) +
            C[0] * u2[cp] +
            C[1] * ( u2[cp - 1] + u2[cp + 1] + u2[cp - (*Ny)] + u2[cp + (*Ny)] );
        }
        else if ((cp > (*Ny) &&
                  cp < ((*Ny)*((*Nx)-1)) &&
                  cpm > 1 &&
                  cpm < ((*Ny)-2)) ||
                 (cp > (2*(*Ny)) &&
                  cp < ((*Ny)*((*Nx)-2)) &&
                  (cpm == 1 ||
                   cpm == (*Nx)-2)))
        {
            u[cp]  = B[4]*u1[cp] +
            B[1] * ( u1[cp - 1] + u1[cp + 1] + u1[cp - (*Ny)] + u1[cp + (*Ny)]) +
            B[2] * ( u1[cp - 2] + u1[cp + 2] + u1[abs((cp - (2 * (*Ny))) % (*ss))] + u1[abs((cp + (2 * (*Ny)))%(*ss))]) +
            B[3] * ( u1[cp - 1 - (*Ny)] + u1[cp + 1 - (*Ny)] + u1[cp + 1 + (*Ny)] + u1[cp - 1 + (*Ny)] ) +
            C[0] * ( u2[cp] ) +
            C[1] * ( u2[cp - 1] + u2[cp + 1] + u2[cp + (*Ny)] + u2[cp - (*Ny)] );
        }
        else if ((cp == ((*Ny) + 1)) ||
                 (cp == (2*(*Ny) - 2)) ||
                 (cp == (((*Nx)-2)*(*Ny) + 1)) ||
                 (cp == (((*Nx)-1)*(*Ny) - 2)))
        {
            u[cp] = B[5]*u1[cp] +
            B[1]*( u1[cp-1] + u1[cp+1] + u1[cp-(*Ny)] + u1[cp+(*Ny)] ) +
            B[2]*( u1[cp+2] + u1[abs((cp+(2*(*Ny)))%(*ss))] + u1[cp-2] + u1[abs((cp-(2*(*Ny)))%(*ss))])   +
            B[3]*( u1[cp-1-(*Ny)] + u1[cp+1-(*Ny)] +u1[cp+1+(*Ny)] + u1[cp-1+(*Ny)] ) +
            C[0]*u2[cp] +
            C[1]*( u2[cp-1] + u2[cp+1] + u2[cp-(*Ny)] + u2[cp+(*Ny)] );
        }
    }
    else
    {
        u[i] = 0.f;
    }
}
