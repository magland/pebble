#include <mex.h>
#include "jisotonic.h"

 
// [B,MSE]=jisotonic_mex(A,weights)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs<2) {
            mexPrintf("Not enough inputs.\n");
            return;
    }
    if (nlhs<2) {
            mexPrintf("Not enough outputs.\n");
            return;
    }
    int M=mxGetM(prhs[0]);
    int N=mxGetN(prhs[0]);
    if (M!=1) {
            mexPrintf("Input must be a row vector.\n");
            return;
    }
    if (mxGetN(prhs[1])!=N) {
            mexPrintf("Inconsistent dimensions between A and W.\n");
            return;
    }

    double *AA=mxGetPr(prhs[0]);
    double *WW=mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);

    double *BB=mxGetPr(plhs[0]);
    double *MSE=mxGetPr(plhs[1]);

    jisotonic(N,BB,MSE,AA,WW);
}
