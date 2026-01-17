#include "mex.h"
#include "matrix.h"
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


void call_T04_s(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]);
void call_RECALC_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]);
void call_IGRF_GSW_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]);
void call_TRANS_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[], 
                 void(*transfun)(double*,double*,double*,double*,double*,double*,int*));
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]);
int mat2int(const mxArray *mat);
int matchecksize(int narray, mwSize ndim, const mwSize *sz, ...);