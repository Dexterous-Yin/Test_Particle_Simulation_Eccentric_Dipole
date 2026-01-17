#include "igrfmex.h"

// COMMON BLOCKS in geopack_08
extern struct {
    double cb[34];
}GEOPACK1;
extern struct {
    double cb[315];
}GEOPACK2;

// SUBROUTINE T04_S (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
extern void t04_s_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
// SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)
extern void recalc_08_(int*,int*,int*,int*,double*,double*,double*,double*);
// SUBROUTINE IGRF_GSW_08 (R,THETA,PHI,BR,BTHETA,BPHI)
extern void igrf_gsw_08_(double*,double*,double*,double*,double*,double*);
// SUBROUTINE XXXXX_08 (XSM,YSM,ZSM,XGSW,YGSW,ZGSW,J)
extern void smgsw_08_(double*,double*,double*,double*,double*,double*,int*);
extern void magsm_08_(double*,double*,double*,double*,double*,double*,int*);
extern void geigeo_08_(double*,double*,double*,double*,double*,double*,int*);
extern void geomag_08_(double*,double*,double*,double*,double*,double*,int*);
extern void gswgse_08_(double*,double*,double*,double*,double*,double*,int*);
extern void geogsw_08_(double*,double*,double*,double*,double*,double*,int*);

// define the index for functions
#define T04_S_IDX 0
#define RECALC_08_IDX 1
#define IGRF_GSW_08_IDX 2
#define SMGSW_08_IDX 3
#define MAGSM_08_IDX 4
#define GEIGEO_08_IDX 5
#define GEOMAG_08_IDX 6
#define GSWGSE_08_IDX 7
#define GEOGSW_08_IDX 8

// // Function Mapping
// typedef struct {
//     const char *func_name;
//     int func_idx;
// } func_map;
// // construct mapping
// func_map func_table[] = {
//     {"T04_s", T04_S_IDX},
//     {"RECALC_08", RECALC_08_IDX},
//     {"IGRF_GSW_08", IGRF_GSW_08_IDX},
//     {"SMGSW_08", SMGSW_08_IDX},
//     {"MAGSM_08", MAGSM_08_IDX},
//     {"GEIGEO_08", GEIGEO_08_IDX},
//     {"GEOMAG_08", GEOMAG_08_IDX},
//     {"GSWGSE_08", GSWGSE_08_IDX},
//     {"GEOGSW_08", GEOGSW_08_IDX}
// };
// // get the index of a certain function name
// int get_func_idx(const char* func_name) {
//     for (int i = 0; i < sizeof(func_table)/sizeof(func_table[0]); i++) {
//         if (strcmp(func_table[i].func_name, func_name) == 0) {
//             return func_table[i].func_idx;
//         }
//     }
//     return -1;  // return -1 if an invalid name is given
// }

void call_T04_s(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
    // input: prhs = [PARMOD,PS,X,Y,Z]
    // declaration
    int iopt = 0;
    size_t i;
    mxArray *oBX,*oBY,*oBZ;
    mxDouble *parmod,*ps,*x,*y,*z,*bx,*by,*bz;
    // check nlhs/nrhs
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("mex:igrfmex","T04: 5 inputs required: PARMOD,PS,X,Y,Z");
    }
    if(nlhs>3) {
        mexErrMsgIdAndTxt("mex:igrfmex","T04: too many outputs. Max 3.");
    }
    // check parmod
    if(mxGetNumberOfElements(prhs[0]) != 10){
        mexErrMsgIdAndTxt("mex:igrfmex","T04: len(Parmod) != 10");
    }
    // check XYZ dimensions
    mwSize dim_in = mxGetNumberOfDimensions(prhs[2]);
    if(dim_in != mxGetNumberOfDimensions(prhs[3]) || dim_in != mxGetNumberOfDimensions(prhs[4])){
        mexErrMsgIdAndTxt("mex:igrfmex","T04: X,Y,Z dimension dismatch");
    }
    const mwSize *sx = mxGetDimensions(prhs[2]);
    const mwSize *sy = mxGetDimensions(prhs[3]);
    const mwSize *sz = mxGetDimensions(prhs[4]);
    const size_t totnum = mxGetNumberOfElements(prhs[2]);
    for(i=0;i<dim_in;i++){
        if(sx[i]!=sy[i] || sx[i]!=sz[i]){
            mexErrMsgIdAndTxt("mex:igrfmex","T04: X,Y,Z dimension dismatch");
        }
    }
    // init BXBYBZ, to be assigned to plhs later
    oBX = mxDuplicateArray(prhs[2]);
    oBY = mxDuplicateArray(prhs[2]);
    oBZ = mxDuplicateArray(prhs[2]);
    // link in-data* to mxarray*
    parmod = mxGetDoubles(prhs[0]);
    ps = mxGetDoubles(prhs[1]);
    x = mxGetDoubles(prhs[2]);
    y = mxGetDoubles(prhs[3]);
    z = mxGetDoubles(prhs[4]);
    bx = mxGetDoubles(oBX);
    by = mxGetDoubles(oBY);
    bz = mxGetDoubles(oBZ);
    // CALCULATION
    for(i=0;i<totnum;i++){
        t04_s_(&iopt, parmod, ps, x+i, y+i, z+i, bx+i, by+i, bz+i);
    }
    // return process
    switch(nlhs){
        case 2:
            plhs[0] = oBX;
            plhs[1] = oBY;
            mxDestroyArray(oBZ);
            break;
        case 3:
            plhs[0] = oBX;
            plhs[1] = oBY;
            plhs[2] = oBZ;
            break;
    default: // 1 & 0
        plhs[0] = oBX;
        mxDestroyArray(oBY);
        mxDestroyArray(oBZ);
    }
}

void call_RECALC_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
    // input: prhs = [IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ]
    if(nrhs!=8){
        mexErrMsgIdAndTxt("mex:igrfmex","Recalc: Too few inputs => 8");
    }
    if(nlhs!=0){
        mexErrMsgIdAndTxt("mex:igrfmex","Recalc: Too much outputs => 0");
    }
    int yy,dd,hh,mm;
    double *ss;
    yy = mat2int(prhs[0]);
    dd = mat2int(prhs[1]);
    hh = mat2int(prhs[2]);
    mm = mat2int(prhs[3]);
    ss = mxGetDoubles(prhs[4]); //mat2int(prhs[4]);
    double *vx,*vy,*vz,ansout=0.;
    vx = mxGetDoubles(prhs[5]);
    vy = mxGetDoubles(prhs[6]);
    vz = mxGetDoubles(prhs[7]);
    
    recalc_08_(&yy,&dd,&hh,&mm,ss,vx,vy,vz);
    plhs[0] = mxCreateDoubleScalar(ansout);
}

void call_IGRF_GSW_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
    //input: prhs = [R,THETA,PHI]
    int i;
    if(nrhs!=3){
        mexErrMsgIdAndTxt("mex:igrfmex","IGRF_GSW: Too few inputs => 3");
    }
    if(nlhs>3){
        mexErrMsgIdAndTxt("mex:igrfmex","IGRF_GSW: num of outputs => 1-3");
    }
    // check XYZ dimensions
    mwSize dim_in = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *sz = mxGetDimensions(prhs[0]);
    const size_t totnum = mxGetNumberOfElements(prhs[0]);
    if(!matchecksize(2,dim_in,sz,prhs[1],prhs[2])){
        mexErrMsgIdAndTxt("mex:igrfmex","Recalc: X,Y,Z: dimension mismatch");
    }
    mxArray *oBX,*oBY,*oBZ;
    oBX = mxDuplicateArray(prhs[0]);
    oBY = mxDuplicateArray(prhs[0]);
    oBZ = mxDuplicateArray(prhs[0]);
    // get data
    double *x,*y,*z,*bx,*by,*bz;
    x = mxGetDoubles(prhs[0]);
    y = mxGetDoubles(prhs[1]);
    z = mxGetDoubles(prhs[2]);
    bx = mxGetDoubles(oBX);
    by = mxGetDoubles(oBY);
    bz = mxGetDoubles(oBZ);
    //CALCULATION
    for(i=0;i<totnum;i++){
        igrf_gsw_08_(x+i, y+i, z+i, bx+i, by+i, bz+i);
    }
    //arrange output
    switch(nlhs){
        case 2:
            plhs[0] = oBX;
            plhs[1] = oBY;
            mxDestroyArray(oBZ);
            break;
        case 3:
            plhs[0] = oBX;
            plhs[1] = oBY;
            plhs[2] = oBZ;
            break;
    default: // 1 & 0
        plhs[0] = oBX;
        mxDestroyArray(oBY);
        mxDestroyArray(oBZ);
    }
}

void call_TRANS_08(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[], 
                 void(*transfun)(double*,double*,double*,double*,double*,double*,int*)){
    // SUBROUTINE XXXXX_08 (XSM,YSM,ZSM,XGSW,YGSW,ZGSW,J)
    int i;
    int direction=1;
    if(nrhs==4){
        direction = mat2int(prhs[3]);
    }else if(nrhs!=3){
        mexErrMsgIdAndTxt("mex:igrfmex","TRANS: Too few inputs => 3,4");
    }
    if(nlhs>3){
        mexErrMsgIdAndTxt("mex:igrfmex","TRANS: num of outputs => 1-3");
    }
    // check XYZ dimensions
    mwSize dim_in = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *sz = mxGetDimensions(prhs[0]);
    const size_t totnum = mxGetNumberOfElements(prhs[0]);
    if(!matchecksize(2,dim_in,sz,prhs[1],prhs[2])){
        mexErrMsgIdAndTxt("mex:igrfmex","SMGSW: X,Y,Z: dimension mismatch");
    }
    // get data
    double *x1,*x2,*y1,*y2,*z1,*z2;
    x1 = mxGetDoubles(prhs[0]);
    y1 = mxGetDoubles(prhs[1]);
    z1 = mxGetDoubles(prhs[2]);
    
    mxArray *oBX,*oBY,*oBZ;
    oBX = mxDuplicateArray(prhs[0]);
    oBY = mxDuplicateArray(prhs[0]);
    oBZ = mxDuplicateArray(prhs[0]);
    
    x2 = mxGetDoubles(oBX);
    y2 = mxGetDoubles(oBY);
    z2 = mxGetDoubles(oBZ);
    
    //CALCULATION
    if(direction > 0){
        for(i=0;i<totnum;i++){
            transfun(x1+i, y1+i, z1+i, x2+i, y2+i, z2+i, &direction);
        }
    }else{
        for(i=0;i<totnum;i++){
            transfun(x2+i, y2+i, z2+i, x1+i, y1+i, z1+i, &direction);
        }
    }
    //arrange output
    switch(nlhs){
        case 2:
            plhs[0] = oBX;
            plhs[1] = oBY;
            mxDestroyArray(oBZ);
            break;
        case 3:
            plhs[0] = oBX;
            plhs[1] = oBY;
            plhs[2] = oBZ;
            break;
    default: // 1 & 0
        plhs[0] = oBX;
        mxDestroyArray(oBY);
        mxDestroyArray(oBZ);
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
    // first input is the function name
    if(nrhs<2){
        mexErrMsgIdAndTxt("mex:igrfmex","mexFunction: Too few params");
    }
    // // get function name
    // char* func_name = mxArrayToString(prhs[0]);
    // // get function index
    // int func_idx = get_func_idx(func_name);
    const int func_idx = mat2int(prhs[0]); 
    if (func_idx == -1) {
        mexErrMsgIdAndTxt("mex:igrfmex", "mexFunction: Invalid function name.");
    }

    // execution
    switch(func_idx){
        case T04_S_IDX:
            call_T04_s(nlhs,plhs,nrhs-1,prhs+1);
            break;
        case RECALC_08_IDX:
            call_RECALC_08(nlhs,plhs,nrhs-1,prhs+1);
            break;
        case IGRF_GSW_08_IDX:
            call_IGRF_GSW_08(nlhs,plhs,nrhs-1,prhs+1);
            break;
        case SMGSW_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,smgsw_08_);
            break;
        case MAGSM_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,magsm_08_);
            break;
        case GEIGEO_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,geigeo_08_);
            break;
        case GEOMAG_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,geomag_08_);
            break;
        case GSWGSE_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,gswgse_08_);
            break;
        case GEOGSW_08_IDX:
            call_TRANS_08(nlhs,plhs,nrhs-1,prhs+1,geogsw_08_);
            break;
        default:
        mexErrMsgIdAndTxt("mex:igrfmex","mexFunction: Function not recognized");
    }

    // mxFree(func_name);
}

int mat2int(const mxArray *mat){
    // convenient get int args from double/int64 mxArray
    void *modeptr;
    int mode;
    if(mxIsInt64(mat)){
        modeptr = mxGetInt64s(mat);
        mode = (int)*(mxInt64*)modeptr;
    }else if(mxIsDouble(mat)){
        modeptr = mxGetDoubles(mat);
        mode = (int)*(mxDouble*)modeptr;
    }else if(mxIsSingle(mat)){
        modeptr = mxGetSingles(mat);
    mode = (int)*(mxSingle*)modeptr;
    }
    else{
        mexErrMsgIdAndTxt("mex:igrfmex","mat2int: cannot convert to int");
    }
    return mode;
}

int matchecksize(int narray, mwSize ndim, const mwSize *sz, ...){
    // fast routine to check whether multiple matrices have same size
    va_list args;
    va_start(args, sz);
    const mxArray *mat;
    int i=0,j=0;
    mwSize ndim_;
    const mwSize *sz_;
    for(i=0;i<narray;i++){
        mat = va_arg(args, mxArray*);
        ndim_ = mxGetNumberOfDimensions(mat);
        if(ndim != ndim_){
            va_end(args);
            return 0;
        }
        sz_ = mxGetDimensions(mat);
        for(j=0;j<ndim;j++){
            if(sz[j]!=sz_[j]){
                va_end(args);
                return 0;
            }
        }
    }
    va_end(args);
    return 1;
}
