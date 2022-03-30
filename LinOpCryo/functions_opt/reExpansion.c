/*/ g = fprojectionAdjoint3D(proj, imageSize, r1,r2, Ti, Tp, LT, supp, shiftX, shiftY)
 * //
 * // 3D projection in parallel-beam geometry.
 * //
 * // cK:    Input image.
 * // r1, r2:projection coordinate.
 * // Ti:    sampling step in image domain Ti[3].
 * // Tp:    sampling step in the projection domain Tp[2].
 * // LT:    lookup table.
 * // numberOfSample: number Of sample in lookup table.
 * // supp : supprot of basis function in the projection domain*/
#include <pthread.h> /* for threading */
#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double *LT;
double supp, step;
double cx1, cx2, cx3, cy1, cy2, cy3, scale;
int numberOfSample;

typedef struct
{
    double *Ck, *Ck_s;
    int    *N,*M       ;
} tdata;


/*// function to read lookup table*/
double readLT(double u, double v, double w) {
    double r = sqrt(u*u+v*v+w*w);
    if (r > supp) {
        return 0;
    } else {
        r = r/step ;
        int rmin = (int)(r);
        rmin     = (rmin>r ? (rmin-1):rmin);
        int rmax = rmin+1    ;
        double p = r-rmin    ;
        return p * LT[rmax] + (1-p) * LT[rmin];
        /*//int rround = round(r);
         * //return LT[rmin];*/
    }
}
int myFloor(double u){
    
    int x = (int) u;
    if (x<=u)
        return x  ;
    else
        return x+1;
    
}

void *reExpand(void * t_data)
{
    tdata *th_data = (tdata *) t_data;
    double *Ck, *Ck_s;
    int    *N, *M      ;
    
    Ck    = (*th_data).Ck;
    Ck_s  = (*th_data).Ck_s;
    N     = (*th_data).N;
    M     = (*th_data).M;
    /*// mexPrintf("\n index is %d ", index);*/
    
    int n1, n2, n3, ns1, ns2, ns3;
    int i;
    double x1n,x2n,x3n;
    
    for (n3=0; n3<N[2]; n3++) {
        /*/ initial computation*/
        x3n  = (n3-cx3);
        
        for (n1=0; n1<N[1]; n1++) {
            /*/ initial computation*/
            x1n  = (n1-cx1);
            
            for (n2=0; n2<N[0]; n2++) {
                /*/ initial computation*/
                x2n  = (n2-cx2);
                
                for (ns1=0; ns1<M[1]; ns1++) {
                    for (ns2=0; ns2<M[0]; ns2++) {
                        for (ns3=0; ns3<=M[2]; ns3++) {
                            Ck[(n3*((int) N[1])+n1)*((int) N[0])+n2] += Ck_s[(ns3*((int) M[1])+ns1)*((int) M[0])+ns2] * readLT((x1n/scale)-(ns1-cy1), (x2n/scale)-(ns2-cy2), (x3n/scale)-(ns3-cy3));
                        }
                    }
                }
            }
        }
        
    }
}
/*/ mex file to implent the forward projection*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize *M;
    int N[3];
    double *Ck, *Ck_s;
    /*/ Check for proper number of input arguments*/
    if (nrhs!=4) mexErrMsgTxt("There must be exactly 4 input arguments.");
    if (nlhs!=1) mexErrMsgTxt("There must be exactly 1 output argument.");
    
    /*/ Inputs*/
    Ck_s = mxGetPr(prhs[0]);
    LT = mxGetPr(prhs[1]);
    scale = mxGetScalar(prhs[2]);
    supp = mxGetScalar(prhs[3]);
    step = supp/((double)(numberOfSample-1));
    
    numberOfSample = mxGetNumberOfElements(prhs[1]);
    /*/ Check for proper type of input arguments*/
    
    /* larger scale*/
    M = mxGetDimensions(prhs[0]);
    
    /* smaller scale*/
    N[0] = ceil(M[0]*scale);
    N[1] = ceil(M[1]*scale);
    N[2] = ceil(M[2]*scale);
    

    /*/ Allocate output variable*/
    plhs[0] = mxCreateNumericArray(3, M, mxDOUBLE_CLASS, mxREAL);
    if (!plhs[0]) mexErrMsgTxt("Could not allocate memory for output variable.");
    Ck = mxGetPr(plhs[0]);
    
    /*/ Center of object smaller order scale*/
    cx1 = ((double) N[1]-1)/2;
    cx2 = ((double) N[0]-1)/2;
    cx3 = ((double) N[2]-1)/2;
    
    /*/ Center of object higher order scale*/
    cy1 = ((double) M[1]-1)/2;
    cy2 = ((double) M[0]-1)/2;
    cy3 = ((double) M[2]-1)/2;

    int i;
    /*  call the C subroutine */
    tdata *data;
    data = (tdata*)  mxCalloc( 1,sizeof(tdata) );
    (*data).Ck     =  Ck  ;
    (*data).Ck_s   =  Ck_s  ;
    (*data).N      =  (int *)N;
    (*data).M      =  (int *)M;
    for (i=0 ; i<3 ; i++){
        (*data).M[i] = M[i];
    }
    reExpand(data);
    
}