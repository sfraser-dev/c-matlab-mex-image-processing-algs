/* MED2_M.cmex */

/* MATLAB USAGE: matrixOut=MED2_M(matrixIn,windowSize) */

#include <math.h>
#include "mex.h"

/* prototypes */
void filter(double**, double**, int, int, int);
double fill(double**,int,int,int,int,int,int);
double median(double*,int);

/* mex 'main' function */
void mexFunction(	int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])	{
				
	int no_rows, no_cols;	/* number of rows and columns (argument 1) */
	int r,c;		/* for loop variables */
	int trv;		/* traverse mex kernelays */
	int size_element;	/* size of kernelay element (for malloc) */
	double *in_handle, *out_handle;	/* copying mex 'fortran' kernelays */ 
	double *in_kernelay, *out_kernelay;	/* dynamic kernelays (1-D) */
	double **in, **out;		/* dynamic kernelays (quasi 2-D) */
	char *err_msg;			/* an error message string */
	int ws;				/* window size */

	/* checking number of inputs */
	if(nrhs !=2)
		mexErrMsgTxt("Must have three input arguments");
	if(nlhs !=1)
		mexErrMsgTxt("Must have one output argument");
		
	/* preventing sparse, complex and string matrices */
	if( mxIsComplex(prhs[0])||!(mxIsClass(prhs[0],"double"))
	    || mxIsClass(prhs[0],"sparse") || mxIsChar(prhs[0]) ){
		err_msg="input must be real, double, full and non-string";
		mexErrMsgTxt(err_msg);
	}
								
	/* getting the number of rows and columns from input matrix */
	no_rows=mxGetM(prhs[0]);
	no_cols=mxGetN(prhs[0]);
	
	/* getting argument two (window size) */
	ws=(int)mxGetScalar(prhs[1]);
		
	/* creating an output kernelay, giving it a handle */
	plhs[0]=mxCreateDoubleMatrix(no_rows,no_cols,mxREAL);
	out_handle=mxGetPr(plhs[0]);
	
	/* creating dynamic 2d kernelays */
	in_kernelay = (double*) mxMalloc (no_rows*no_cols*sizeof(double));
	in = (double**) mxMalloc (no_rows*sizeof(double*));
	out_kernelay = (double*) mxMalloc (no_rows*no_cols*sizeof(double));
	out = (double**) mxMalloc (no_rows*sizeof(double*));
	for(r=0; r<no_rows; r++){
		in[r]=&(in_kernelay[r*no_cols]);
		out[r]=&(out_kernelay[r*no_cols]);
	}
	
	/* creating an input kernelay, giving it a handle */
	in_handle=mxGetPr(prhs[0]);
	
	/*************************************************************************
	** copying 1-d input kernelay into 2-d matrix, note that since MATLAB	**		
	** was originally written in fortran, it treats kernelays in a manner 	**
	** similiar to fortran. Thus, in MATLAB, 2-d matrices are stored 	**
	** as 1-d matrices in the following way:				**
	**									**
	** (2-d)                   (1-d)					**
	** 1 2 3								**
	** 4 5 6	= 	1 4 2 5 3 6					**
	*************************************************************************/
	trv=0;	
	for(c=0; c<no_cols; c++){
		for(r=0; r<no_rows; r++)
			in[r][c]=in_handle[trv++];
	}

	/* PROCESS THE MATRIX/ARRAY HERE */
	filter(in,out,no_rows,no_cols,ws);		
	
	/* copying the dynamic matrix to the output handle */
	trv=0;
	for(c=0; c<no_cols; c++){
		for(r=0; r<no_rows; r++)
			out_handle[trv++]=out[r][c];
	}

	/* freeing dynamic memory */
	mxFree(in[0]); in_kernelay=0;
	mxFree(in); in=0;
	mxFree(out[0]); out_kernelay=0;
	mxFree(out); out=0;
	
	mexUnlock();		/* allows for re-compiling */
}

/* perform filtering */
void filter(double **m_in, double **m_out, int no_rows, int no_cols, int ws){
	int curRow, curCol;
	int side, scale;
	
	side=ws; 		/* size of kernel side (assumed square) */
	scale=(int)(side-1)/2;	/* "width" of kernel */
		
	/* filling output matrix by first filling kernel values and then
	 * processing these kernel values */
	for(curRow=0; curRow<no_rows; curRow++){
		for(curCol=0; curCol<no_cols; curCol++)		
			m_out[curRow][curCol]=
				fill(m_in,curRow,curCol,side,scale,no_rows,no_cols);
	}
}

/* "fill" kernel */
double fill(double **m_in,int curRow,int curCol,
	int side,int scale,int no_rows,int no_cols){
				
	int length=side*side;
	int rPos, cPos, rDiff, cDiff, r, c;
	double *kernel_kernelay;
	int retVal;
		
	/* creating dynamic 1-D kernelay */
	kernel_kernelay = (double*) mxMalloc (no_rows*no_cols*sizeof(double));
		
	/* filling the kernel_kernelay */
	for(r=0; r<side; r++){
	    	rPos=curRow-scale+r;
	    	for(c=0; c<side; c++){
			cPos=curCol-scale+c;
			/* mirror padding */	
			if(rPos<0)			
			  	rPos=(sqrt(pow(rPos,2)))-1;
			if(rPos>=no_rows){
		   		rDiff=rPos-no_rows;
		   		rPos=no_rows-rDiff-1;
			}
			if(cPos<0)
		     		cPos=(sqrt(pow(cPos,2)))-1;
	 		if(cPos>=no_cols){
		   		cDiff=cPos-no_cols;
				cPos=no_cols-cDiff-1;
			}
			kernel_kernelay[(side*r)+c]=m_in[rPos][cPos];							
	    	}
	}
				
	/* processing the values within the kernel */	
	retVal = median(kernel_kernelay,length);
	mxFree(kernel_kernelay); kernel_kernelay=0;
	return retVal;
}
	
#define ELEM_SWAP(a,b) { double temp =(a); (a)=(b); (b)=temp; }

/* "process" the kernel values with quickselect routine */
double median(double* kernel, int length){
        int low, high;
        int median;
        int middle, ll, hh;

        low=0; high=length-1; median=((low+high)/2);
        for(;;) {
                if(high <= low) return kernel[median];
                if(high == low+1) {
                        if(kernel[low] > kernel[high])
                                ELEM_SWAP(kernel[low],kernel[high]);
                        return kernel[median];
                }
                middle = (low+high)/2;
                if(kernel[middle] > kernel[high]) ELEM_SWAP(kernel[middle],kernel[high]);
                if(kernel[low] > kernel[high]) ELEM_SWAP(kernel[low],kernel[high]);
		if(kernel[middle] > kernel[low]) ELEM_SWAP(kernel[middle],kernel[low]);
                ELEM_SWAP(kernel[middle],kernel[low+1]);
                ll=low+1;
                hh=high;
                for(;;) {
                        do ll++; while (kernel[low] > kernel[ll]);
                        do hh--; while (kernel[hh] > kernel[low]);
                        if(hh < ll) break;
                        ELEM_SWAP(kernel[ll],kernel[hh]);
                }
                ELEM_SWAP(kernel[low],kernel[hh]);
                if(hh <= median) low=ll;
                if(hh >= median) high=hh-1;
        }
}
