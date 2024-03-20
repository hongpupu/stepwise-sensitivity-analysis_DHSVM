/*
Code Reference:
Tang, Yong, Reed, P. M., Wagener, T., and van Werkhoven, K., "Comparing 
sensitivity analysis methods to advance lumped watershed model identification 
and evaluation." Hydrology and Earth System Sciences, 11, 793-817, 2007.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int *IntVector(int n)
{
    int *v;

    v=(int *) calloc( (size_t) n, (size_t) sizeof(int));
    return v;
}

int FreeIntVector(int *v)
{
   free((char *) v);
   return 0;
}

double *DoubleVector(int n)
{
    double *v;
   
    v = (double *) calloc((size_t) n, (size_t) sizeof(double));
    return v;
}

int FreeDoubleVector(double *v)
{
   free((char *) v);
   return 0;
}

float *FloatVector(int n)
{
    float *v;
   
    v = (float *) calloc((size_t) n, (size_t) sizeof(float));
    return v;
}

int FreeFloatVector(float *v)
{
   free((char *) v);
   return 0;
}

int **IntMatrix(int m, int n)
{
   register int i;
   int **x;

   x=(int **) calloc((size_t) m, (size_t) sizeof(int *));
   for(i=0;i<m;i++) {
      x[i]=(int *) calloc((size_t) n, (size_t) sizeof(int));
   }
   return x;
}

int FreeIntMatrix (int **x, int m)
{
   register int i;

   for (i=0;i<m;i++) free((char *) x[i]);
   free((char *) x);
   return 0;
}

float **FloatMatrix(int m, int n)
{
   register int i;
   float **x;

   x=(float **) calloc((size_t) m, (size_t) sizeof(float *));
   for(i=0;i<m;i++) {
      x[i]=(float *) calloc((size_t) n, (size_t) sizeof(float));
   }
   return x;
}

int FreeFloatMatrix (float **x, int m)
{
   register int i;

   for (i=0;i<m;i++) free((char *) x[i]);
   free((char *) x);
   return 0;
}

double **DoubleMatrix(int m, int n)
{
   register int i;
   double **x;

   x=(double **) calloc((size_t) m, (size_t) sizeof(double *));
   for(i=0;i<m;i++) {
      x[i]=(double *) calloc((size_t) n, (size_t) sizeof(double));
   }
   return x;
}

int FreeDoubleMatrix (double **x, int m)
{
   register int i;

   for (i=0;i<m;i++) free((char *) x[i]);
   free((char *) x);
   return 0;
}

double Round(double x)
{
  double y;
  y = ceil(x);
  if ((ceil(x)-x) <= 0.5) return(y);
  return(floor(x));
}

double RandomDbl()
{
	double x;
	x=1.0/(RAND_MAX)*(double) rand();
	return (x);
}

double RandomDblRge(double lb, double ub)
{
	double x;
	x=lb+(ub-lb)/(RAND_MAX)*(double) rand();
	return (x);
}

//used to generate hadamard matrix n*n
//The size n must be of the form 2^k*p for p=1, 12 or 20
//The code is translated from matlab
void hadamard(int **H, int n)
{
	int i,j,k;
	int e,p;
	//base hadamard
	int h1=1;
	int h12[12][12]={{1,1,1,1,1,1,1,1,1,1,1,1},
	{1,-1, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1},
	{1,-1,-1, 1,-1, 1, 1, 1,-1,-1,-1, 1},
	{1, 1,-1,-1, 1,-1, 1, 1, 1,-1,-1,-1},
	{1,-1, 1,-1,-1, 1,-1, 1, 1, 1,-1,-1},
	{1,-1,-1, 1,-1,-1, 1,-1, 1, 1, 1,-1},
	{1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1, 1},
	{1, 1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1},
	{1, 1, 1,-1,-1,-1, 1,-1,-1, 1,-1, 1},
	{1, 1, 1, 1,-1,-1,-1, 1,-1,-1, 1,-1},
	{1,-1, 1, 1, 1,-1,-1,-1, 1,-1,-1, 1},
	{1, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1,-1}};
	int h20[20][20]={{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
	{1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1},
	{1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1},
	{1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1},
	{1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1},
	{1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1},
	{1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1},
	{1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1},
	{1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1},
	{1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1},
	{1,-1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1},
	{1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1},
	{1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1},
	{1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1},
	{1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1},
	{1, 1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1},
	{1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1},
	{1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1},
	{1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1},
	{1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1}};
	
	//find k if n = 2^k*p
	e = 0;
	while (n>1 && (n/2)*2 == n)
	{
		e++;
		n=n/2;
	}
	
	if (n!=1) e-=2; // except for n=2^k, need a multiple of 4
	if (e<0) n=-1; // trigger error if not a multiple of 4
	
	//Kronecker product construction
	if (n==1)
	{
	    
		p=1;
		H[0][0]=h1;
		for(i=0;i<e;i++)
		{
			for(j=0;j<(int)pow(2,i);j++)
				for(k=0;k<(int)pow(2,i);k++)
				{
					
					H[j][k+p*(int)pow(2,i)]=H[j][k];
					H[j+p*(int)pow(2,i)][k]=H[j][k];
					H[j+p*(int)pow(2,i)][k+p*(int)pow(2,i)]=-H[j][k];
				}

		}
	}
	else if(n==3)
	{
		p=12;
		for(j=0;j<p;j++)
			for(k=0;k<p;k++)
				H[j][k]=h12[j][k];
		for(i=0;i<e;i++)
		{
			for(j=0;j<(int)pow(2,i);j++)
				for(k=0;k<(int)pow(2,i);k++)
				{
					
					H[j][k+p*(int)pow(2,i)]=H[j][k];
					H[j+p*(int)pow(2,i)][k]=H[j][k];
					H[j+p*(int)pow(2,i)][k+p*(int)pow(2,i)]=-H[j][k];
				}

		}
	}
	else if(n==5)
	{
		p=20;
		for(j=0;j<p;j++)
			for(k=0;k<p;k++)
				H[j][k]=h20[j][k];
		for(i=0;i<e;i++)
		{
			for(j=0;j<(int)pow(2,i);j++)
				for(k=0;k<(int)pow(2,i);k++)
				{
					
					H[j][k+p*(int)pow(2,i)]=H[j][k];
					H[j+p*(int)pow(2,i)][k]=H[j][k];
					H[j+p*(int)pow(2,i)][k+p*(int)pow(2,i)]=-H[j][k];
				}

		}
	}
	else
	{
		printf("n must be 2^e*p, for p = 1, 12, 20 ");
		exit(-1);
	}

}

//Permuation, random select n numbers from 1,...m, m>=n, without replacement
void permutate(int *P, int m, int n)
{
	int i,j,k;
	int *index1, *index2;
    int element, nPos;
	
	index1 = IntVector(m);
    index2 = IntVector(m);
    
    for (i=0; i<m; i++)   
	{
		index1[i] = 1;
		index2[i] = 0;
    }
	for(i=0;i<n;i++)
	{
		nPos = 0;
		for (j=0; j<m; j++)   
			if (index1[j] > 0)   
			{
				index2[nPos] = (j+1);
				nPos = nPos+1;
			}
		
		k = -1;
		while (k < 1)  k = (int)Round(0.5+nPos*RandomDbl());
		element = index2[k-1];
		index1[element-1] = 0;
		P[i] = element;
	}
	FreeIntVector(index1);
	FreeIntVector(index2);
}

//group parameters, return a nrep*nparm matrix
void groupMap(int **group, int nrep, int nparm, int ncol,int nzero)
{
	int i,j,k,kk;
	int *P;

	P=IntVector(nzero); 

	for(i=0;i<nparm;i++)
	{
		//select the zero replicates
		permutate(P,nrep,nzero);
		for(j=0;j<nzero;j++)
		{
			group[P[j]-1][i]=0;
		}
		//assign group and orientation randomly to the remaining replicates
		for(j=0;j<nrep;j++)
		{
			kk=0;
			for(k=0;k<nzero;k++)
				if(P[k]==j+1) kk=-1;
			if(kk==0)
			{
				kk=-1;
				while (kk < 1)  kk = (int)Round(0.5+ncol*RandomDbl());
				if(RandomDbl()>0.5)
					group[j][i]=kk;
				else
					group[j][i]=-kk;
			}
		}
	}
	
    FreeIntVector(P);
}

//translate -1 and 1 in a Hadamard matrix to the levels 0 and 1, lower, upper level
void FFLevels(int **H,int n)
{
	int i,j;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			if(H[i][j]==-1) H[i][j]=0;	
}

//fold samples with integer type
void foldIntSample(int **oldSample, int **newSample, int m, int n)
{
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			newSample[i][j]=oldSample[i][j];
	for(i=n;i<2*m;i++)
		for(j=0;j<n;j++)
			newSample[i][j]=1-oldSample[i-m][j];
}

//Iterated Fractional Factorial Design 
//nrep-number of replicates, repeat LHS process nrep times
//ncol-size of Hadamard matrix, define the number of levels, must be of the form 2^k*p for p=1, 12 or 20
//nzero-number of replicates that take middle level, nzero<=nrep
//npar-number of parameters
//lb,ub-lower and upper bound of parameter values
//LL is the factor level matrix, can be used in ANOVA
void IFFD(int nrep,int nparm,int ncol,int nzero,double *lb, double *ub,double **data,int **LL)
{
	int i,j,k;
	double l,u,midwid;
	int **group;
	int *Levels;
	int **H1,**H2;
	int nrow;

	nrow=ncol*2;
	group = IntMatrix(nrep, nparm);
	Levels=IntVector(2*ncol);

	H1=IntMatrix(ncol,ncol);
	H2=IntMatrix(nrow,ncol);
	
	//assign parameter groups
	groupMap(group,nrep,nparm,ncol,nzero);
	//create Hadamard matrix
	hadamard(H1,ncol);
	//translate to FF levels
	FFLevels(H1,ncol);
	//fold samples
	foldIntSample(H1,H2,ncol,ncol);
		
	//Generate Samples
	for(i=0;i<nrep;i++)
	{
		for(j=0;j<nparm;j++)
		{
			midwid=nzero/((double)nrep)*(ub[j]-lb[j]);
			if(group[i][j]==0)
			{
				l=(lb[j]+ub[j])/2-midwid/2;
				u=(lb[j]+ub[j])/2+midwid/2;
				for(k=0;k<nrow;k++)
				{
					data[i*nrow+k][j]=RandomDblRge(l,u);
					LL[i*nrow+k][j]=0;
				}
			}
			else if(group[i][j]>0)
			{
				for(k=0;k<nrow;k++)
				{
					Levels[k]=H2[k][group[i][j]-1];
				}
				//generate sampels
				for(k=0;k<nrow;k++)
				{
					if(Levels[k]==0)
					{
						l=lb[j];
						u=(lb[j]+ub[j])/2-midwid/2;
						LL[i*nrow+k][j]=-1;
					}
					else if(Levels[k]==1)
					{
						l=(lb[j]+ub[j])/2+midwid/2;
						u=ub[j];
						LL[i*nrow+k][j]=1;
					}
					data[i*nrow+k][j]=RandomDblRge(l,u);
				}
			}
			else
			{
				for(k=0;k<nrow;k++)
				{
					Levels[k]=1-H2[k][abs(group[i][j])-1];
				}
				
				//generate sampels
				for(k=0;k<nrow;k++)
				{
					if(Levels[k]==0)
					{
						l=lb[j];
						u=(lb[j]+ub[j])/2-midwid/2;
						LL[i*nrow+k][j]=-1;
					}
					else if(Levels[k]==1)
					{
						l=(lb[j]+ub[j])/2+midwid/2;
						u=ub[j];
						LL[i*nrow+k][j]=1;
					}
					data[i*nrow+k][j]=RandomDblRge(l,u);
				}
			}
		}
	}
	
	FreeIntMatrix(group,nrep);
	FreeIntVector(Levels);
}

//Caculate F values for main effects and 2-way interactions based on anova
//objs is the array storing the output responses (objectives) for nsample
//LL is the factor Level matrix, -1, 0, and +1
//F is the array storing F values
void ANOVA(float *objs,int **LL, int nsample, int nparam, float *F, float *R2)
{
	int i,j,k,l;
	float Yi[3],Yj[3],GY;//Level mean, grand mean
	int ni[3],nj[3];//number of points at each level
	float Yij[3][3];
	int nij[3][3];
	float SSTR;//treatment sum of squares
	float SSTO;
	float SSE;//error sum of squares
	float MSTR;//treatmeant mean squre 
	float MSE;//error eman squre
	float SSAB;//cross-factor sum of squares
	float SSTRAB;
	float SSEAB;
	float MSAB;
	float MSEAB;
	float s; //sum of objs()
	int index;
	
	s=0;
	for(i=0;i<nsample;i++)
	{
		s=s+objs[i];
	}
	GY=s/nsample;
	
	SSTO=0;
	for(i=0;i<nsample;i++)
	{
		SSTO=SSTO+(objs[i]-GY)*(objs[i]-GY);
	}
	R2[0]=0;
	R2[1]=0;
	//main effects
	for(i=0;i<nparam;i++)
	{
		for(j=0;j<3;j++) 
		{
			Yi[j]=0;
			ni[j]=0;
		}
		for(j=0;j<nsample;j++)
		{
			Yi[LL[j][i]+1]=Yi[LL[j][i]+1]+objs[j];
			ni[LL[j][i]+1]=ni[LL[j][i]+1]+1;
		}
		SSTR=0;
		for(j=0;j<3;j++) 
		{
			if(ni[j]>0) Yi[j]=Yi[j]/ni[j]; 
			SSTR=SSTR+ni[j]*(Yi[j]-GY)*(Yi[j]-GY); 
		}
		SSE=0;
		for(j=0;j<nsample;j++)
		{
			SSE=SSE+(objs[j]-Yi[LL[j][i]+1])*(objs[j]-Yi[LL[j][i]+1]);
		}
		MSTR=SSTR/(3-1);
		MSE=SSE/(nsample-3);
		if(MSE>0) F[i]=MSTR/MSE;
		R2[0]=R2[0]+SSTR/SSTO;
		R2[1]=R2[1]+SSTR/SSTO;
	}
	//Interactions
	index=nparam;
	for(i=0;i<nparam-1;i++)
	{
		for(j=0;j<3;j++) 
		{
			Yi[j]=0;
			ni[j]=0;
		}
		for(j=0;j<nsample;j++)
		{
			Yi[LL[j][i]+1]=Yi[LL[j][i]+1]+objs[j];
			ni[LL[j][i]+1]=ni[LL[j][i]+1]+1;
		}
		for(j=0;j<3;j++) 
		{
			if(ni[j]>0) Yi[j]=Yi[j]/ni[j]; 
		}

		for(j=i+1;j<nparam;j++)
		{
			for(k=0;k<3;k++) 
			{
				Yj[k]=0;
				nj[k]=0;
			}
			for(k=0;k<nsample;k++)
			{
				Yj[LL[k][j]+1]=Yj[LL[k][j]+1]+objs[k];
				nj[LL[k][j]+1]=nj[LL[k][j]+1]+1;
			}
			for(k=0;k<3;k++) 
			{
				if(nj[k]>0) Yj[k]=Yj[k]/nj[k]; 
			}

			for(k=0;k<3;k++)
			{
				for(l=0;l<3;l++)
				{
					Yij[k][l]=0;
					nij[k][l]=0;
				}
			}
			for(k=0;k<nsample;k++)
			{
				Yij[LL[k][i]+1][LL[k][j]+1]=Yij[LL[k][i]+1][LL[k][j]+1]+objs[k];
				nij[LL[k][i]+1][LL[k][j]+1]=nij[LL[k][i]+1][LL[k][j]+1]+1;
			}
			SSTRAB=0;
			for(k=0;k<3;k++)
			{
				for(l=0;l<3;l++)
				{
					if(nij[k][l]>0) Yij[k][l]=Yij[k][l]/nij[k][l];
					SSTRAB=SSTRAB+nij[k][l]*(Yij[k][l]-GY)*(Yij[k][l]-GY);
				}
			}
			SSAB=0;
			for(k=0;k<3;k++)
			{
				for(l=0;l<3;l++)
				{
					SSAB=SSAB+nij[k][l]*(Yij[k][l]-Yi[k]-Yj[l]+GY)*(Yij[k][l]-Yi[k]-Yj[l]+GY);
				}
			}
			MSAB=SSAB/((3-1)*(3-1));
			SSEAB=0;
			for(k=0;k<nsample;k++)
			{
				SSEAB=SSEAB+(objs[k]-Yij[LL[k][i]+1][LL[k][j]+1])*(objs[k]-Yij[LL[k][i]+1][LL[k][j]+1]);
			}
			MSEAB=SSEAB/(nsample-3*3);
			if(MSEAB>0) F[index]=MSAB/MSEAB;
			R2[1]=R2[1]+SSAB/SSTO;
			index=index+1;
		}
	}
	
}

void Bt_ANOVA(int nsample, int nparam, int np, int nrspl, int **LL, float *stobj, float *FCI)
{
	int i,j,k;
	float *obj;
	int rspl;
	int **LL1;
	float **s;
	float *ss,*sss;
	float R2[2];
	int ii,kk;

	s=FloatMatrix(nrspl,np);
	ss=FloatVector(np);
	sss=FloatVector(np);
	LL1=IntMatrix(nsample,nparam);
	obj=FloatVector(nsample);

			for(ii=0;ii<np;ii++)
			{
				ss[ii]=0;
				sss[ii]=0;
			}
			for(rspl=0;rspl<nrspl;rspl++)
			{
				for(k=0;k<nsample;k++)
				{
					kk=(int)Round((nsample-1)*RandomDbl());
					obj[k]=stobj[kk];
					for(ii=0;ii<nparam;ii++) LL1[k][ii]=LL[kk][ii];
				}
				ANOVA(obj,LL1,nsample,nparam,s[rspl],R2);
				for(ii=0;ii<np;ii++) ss[ii]=ss[ii]+s[rspl][ii];
			}
			for(ii=0;ii<np;ii++) ss[ii]=ss[ii]/nrspl;
			for(rspl=0;rspl<nrspl;rspl++)
			{
				for(ii=0;ii<np;ii++)
					sss[ii]=sss[ii]+(s[rspl][ii]-ss[ii])*(s[rspl][ii]-ss[ii]);
			}
			for(ii=0;ii<np;ii++)
			{
				FCI[ii]=1.96*sqrt(sss[ii]/(nrspl-1));
			}

	FreeFloatMatrix(s,nrspl);
	FreeIntMatrix(LL1,nsample);
	FreeFloatVector(ss);
	FreeFloatVector(sss);
	FreeFloatVector(obj);
}

int main(int argc, char* argv[]) {
	int nrep = 100; //number of replications
	int ncol = 20; //Hadamard matrix size (2, 12 or 20)
	int nzero = 15; //number of replicates that take middle level, useful for IFFD
	int nsample=nrep*(2*ncol); //# samples
	int nparam=2; //# parameters
	int np=nparam+nparam*(nparam-1)/2; //# of parameters and parameter interactions
	int nrspl=1000; //bootstrap resamples
	double* lb = DoubleVector(nparam); //the parameter lower bounds
	double* ub = DoubleVector(nparam); //the parameter upper bounds
	double** par=DoubleMatrix(nsample,nparam); //parameter matrix
	int** LL=IntMatrix(nsample,nparam); //parameter Level matrix
	float* obj=FloatVector(nsample); //objective (model response) vector
	float* F = FloatVector(np); //the F values assigned by ANOVA
	float* R2 = FloatVector(2); //the R2 values assigned by ANOVA
	float* FCI = FloatVector(np); //the confidence interval values assigned by ANOVA bootstrapping
	int i, j, k;

	//initialize random seed
	srand((unsigned)time( NULL )); //for time-based seeding, use 

	//initialize parameter lower and upper bounds
	for (i=0; i<nparam; i++) {
		lb[i] = 0.0;
		ub[i] = 1.0;
	}
	
	//generate IFFD samples
	IFFD(nrep,nparam,ncol,nzero,lb,ub,par,LL);

	//***insert code here to run IFFD parameters through model and assign reponses to obj***
	for (i=0; i<nsample; i++) {
		obj[i] = par[i][0]*par[i][1];
	}

	//calculate F values according ANOVA
	ANOVA(obj,LL,nsample,nparam,F,R2);

	//bootstrap ANOVA
	Bt_ANOVA(nsample,nparam,np,nrspl,LL,obj,FCI);

	//output results
	printf("Main Effects\n");
	for (i=0; i<nparam; i++) {
		printf("  Parameter %d:\t%10.4f [ %6.4f ]\n", i, F[i], FCI[i]);
	}
	printf("Interaction Effects\n");
	for(j=0;j<nparam-1;j++) {
		for(k=j+1;k<nparam;k++) {
			printf("  Parameters %d & %d:\t%10.4f [ %6.4f ]\n", j, k, F[i], FCI[i]);
			i++;
		}
	}
	printf("R2-1: %6.4f\n", R2[0]);
	printf("R2-2: %6.4f\n", R2[1]);

	//free objects
	FreeFloatVector(F);
	FreeFloatVector(FCI);
	FreeFloatVector(R2);
	FreeFloatVector(obj);
	FreeDoubleMatrix(par,nsample);
	FreeIntMatrix(LL,nsample);
	FreeDoubleVector(lb);
	FreeDoubleVector(ub);

	return 0;
}
