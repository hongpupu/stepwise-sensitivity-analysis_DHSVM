/*
Reference:
Tang, Yong, Reed, P. M., Wagener, T., and van Werkhoven, K., "Comparing 
sensitivity analysis methods to advance lumped watershed model identification 
and evaluation." Hydrology and Earth System Sciences, 11, 793-817, 2007.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

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
	float s;
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

void ANOVAA(float *objs,int **LL, int nsample, int nparam, float *F, float *R2)
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
	float s;
	int index;
	
    
	s=0;
	for(i=0;i<nsample;i++)
	{
	//%printf("%lf", objs[i]);
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

        ///////////////////////////////////////////
        FILE*main;
        main = fopen("main01","w");
        if (main==NULL)
        { 
            printf("error main\n");
            exit(1);
        }
        ///////////////////////////////////////////

	//main effects
	for(i=0;i<nparam;i++)
	{
		printf("%d", nparam);
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
                //printf("parameter%d: %f\t%f\n",i,SSTR,SSE);
                fprintf(main,"%s%d\t%f\t%f\n","Parameter",i,SSTR,SSE);///////////////////////////

		MSTR=SSTR/(3-1);
		MSE=SSE/(nsample-3);
		if(MSE>0) F[i]=MSTR/MSE;
		R2[0]=R2[0]+SSTR/SSTO;
		R2[1]=R2[1]+SSTR/SSTO;
	}
        fclose(main);///////////////////////////////////

        ///////////////////////////////////////////
        FILE*interactions;
        interactions = fopen("interactions01","w");
        if (interactions==NULL)
        { 
            printf("error interactions\n");
            exit(1);
        }
        ///////////////////////////////////////////

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

                        fprintf(interactions,"%s%d%s%d\t\t%f\t%f\n","Parameter",i,"&",j,SSAB,SSEAB);///////////////////////////
                        //printf("parameter%d&%d: %f\t%f\n",i,j,SSAB,SSEAB);

			MSEAB=SSEAB/(nsample-3*3);
			if(MSEAB>0) F[index]=MSAB/MSEAB;
			R2[1]=R2[1]+SSAB/SSTO;
			index=index+1;
		}
	}
        fclose(interactions);///////////////////////////////////
	
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

int main(int argc, char* argv[]) 
{
	int nrep = 500; //number of replications
	int ncol = 20; //Hadamard matrix size (2, 12 or 20)
	int nzero = 20; //number of replicates that take middle level, useful for IFFD
	int nsample=12916; //# samples
	int nparam=104; //# parameters,you can change depends on specific situation
	int np=nparam+nparam*(nparam-1)/2; //# of parameters and parameter interactions
	int nrspl=1000; //bootstrap resamples 1000
	double* lb = DoubleVector(nparam); //the parameter lower bounds
	double* ub = DoubleVector(nparam); //the parameter upper bounds
	double** par=DoubleMatrix(nsample,nparam); //parameter matrix
	int** LL=IntMatrix(nsample,nparam); //parameter Level matrix
	float* obj=FloatVector(nsample); //objective (model response) vector
	float* F = FloatVector(np); //the F values assigned by ANOVA
	float* R2 = FloatVector(2); //the R2 values assigned by ANOVA
	float* FCI = FloatVector(np); //the confidence interval values assigned by ANOVA bootstrapping
	int i, j, k;

        if (argc != 2)
        {
            fprintf(stderr, "\nUsage: ./anova-final  input_filename\n\n");
            exit(EXIT_FAILURE);
        }

        //read samplefile "par"
        FILE*samplefile;
        samplefile=fopen("final_par_no_runnumber","r");
        if (samplefile==NULL)
        {
                printf("fail to open samplefile\n");
                exit(1);
        }
        for (i=0; i<nsample; i++) 
        {
             for (j=0; j<nparam; j++)
             {                                            
                  fscanf(samplefile,"%lf",&par[i][j]);             
             }
        } 
        fclose(samplefile);
        //read parameter Level matrix "LL"
        FILE*LL_file;
        LL_file=fopen("final_LL","r");
        if (LL_file==NULL)
        {
                printf("fail to open LL_file\n");
                exit(1);
        }
        for (i=0; i<nsample; i++) 
        {
             for (j=0; j<nparam; j++)
             {                                            
                  fscanf(LL_file,"%d",&LL[i][j]);                
             }
        } 
        fclose(LL_file);
printf("2");
        //read objectives "obj"
        FILE*obj_file;
        char objective_filename[50];
        sprintf(objective_filename,"%s",argv[1]);
        obj_file=fopen(objective_filename,"r");
        if (obj_file==NULL)
        {
                printf("fail to open %s\n",objective_filename);
                exit(1);
        }
        for (i=0; i<nsample; i++) 
        {
             fscanf(obj_file,"%f",&obj[i]);                              
        } 
        fclose(obj_file);
	
printf("read file successfully");
		    
	//calculate F values according ANOVA
	ANOVAA(obj,LL,nsample,nparam,F,R2);

	//bootstrap ANOVA
	Bt_ANOVA(nsample,nparam,np,nrspl,LL,obj,FCI);

	//output results
        FILE*anova_results;
        anova_results=fopen("anovaresults01","w"); 
        if (anova_results==NULL)
        {
            printf("fail to open anovaresults\n");
            exit(1);
        }

            fprintf(anova_results,"%s","Main Effects\n");	
	    for (i=0; i<nparam; i++) 
            {
		 fprintf(anova_results,"%s%d%s\t%10.4f%s%6.4f%s\n","Parameter",i,":",F[i],"[",FCI[i],"]");
	    }

            fprintf(anova_results,"%s","Interaction Effects\n");	
	    for (j=0;j<nparam-1;j++)
            {
		 for (k=j+1;k<nparam;k++)
                 {
                      fprintf(anova_results,"%s%d%s%d%s\t%10.4f%s%6.4f%s\n","Parameters",j,"&",k,":",F[i],"[",FCI[i],"]");
		      i++;
		 }
	    }
                 
            fprintf(anova_results,"%s%6.4f\n","R2-1:",R2[0]);
            fprintf(anova_results,"%s%6.4f\n","R2-2:",R2[1]);
            printf("output results in anova_results\n");   
	//free objects
	FreeFloatVector(F);
	FreeFloatVector(FCI);
	FreeFloatVector(R2);
	FreeFloatVector(obj);
	FreeDoubleMatrix(par,nsample);
	FreeIntMatrix(LL,nsample);
	FreeDoubleVector(lb);
	FreeDoubleVector(ub);
        printf("End of Anova-DHSVM run.\n");

	return 0;
}
