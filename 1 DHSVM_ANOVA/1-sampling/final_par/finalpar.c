/*In initial sample file(par), Soil Porosity, Field Capacity, Wilting Point may mismatch. This will result in the model DHSVM cannot run, therefore, we should delete this samples.the corresponding parameter level matrix in LL file also should be deleted.*/

#include <stdio.h>
#include <stdlib.h>
#define NSAMPLE 20000 //number of generated parameter sets
#define NPARAM 104 // number of parameters

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

int main()
{
    double** par=DoubleMatrix(NSAMPLE,NPARAM); //parameter matrix
    double** finalpar=DoubleMatrix(NSAMPLE,NPARAM); //final parameter matrix
    int** LL=IntMatrix(NSAMPLE,NPARAM); //parameter Level matrix
    int** finalLL=IntMatrix(NSAMPLE,NPARAM); //final parameter Level matrix
    int i,j,n=0,m=0;
    double porosity[4],field_capacity[4],wilting_point[4];
    
    //read parameter sample "par"
    FILE*parfp;
    parfp=fopen("par","r");
    if (parfp==NULL)
    {
        printf("fail to open par file\n");
        exit(1);
    }
    for (i=0; i<NSAMPLE; i++) 
    {
         for (j=0; j<NPARAM; j++)
         {                                            
              fscanf(parfp,"%lf",&par[i][j]);                
         }
    } 
    fclose(parfp);

    //read parameter level matrix "LL"
    FILE*LLfp;
    LLfp=fopen("LL_file","r");
    if (LLfp==NULL)
    {
        printf("fail to open LL file\n");
        exit(1);
    }
    for (i=0; i<NSAMPLE; i++) 
    {
         for (j=0; j<NPARAM; j++)
         {                                            
              fscanf(LLfp,"%d",&LL[i][j]);                
         }
    } 
    fclose(LLfp);

    //delete samples which is not satisfied for DHSVM.
    for(i=0; i<NSAMPLE; i++)
    { 
        porosity[0]=par[i][4]; //clay
        field_capacity[0]=par[i][7];
        wilting_point[0]=par[i][8];

        porosity[1]=par[i][17]; //SILTY LOAM
        field_capacity[1]=par[i][20];
        wilting_point[1]=par[i][21];
 
        porosity[2] = par[i][30]; //LOAM
        field_capacity[2] = par[i][33];
        wilting_point[2] = par[i][34];

        porosity[3] = par[i][43]; //SANDY CLAY LOAM
        field_capacity[3] = par[i][46];
        wilting_point[3] = par[i][47];

        if((porosity[0]>field_capacity[0])&&(field_capacity[0]>wilting_point[0])
           &&(porosity[1]>field_capacity[1])&&(field_capacity[1]>wilting_point[1])
            &&(porosity[2]>field_capacity[2])&&(field_capacity[2]>wilting_point[2])
            &&(porosity[3]>field_capacity[3])&&(field_capacity[3]>wilting_point[3]))
        {
           for(j=0; j<NPARAM; j++)
           {
               finalpar[n][j]=par[i][j];
               finalLL[n][j]=LL[i][j];
           }
           n++; 
           //printf("n=%d\n",n);          
        }  
        else
        {
            m++;
            printf("m=%d\n",m);            
        }                                      
    }

    //write out the final parameter matrix in "final_par"
    FILE*finalparfp;
    finalparfp=fopen("final_par_no_runnumber","w");
    if (finalparfp==NULL)
    {
        printf("fail to open final_par\n");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        for(j=0; j<NPARAM; j++)
        {
            if (j==NPARAM-1)
            {
                fprintf(finalparfp,"%lf\n",finalpar[i][j]);//if j=nparam-1, newline start. i is RunNumber(new added).
            }   
            else   
                fprintf(finalparfp,"%lf\t",finalpar[i][j]);
        }
    }
    fclose(finalparfp);

    //write out the final parameter matrix in "final_LL"
    FILE*finalLLfp;
    finalLLfp=fopen("final_LL","w");
    if (finalLLfp==NULL)
    {
        printf("fail to open final_LL\n");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        for(j=0; j<NPARAM; j++)
        {
            if (j==NPARAM-1)
            {
                fprintf(finalLLfp,"%d\n",finalLL[i][j]);//if j=nparam-1, newline start. 
            }   
            else   
                fprintf(finalLLfp,"%d\t",finalLL[i][j]);
        }
    }
    fclose(finalLLfp);

    FreeDoubleMatrix(par,NSAMPLE); 
    FreeDoubleMatrix(finalpar,NSAMPLE); 
    FreeIntMatrix(LL,NSAMPLE);
    FreeIntMatrix(finalLL,NSAMPLE);
    return 0;
}
