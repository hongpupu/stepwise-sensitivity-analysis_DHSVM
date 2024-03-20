#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#define NSAMPLE 12916 //number of final parameter sets
#define NPARAM 105 //NPARAM+RunNumber=104+1=105
int min(int x,int y)
{
   if (x>=y)
     return 1;
   else
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

int main(int argc,char *argv[])
{
    double** par=DoubleMatrix(NSAMPLE,NPARAM); //parameter matrix
    int i,j,k;
    int myid,numprocs;
    double buffer[NPARAM]={0};
    
    //read parameter sample "par"
    FILE*parfp;
    parfp=fopen("final_par","r");
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


    MPI_Init(&argc,&argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    int send_cycle_number=(int)(NSAMPLE/(numprocs-1));
    int remainder=NSAMPLE%(numprocs-1);
    if (myid==0)//send par from master processor.
    {
        printf("send_cycle_number= %d,remainder=%d\n",send_cycle_number,remainder);
        for (k=0; k<send_cycle_number; k++)//send_cycle_number
        {
             for (i=0; i<numprocs-1; i++)
             {
                  for (j=0; j<NPARAM; j++)
                  {
                       buffer[j]=par[k*(numprocs-1)+i][j];
                  }
                  MPI_Send(buffer,NPARAM,MPI_DOUBLE,i+1,k*(numprocs-1)+i,MPI_COMM_WORLD);
             }
        }
        for (i=0; i<remainder; i++)
        {
             for (j=0; j<NPARAM; j++)
             {
                  buffer[j]=par[send_cycle_number*(numprocs-1)+i][j];
             }
             MPI_Send(buffer,NPARAM,MPI_DOUBLE,i+1,send_cycle_number*(numprocs-1)+i,MPI_COMM_WORLD);
        }
    }
    else
    {
        //receive message.
        char Filename[100];//tem_file
        FILE*temfp[numprocs];//the temporary file is used to deposit a sample which will be runned by DHSVM soon.
        char name_for_system[400];

        for (i=0; i<send_cycle_number+min(remainder,myid); i++)//if remainder-1>=myid,plus 1; else 0.
        {
             MPI_Recv(buffer,NPARAM,MPI_DOUBLE,0,i*(numprocs-1)+myid-1,MPI_COMM_WORLD,&status);

             sprintf(Filename,"tem_file[%d]",myid);
             printf("Filename = %s\n",Filename);
             temfp[myid]=fopen(Filename,"w");
             if (temfp[myid]==NULL)
             {
                 printf("fail to open tem_file[%d]\n",myid);
                 exit(1);
             }
             else
                 printf("success to open tem_file[%d]\n",myid);
             for (j=0; j<NPARAM; j++)
             {
                  fprintf(temfp[myid],"%lf\t",buffer[j]);
             }
             fclose(temfp[myid]);
             
             sprintf(name_for_system,"./DHSVM/sourcecode/DHSVM3.2 ./DHSVM/config/Shaduan_modified2.txt %d",myid);
             int b = system(name_for_system); //run another exe

             if (b!=0)
             {
                 printf("fail to run DHSVM\n");
                 exit(1);
             }
             else
                 printf("success to run DHSVM\n");
        }
    }
    
    MPI_Finalize();
    FreeDoubleMatrix(par,NSAMPLE);
    
    return 0;
}
