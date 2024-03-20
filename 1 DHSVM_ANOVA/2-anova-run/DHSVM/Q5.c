#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define NSTEP 7670
#define NSAMPLE 12916


void del_date_sim(FILE *fp)
{
        fgetc(fp);// 0
        fgetc(fp);// 1
        fgetc(fp);// .
        fgetc(fp);// 0
        fgetc(fp);// 2
        fgetc(fp);// .
        fgetc(fp);// 2
        fgetc(fp);// 0
        fgetc(fp);// 0
        fgetc(fp);// 0
        fgetc(fp);// -
        fgetc(fp);// 0
        fgetc(fp);// 0
        fgetc(fp);// :
        fgetc(fp);// 0
        fgetc(fp);// 0
        fgetc(fp);// :
        fgetc(fp);// 0
        fgetc(fp);// 0
}
//for DHSVM to loadinput observed outlet discharge and loadoutput simulated outlet discharge
void LoadInput(double Obs[NSTEP])
{
        int i = 0;
        FILE *fp;
        fp=fopen("obsfile.txt","r");//obsfile includes observed outlet discharge
                                                                     //(the first date:1962/01/02,because discharge 
                                                                     //on 01.01.1962 had not been calculated in simfile)
        if (fp==NULL)
        {
                printf("can not open obsfile.txt.\n");
                exit(1);
        }

        for(i=0;i<NSTEP;i++)
        {
                  //del_date_obs(fp);
                  fscanf(fp,"%lf",&Obs[i]);
                  fgetc(fp);//read the end of line, EOF!!
                  //printf("%f\n", Obs[i]);
        }
        fclose(fp);
}

void LoadOutput(double Sim[NSTEP], char filename[200])
{
       char string[50];
       int j = 0;
       FILE *fp;
       fp = fopen(filename,"r");
       if(fp==NULL)
       {
              printf("can not open streamflow.only.\n");
              exit(1);
       }

       fgets(string,50,fp);//delete "DATE OUTLET "
       fgets(string,50,fp);//delete "01.01.1962-00:00:00 "
       for(j=0;j<NSTEP;j++)
       {
                 del_date_sim(fp);
                 fscanf(fp,"%lf",&Sim[j]);
                 Sim[j]=Sim[j]/3600/24;// conversion of unit:m^3/day to m^3/s (same as obsfile).
                 fgetc(fp);//read the end of line, EOF!!!
                 fgetc(fp);
                 //printf("%f\n", Sim[j]);
       }
       fclose(fp);
}

int main(int argc, char* argv[]) 
{
	int i, j, k;
        double sum_obs, average_obs;
        double obj[NSAMPLE], obj1[NSAMPLE], obj2[NSAMPLE];

        if (argc != 2)
        {
            fprintf(stderr, "\nUsage: ./NS output_filename\n\n");
            exit(EXIT_FAILURE);
        }


        //Load Obs_file: observed outlet discharge.
        static double Obs[NSTEP];
        LoadInput(Obs);
   
        for (k=0; k<NSAMPLE; k++) 
        { 
             //Load Sim_file: simulated outlet discharge streamflow.Only (yield by DHSVM).
	     double Sim[NSTEP];
	     char filename[300];
		int NSTEPlow=0;
	     sprintf(filename,"./output/Streamflow.Only.[%d]", k);
             LoadOutput(Sim, filename);  

             //calculate obj: NS.
             sum_obs = 0;
             average_obs = 0;
             obj[k] = 0;
             obj1[k] = 0;
             obj2[k] = 0;
		
             for (j=0; j<NSTEP; j++) 
             { 
			
			NSTEP;
		  sum_obs = sum_obs + Obs[j];
             }	
             average_obs = sum_obs / NSTEP;
             	
             for (j=30; j<NSTEP; j++)
             { 
                 if (Obs[j] >= 140) 
                 {        
                     obj1[k] = obj1[k] + Sim[j];
                 }
             }	
             obj[k] = obj1[k];
         }

         FILE *NSfile;
         char obj_filename[100];
         sprintf(obj_filename,"%s",argv[1]);
         NSfile = fopen(obj_filename,"w");
         if (NSfile==NULL)
         {
             printf("fail to open %s.\n",obj_filename);
             exit(1);
         }
         for (k=0; k<NSAMPLE; k++) 
         {
             fprintf(NSfile,"%lf\n",obj[k]);
         }
         fclose(NSfile);

	return 0;
}

