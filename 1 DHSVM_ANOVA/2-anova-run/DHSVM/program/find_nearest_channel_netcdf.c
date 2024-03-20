/*
 * Program:     find_nearest_channel_netcdf.c 
 * USAGE:       find the neareste channel for each cell
 *
 * AUTHOR:       Ning Sun
 * E-MAIL:       ning@hydro.washington.edu
 * ORIG-DATE:    Feb-2013
 * Last Change:  
 *               
 * Comment:
 * $Id: find_nearest_channel_netcdf.c,v 3.1.1 2013/02/4 Ning Exp $  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fifoNetCDF.h"
#include "sizeofNetCDF.h"
#include "data.h"
#include "fileio.h"
#include "init.h"

/* This code is used to find out the location of the nearest channel to ALL in-basin cells (mask > 0)
   
   ************ One potential issue with this code is that if no channel is found in 100 cell-moves or loops,
    the location of the nearest channel is wherevere at the time of the bail ************

   ****Usage: nrows ncols binary_flowd_file binary_mask_file stream_map_file n_header_map_file
   Note that n_header_map_file = the number of header lines in the stream map file ********

   Algorithms used to find the nearest channel:
1) Recast the flow direction from 1~128
   |-----|-----|-----|      
   | 32  | 64  | 128 |      The central cell searches its nearest channel. The search radius/path coincides with
   |-----|-----|-----|      the flow direction.The searching radius < 100 grid cells. 
   | 16  |     |  1  |       
   |-----|-----|-----|      
   |  8	 |  4  |  2  |
   |-----|-----|-----|

   The flow direction is recasted to 0~7:
   |-----|-----|-----|      
   |  5  |  6  |  7  |      When recasted flowdir =0, the cell moves to the parallel right cell, X axis +1 cell, Y axis + 0 cell
   |-----|-----|-----|      Likewise, when flowdir = 7, X move +1 cell, Y move -1 cell;
   |  4  |     |  0  |      
   |-----|-----|-----|      **(The confusion could be the increasing y is to the south **
   |  3	 |  2  |  1  |        
   |-----|-----|-----|      xneighbor[16] = {1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1}
                            yneighbor[16] = {0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1}

Code Credit:   Original code believed to be written by Pascal Storck, 3/13/2001 for the PRISM project Modifed to deal with 
imperfect basin masks (in coastal areas) and 8-direction flow by Matthew Wiley 12/15/2004.
*********************************************************************************************************************************/

int GetNumber(char *numberStr);
float GetFloat(char *numberStr);
int CopyDouble(double *Value, char *Str, const int NValues);

/* argc stands for "argument count"; argc contains the number of arguments passed to the program. 
   The name of the variable argv stands for "argument vector". argv is a one-dimensional array of strings. 
   Each string is one of the arguments that was passed to the program.

   The first argument (argv[0]) is the name by which the program was called, in this case gcc. 
   This argv[0] plus the 6 argument vectors make 7 as the argc.*/

int main(int argc, char **argv)
{
  FILE *channelfile, *outfile;
  int  ncols,nrows;
  float dx;
  int  NElements;
  unsigned char **mask,**flowd;
  unsigned char **has_channel;
  float *tempf;
  unsigned char *tempuc, *tempuc2;
  int y,x;
  int i,nskip;
  char line[255],flowdname[255],maskname[255],mapname[255], outputpath[255];
  int miny,minx;
  int icol,irow;
  int dist,mindist,sx,sy,mx,my,tempy,tempx,tx,ty;
  int gotchannel;
  int xneighbor[16] = {1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1};
  int yneighbor[16] = {0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1};
  int err = 1;
  MAPSIZE Map;
  MAPDUMP DMap;
  char VarName[255];
  int flag;

  /* Note that the arrays are read in as if the northwest corner is the origin
  increasing x is to the east, increasing y is to the south */

  if(argc < 11) {
    printf("usage: <nrows> <ncols> <cell_size> <Xorig> <Yorig> <flowd_file> \n"); 
	printf("usage: <mask_file> <stream_map_file> <n_header>\n");
    printf("where: flowd_file is a netCDF flow direction file in the same format as mask file\n");
    printf("       make sure that the flowd_file is free of sinks, etc\n");
    printf("       flowdirection is assumed to be from ARC-INFO, i.e. 1 to 128\n");
    printf("       netCDF mask_file and stream_map_file are the DHSVM specific \n");
    printf("       input files for mask and stream_map file, respectively \n");
    printf("       n_header are the number of header lines in the stream map file\n");
    printf("       enter 0 if there are no header lines, i.e. lines starting with #\n");
    printf("       caution: make sure you are referring to the map file not the network file\n");
    exit(-1);
  }

  nrows = GetNumber(argv[1]);
  ncols = GetNumber(argv[2]);
  dx = GetFloat(argv[3]); /* the cellsize of the dem */   
  
  /* extreme west coordinatee */
  if (!(CopyDouble(&(Map.Xorig), argv[4], 1)))
	  exit (-1);;
  /* exterme north coordinate */
  if (!(CopyDouble(&(Map.Yorig), argv[5], 1)))
	  exit (-1);;
  strcpy(flowdname,argv[6]);
  strcpy(maskname,argv[7]);
  strcpy(mapname,argv[8]);
  strcpy(outputpath, argv[9]);
  nskip = GetNumber(argv[10]);


  Map.X = 0;
  Map.Y = 0;
  Map.OffsetX = 0;
  Map.OffsetY = 0;
  Map.NX = ncols;
  Map.NY = nrows;
  Map.DX = dx;
  Map.DY = dx;
  Map.DXY = (float) sqrt(Map.DX * Map.DX + Map.DY * Map.DY);

  /**************************** Open the files ***********************************/
  /*  Note that in order to open a file as a binary file,  
      a "b" character has to be included in the mode string.*/
  if (!(channelfile = fopen(mapname, "rt"))){
    printf("input file not opened \n");
    exit(-1);
  }
  if (!(outfile = fopen(outputpath, "wt"))){
    printf("output file not opened \n");
    exit(-1);
  }
  printf("opened all output files \n");

  /*****Allocate memories for the matrix (row, column)that stores input data*****/

  /*****Allocate memories for the matrix (row, column)that stores input data*****/

   if (!((flowd) = (unsigned char**) calloc(nrows, sizeof(unsigned char*))))
	   exit(-1); 
   for (y = 0; y < nrows; y++) {
	   if (!((flowd)[y] = (unsigned char*) calloc(ncols, sizeof(unsigned char))))
		   exit(-1);
   }
   if (!((mask) = (unsigned char**) calloc(nrows, sizeof(unsigned char*))))
	   exit(-1);
   for (y = 0; y < nrows; y++) {
	   if (!((mask)[y] = (unsigned char*) calloc(ncols, sizeof(unsigned char))))
       exit(-1);
   }
   if (!((has_channel) = (unsigned char**) calloc(nrows, sizeof(unsigned char*))))
	   exit(-1);
   for (y = 0; y < nrows; y++) {
	   if (!((has_channel)[y] = (unsigned char*) calloc(ncols, sizeof(unsigned char))))
      exit(-1);
   }

   /************************ allocate some memory for the reader 1-d arrays *******************/
   if (!(tempf = (float *) calloc(nrows * ncols, 
				sizeof(float)))) {
     printf("failed to allocate memory \n");
     exit(-1);
   }
   if (!(tempuc = (unsigned char *) calloc(nrows * ncols, 
				sizeof(unsigned char)))) {
     printf("failed to allocate memory \n");
     exit(-1);
   }
   if (!(tempuc2 = (unsigned char *) calloc(nrows * ncols, 
				sizeof(unsigned char)))) {
     printf("failed to allocate memory \n");
     exit(-1);
   }
  printf("assigned all the memory \n");

 /****************Read in the flow direction data and assign to maps********************/
  strcpy(VarName, "Flow.Dir");
  flag = Read2DMatrixNetCDF(flowdname, tempuc, NC_BYTE, Map.NY, Map.NX, 0,
	       VarName, 0);  
  if (flag == 0){
	  for (y = 0, i = 0; y < Map.NY; y++) {
		  for (x = 0; x < Map.NX; x++, i++) {
			  flowd[y][x] = tempuc[i]; }
	  }
  }
  else if (flag == 1){
	  for (y = Map.NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map.NX; x++, i++) {
			  flowd[y][x] = tempuc[i]; 
		  }
	  }
  }
  else exit (-1);

  printf("recasting flowdirections: from 1-128 to 0-7\n");

  /***********************Recasting Flowdirections: from 1-128 to 0-7*************************/
  for (y = 0; y < nrows; y++) {
	  for (x = 0; x < ncols; x++) {  
		  if(flowd[y][x]==1)  
			  flowd[y][x]=0;
		  else if(flowd[y][x]==2)
			  flowd[y][x]=1;
		  else if(flowd[y][x]==4)
			  flowd[y][x]=2;
		  else if(flowd[y][x]==8)
			  flowd[y][x]=3;
		  else if(flowd[y][x]==16)
			  flowd[y][x]=4;
		  else if(flowd[y][x]==32)
			  flowd[y][x]=5;
		  else if(flowd[y][x]==64)
			  flowd[y][x]=6;
		  else if(flowd[y][x]==128)
			  flowd[y][x]=7;
		  else {
			  printf("pixel with undefined flow direction encountered at [%d][%d]:%d %d\n",
              y,x,tempuc[y*ncols+x],flowd[y][x]);
			  exit(-1);
		  }
	  }
  }
  printf("got the flowd \n"); 
  free(tempuc);

  /****************Read in the mask file and assign to maps********************/
  strcpy(VarName, "Basin.Mask");
  flag = Read2DMatrixNetCDF(maskname, tempuc2, NC_BYTE, Map.NY, Map.NX, 0,
	       VarName, 0);  
  if (flag == 0){
	  for (y = 0, i = 0; y < Map.NY; y++) {
		  for (x = 0; x < Map.NX; x++, i++) {
			  mask[y][x] = tempuc2[i]; 
			  has_channel[y][x]=0;
		  }
	  }
  }
  else if (flag == 1){
	  for (y = Map.NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map.NX; x++, i++) {
			  mask[y][x] = tempuc2[i]; 
			  has_channel[y][x]=0;
		  }
	  }
  }
  else exit (-1);
  free(tempuc2);
  printf("got the mask \n"); 

  /****************Read in the channel (stream map) file and assign to maps********************

  Stream file (or channel file as indicated in this code) includes:
  (1) no. of row
  (2) no. of column
  (3) Channel ID
  (4) length: straight-line length of the channel segment lying within the cell in meters (float)	
  (2) Depth: stream bank height in meters (float)	
  (3) Width: effective stream channel width in meters (float)
  ********************************************************************************************/
  
  /* This short script reads the specified lines (=nskip)of the file hearder, and
     the pointer moves to the first line below the header lines at the end ******************/
  for(i=0;i<nskip;i++)
  {
	  fgets(line,255,channelfile);
  }
  /******************************************************************************************/
  while((fgets(line,255,channelfile)!=NULL))
  {
	  /*** sscanf ( ) reads formatted data: in this instance,
	  it reads the first 2 columns of channelfile.**********/
	  sscanf(line,"%d %d ", &icol,&irow); 
      has_channel[irow][icol]=1;
      //printf("%d %d \n",icol,irow);
  } 
 
  /* Trace each pixel in the masked area to the nearest downslope pixel */
  printf("looking for channels \n");
  for (y = 0; y < nrows; y++) 
  {
	  // printf("row %d \n",y);
	  for (x = 0; x < ncols; x++) 
	  {
		  if (mask[y][x]>0) 
		  {
			  /* remember the starting point and reset the total distance (in # of cells of moves)*/
			  sx=x;
		      sy=y;
			  mx=x;
		      my=y;
			  dist = 0;
			  gotchannel=0;
			  gotchannel=has_channel[sy][sx];
			  /* provide a switch to bail if we haven't escaped in 100 moves, uses wherevere it was at the time of the bail */
			  while(dist < 100.0 && gotchannel == 0) 
			  {
				  /*make the move and check to see if it has a channel*/   
				  tempx=xneighbor[flowd[my][mx]];
				  tempy=yneighbor[flowd[my][mx]];
				  tx = mx;
				  ty = my;
				  mx+=tempx;
				  my+=tempy;
				  err = 1;
				  while (mask[my][mx] == 0 ) 
				  {
					  //printf("Routing out of mask! Changing flowd[%d][%d] from %d to %d\n", ty,tx,flowd[ty][tx], flowd[ty][tx]+ err);
					  mx = tx;
					  my = ty;
				      tempx=xneighbor[flowd[ty][tx]+ err];
					  tempy=yneighbor[flowd[ty][tx]+ err];
					  err++;
					  mx+=tempx;
					  my+=tempy;
				  }
				  gotchannel=has_channel[my][mx];
				  dist+=1.0;
				  //printf("start cell(%d, %d) to cell(%d, %d) with flow dir = %d \n",sy,sx,my,mx,flowd[my][mx]);
			  }
			  printf("nearest channel to cell(%d, %d) is cell(%d, %d) \n",y,x,my,mx);
			  fprintf(outfile,"%d %d %d %d \n",sy,sx,my,mx);
		  } 
	  }
  }
  return EXIT_SUCCESS;
}

/*****************************************************************************
  GetNumber()
*****************************************************************************/
int GetNumber(char *numberStr) 
{
char *endPtr;
int number = 0;

/* strtol() converts string to long integer */
number = (int) strtol(numberStr, &endPtr, 0);
if (*endPtr != '\0')
     exit(-1);

     return number;
}
/*****************************************************************************
  GetFloat()
*****************************************************************************/
float GetFloat(char *numberStr) 
{
  char *endPtr;
  float number = 0;

  number = (float) strtod(numberStr, &endPtr);
  if (*endPtr != '\0'){
    printf("problem extracting float from %s \n",numberStr);
    exit(-1);
  }

  return number;
}
/*****************************************************************************
  CopyDouble()
*****************************************************************************/
int CopyDouble(double *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
