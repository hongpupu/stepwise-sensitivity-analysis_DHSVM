/*
 * SUMMARY:      InitMetSources.c - Initialize met sources for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Sat May 31 00:31:29 1997 by  <nijssen@u.washington.edu>
 * DESCRIPTION:  Initialize meteorology options for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    InitMetSources()
 *               InitMM5()
 *               InitRadar()
 *               InitStations()
 *               InitWindModel()
 *
 * COMMENTS:
 * $Id: InitMetSources.c,v3.1.2 2013/10/28 ning Exp $
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"
                          
 /*******************************************************************************
   Function name: InitMetSources()

   Purpose      : Initialize and configure the model to process meteorological
                  data from various different sources
                  Processes the following section in the input file:
                  [METEOROLOGY]

   Required     :
     LISTPTR Input           - Linked list with input strings
     OPTIONSTRUCT *Options   - Structure with different program options
     MAPSIZE *Map            - Coverage and resolution of model area
     int NSoilLayers         - Number of soil layers
     TIMESTRUCT *Time        - Begin and end times, model timestep
     INPUTFILES *InFiles     - Various input filenames
     int *NStats             - Number of meteorological stations
     METLOCATION **Stat      - Station information
     MAPSIZE *Radar          - Areal extent of radar files

   Returns      : void

   Modifies     : Members of Time, InFiles, RadClass, and Radar

   Comments     :
 *******************************************************************************/
void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  TOPOPIX **TopoMap, int NSoilLayers, TIMESTRUCT *Time, INPUTFILES *InFiles,
  int *NStats, METLOCATION **Stat, MAPSIZE *Radar, MAPSIZE *MM5Map,
  GRID *Grid)
{
  const char *Routine = "InitMetSources";

  if (Options->Outside == TRUE && Options->MM5 == FALSE) {
    printf("\nAll met stations in list will be included \n");
    if (Options->Prism == TRUE) {
      printf("WARNING: PRISM Option is also on\n");
      printf("Make sure file .prism files exist\n\n");
    }
  }
  /* The MM5 option overrides all other options, so check that one first */
  if (Options->MM5 == TRUE) {
    InitMM5(Input, NSoilLayers, Time, InFiles, Options, MM5Map, Map);
  }

  /* Use gridded met forcing data */
  if (Options->GRIDMET == TRUE)
    InitGridMet(Options, Input, Map, TopoMap, Grid, Stat, NStats);

  /* otherwise, check and initialize the other options */
  if (Options->QPF == TRUE || (Options->MM5 == FALSE && Options->GRIDMET == FALSE))
    InitStations(Input, Map, Time->NDaySteps, Options, NStats, Stat);

  if (Options->QPF == TRUE || Options->MM5 == FALSE) {
    if (Options->PrecipType == RADAR)
      InitRadar(Input, Map, Time, InFiles, Radar);
    if (Options->WindSource == MODEL)
      InitWindModel(Input, InFiles, *NStats, *Stat);
    if (Options->PrecipLapse == MAP) {
      if (*NStats > 1)
        ReportError((char *)Routine, 54);
      InitPrecipLapse(Input, InFiles);
    }
  }
}

/*******************************************************************************
  Function name: InitStations()

  Purpose      : Read the station information from the options file.  This
                 information is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    MAPSIZE *Map        - Information about the basin area
    int NDaysSteps      - Number of time steps in a day
    int *NStats         - Number of met stations
    METLOCATION **Stat  - Information about each met station

  Returns      : void

  Modifies     : NStats, Stat and members of Stat

  Comments     :
*****************************************************************************/
void InitStations(LISTPTR Input, MAPSIZE *Map, int NDaySteps,
  OPTIONSTRUCT *Options, int *NStats, METLOCATION **Stat)
{
  char *Routine = "InitStations";
  int i;
  int j;
  int k;
  char tempfilename[BUFSIZE + 1];
  char KeyName[station_file + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "STATION NAME",
    "NORTH COORDINATE",
    "EAST COORDINATE",
    "ELEVATION",
    "STATION FILE"
  };
  char *SectionName = "METEOROLOGY";
  char VarStr[station_file + 1][BUFSIZE + 1];
  float East;
  float North;
  FILE *PrismStatFile;

  /* Get the number of different stations */
  GetInitString(SectionName, "NUMBER OF STATIONS", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(NStats, VarStr[0], 1))
    ReportError("NUMBER OF STATIONS", 51);

  if (*NStats <= 0)
    ReportError("Input Options File", 6);

  printf("\nEvaluating %d Met stations for inclusion\n", *NStats);

  /* Allocate memory for the stations */
  if (!(*Stat = (METLOCATION *)calloc(*NStats, sizeof(METLOCATION))))
    ReportError(Routine, 1);

  /* Read key-entry pairs for each station from the input file */
  /* for each potential station, up to NStats, read in the data and */
  /* determine if it is in the current model bounding box */
  /* If it is then put it into memory, otherwise, forget about it */
  /* unless Outside option is TRUE, then include it anyway */
  /* use temp counter k to track number of valid stations */
  k = 0;
  for (i = 0; i < *NStats; i++) {
    for (j = 0; j <= station_file; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
        (unsigned long)BUFSIZE, Input);
    }

    /* Assign the entries to the variables */
    if (IsEmptyStr(VarStr[station_name]))
      ReportError(KeyName[number_of_maps], 51);
    strcpy((*Stat)[k].Name, VarStr[station_name]);

    if (!CopyFloat(&North, VarStr[station_north], 1))
      ReportError(KeyName[station_north], 51);

    if (!CopyFloat(&East, VarStr[station_east], 1))
      ReportError(KeyName[station_east], 51);

    (*Stat)[k].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
    (*Stat)[k].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);

    if (!CopyFloat(&((*Stat)[k].Elev), VarStr[station_elev], 1))
      ReportError(KeyName[station_elev], 51);

    if (IsEmptyStr(VarStr[station_file]))
      ReportError(KeyName[station_file], 51);
    strcpy((*Stat)[k].MetFile.FileName, VarStr[station_file]);

    OpenFile(&((*Stat)[k].MetFile.FilePtr), (*Stat)[k].MetFile.FileName, "r", FALSE);

    /* check to see if the stations are inside the bounding box */
    if (((*Stat)[k].Loc.N >= Map->NY || (*Stat)[k].Loc.N < 0 ||
      (*Stat)[k].Loc.E >= Map->NX || (*Stat)[k].Loc.E < 0)
      && Options->Outside == FALSE){
      k = k;
	  //printf("Station %d outside bounding box: %s ignored\n", i + 1, (*Stat)[k].Name);
	}
    else
      k = k + 1;
  }
  if (Options->Outside == FALSE)
    printf("Final number of stations in bounding box is %d \n\n", k);
  else
    printf("Forced to include all %d stations \n", k);
  *NStats = k;

  if (Options->Outside == TRUE && Options->Prism == TRUE) {

    for (i = 0; i < *NStats; i++) {
      sprintf(tempfilename, "%s.prism", (*Stat)[i].MetFile.FileName);
      /* Options->PrismDataExt); */
      OpenFile(&PrismStatFile, tempfilename, "rt", FALSE);
      for (k = 0; k < 12; k++) {
        fscanf(PrismStatFile, "%f ", &(*Stat)[i].PrismPrecip[k]);
      }
      fclose(PrismStatFile);
    }
  }
}

/*******************************************************************************
Function name: InitGridMet()

Purpose      : Read the gridded met file.  This information
is in the [METEOROLOGY] section
*****************************************************************************/
void InitGridMet(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  TOPOPIX **TopoMap, GRID *Grid, METLOCATION **Stat, int *NStats)
{
  char *Routine = "InitGridMet";
  char KeyName[BUFSIZE + 1];
  char VarStr[BUFSIZE + 1];
  int i, k, m;
  float lat, lon, North, East;
  char tempfilename[BUFSIZE + 1];
  FILE *PrismStatFile;
  char junk[BUFSIZE + 1], infileformat[BUFSIZE + 1];
  DIR *dir;                   
  struct dirent *ent;

  STRINIENTRY StrEnv[] = {
    { "METEOROLOGY", "EXTREME NORTH LAT", "", "" },
    { "METEOROLOGY", "EXTREME SOUTH LAT", "", "" },
    { "METEOROLOGY", "EXTREME EAST LON", "", "" },
    { "METEOROLOGY", "EXTREME WEST LON", "", "" },
    { "METEOROLOGY", "NUMBER OF GRIDS", "", "" },
	{ "METEOROLOGY", "GRID_DECIMAL", "", "" },
    { "METEOROLOGY", "MET FILE PATH", "", "" },
    { "METEOROLOGY", "FILE PREFIX", "", "" },
    { "METEOROLOGY", "UTM ZONE", "", "" },                                      
    { NULL, NULL, "", NULL },
  };

  /* Allocate memory for the stations */
  if (!(Grid = (GRID *)calloc(1, sizeof(GRID))))
    ReportError(Routine, 1);

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);

  /* if the outside is FALSE, only the grids within the bounding box defined
  by the parameters below are included. Else, the grids within the basin mask
  are include. */
  if (Options->Outside == TRUE) {
    if (!CopyFloat(&(Grid->LatNorth), StrEnv[grid_ext_north].VarStr, 1))
      ReportError(StrEnv[grid_ext_north].KeyName, 51);
    if (!CopyFloat(&(Grid->LatSouth), StrEnv[grid_ext_south].VarStr, 1))
      ReportError(StrEnv[grid_ext_south].KeyName, 51);
    if (!CopyFloat(&(Grid->LonEast), StrEnv[grid_ext_east].VarStr, 1))
      ReportError(StrEnv[grid_ext_east].KeyName, 51);
    if (!CopyFloat(&(Grid->LonWest), StrEnv[grid_ext_west].VarStr, 1))
      ReportError(StrEnv[grid_ext_west].KeyName, 51);
  }
  else {
    Grid->LatNorth = NOT_APPLICABLE;
    Grid->LatSouth = NOT_APPLICABLE;
    Grid->LonEast = NOT_APPLICABLE;
    Grid->LonWest = NOT_APPLICABLE;
  }

  /* estimate of total grid cells for the basin (must > actual used grids */
  if (!CopyInt(&(Grid->NGrids), StrEnv[tot_grid].VarStr, 1))
    ReportError(StrEnv[tot_grid].KeyName, 51);
  
  /* Number of digits after decimal point in forcing file names */
  if (!CopyInt(&(Grid->Decimal), StrEnv[decim].VarStr, 1))
    ReportError(StrEnv[decim].KeyName, 51);
  
  /* path to all grid met files */
  if (IsEmptyStr(StrEnv[grid_met_file].VarStr))
    ReportError(StrEnv[grid_met_file].KeyName, 51);
  strcpy(Grid->filepath, StrEnv[grid_met_file].VarStr);

  /* (same) file prefix of all grid met files */
  if (IsEmptyStr(StrEnv[file_prefix].VarStr))
    ReportError(StrEnv[file_prefix].KeyName, 51);
  strcpy(Grid->fileprefix, StrEnv[file_prefix].VarStr);

  /* UTM zone used as reference for all spatial input files */
  if (!CopyInt(&(Grid->utmzone), StrEnv[utm_zone].VarStr, 1))
    ReportError(StrEnv[utm_zone].KeyName, 51);                                              

  /* Allocate memory for the stations */
  if (!(*Stat = (METLOCATION *)calloc(Grid->NGrids, sizeof(METLOCATION))))
    ReportError(Routine, 1);

  /* define the format of file input */
  sprintf(junk, "%s_%%f_%%f", Grid->fileprefix);
  /* define input file format */
  sprintf(infileformat, "%%s%%s_%%.%if_%%.%if", Grid->Decimal, Grid->Decimal);
  
  printf("\nReading the gridded met files ...\n");
  k = 0;
  m = 0;
  if ((dir = opendir(Grid->filepath)) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir(dir)) != NULL) {                        
      if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, "..")) {
        continue;
      }
      else {
        sscanf(ent->d_name, junk, &lat, &lon);

        /* convert lat and lon to utm */
        deg2utm(lat, lon, &East, &North, Grid->utmzone);		
		sprintf((*Stat)[k].Name, "data_%f_%f\n", lat, lon);
    printf("%f, %f, %f, %f\n", lat, lon, East, North);
		
		(*Stat)[k].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
        (*Stat)[k].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);
        m += 1;

        /* met grids must be with the bounding box of the basin */
		if (((*Stat)[k].Loc.N >= Map->NY || (*Stat)[k].Loc.N < 0 ||
          (*Stat)[k].Loc.E >= Map->NX || (*Stat)[k].Loc.E < 0)) {                                   
          //printf("..... Station %d outside the basin bounding box: %s ignored\n", m, (*Stat)[k].Name);           
		  k = k;
		}
        else {
          /* only include grids within the mask */
          if (Options->Outside == FALSE) {
            if (INBASIN(TopoMap[(*Stat)[k].Loc.N][(*Stat)[k].Loc.E].Mask)) {
                                                                         
			  /* open met data file */
              sprintf((*Stat)[k].MetFile.FileName, infileformat, Grid->filepath, Grid->fileprefix, lat, lon);
                                                                  
              if (!((*Stat)[k].MetFile.FilePtr = fopen((*Stat)[k].MetFile.FileName, "r"))) {
                printf("..... %s doesn't exist\n", (*Stat)[k].MetFile.FileName);
                continue;
              }                                        
              printf("..... Station %d: %s is selected\n", m, (*Stat)[k].Name);            
              k = k + 1;             
            }
          }
          else {
            if (lat <= Grid->LatNorth && lat >= Grid->LatSouth && lon >= Grid->LonWest && lon <= Grid->LonWest) {
              if (!((*Stat)[k].MetFile.FilePtr = fopen((*Stat)[k].MetFile.FileName, "r"))) {
                //printf("..... %s doesn't exist\n", (*Stat)[k].MetFile.FileName);
                continue;
              }                                        
              printf("..... Station %d: %s is selected\n", m, (*Stat)[k].Name);               
              k = k + 1;
            }
          }
        }
      }
    }
    closedir(dir);
  }
  else {
    /* could not open directory */
    ReportError(Grid->filepath, 3);
  }

  /* if no grid is found within the mask, exit with error */
  if (k < 1)
    ReportError(Routine, 69);

  *NStats = k;
  printf("Final number of stations in bounding box is %d \n\n", k);

  if (Options->Outside == TRUE && Options->Prism == TRUE) {
    for (i = 0; i < *NStats; i++) {
      sprintf(tempfilename, "%s.prism", (*Stat)[i].MetFile.FileName);
      /* Options->PrismDataExt); */
      OpenFile(&PrismStatFile, tempfilename, "rt", FALSE);
      for (k = 0; k < 12; k++) {
        fscanf(PrismStatFile, "%f ", &(*Stat)[i].PrismPrecip[k]);
      }
      fclose(PrismStatFile);
    }
  }
}
/*******************************************************************************
  Function name: InitMM5()

  Purpose      : Read the MM5 information the options file.  This information
                 is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input options
    int NSoilLayers     - Number of soil layers
    TIMESTRUCT *Time    - Time information
    INPUTFILES *InFiles - Filenames for various input files

  Returns      : void

  Modifies     : Members of Time and InFiles

  Comments     :
*****************************************************************************/
void InitMM5(LISTPTR Input, int NSoilLayers, TIMESTRUCT *Time,
  INPUTFILES *InFiles, OPTIONSTRUCT *Options, MAPSIZE *MM5Map, MAPSIZE *Map)
{
  DATE Start;
  char *Routine = "InitMM5";
  char KeyName[BUFSIZE + 1];
  char VarStr[BUFSIZE + 1];
  int i;
  STRINIENTRY StrEnv[] = {
    {"METEOROLOGY", "MM5 START", "", ""},
    {"METEOROLOGY", "MM5 TEMPERATURE FILE", "", ""},
    {"METEOROLOGY", "MM5 HUMIDITY FILE", "", ""},
    {"METEOROLOGY", "MM5 WIND SPEED FILE", "", ""},
    {"METEOROLOGY", "MM5 SHORTWAVE FILE", "", ""},
    {"METEOROLOGY", "MM5 LONGWAVE FILE", "", ""},
    {"METEOROLOGY", "MM5 PRECIPITATION FILE", "", ""},
    {"METEOROLOGY", "MM5 TERRAIN FILE", "", ""},
    {"METEOROLOGY", "MM5 TEMP LAPSE FILE", "", "none"},
    {"METEOROLOGY", "MM5 TEMP LAPSE FREQUENCY", "", "single"},
    {"METEOROLOGY", "MM5 ROWS", "", ""},
    {"METEOROLOGY", "MM5 COLS", "", ""},
    {"METEOROLOGY", "MM5 EXTREME NORTH", "", ""},
    {"METEOROLOGY", "MM5 EXTREME WEST", "", ""},
    {"METEOROLOGY", "MM5 DY", "", ""},
    {"METEOROLOGY", "MM5 PRECIPITATION DISTRIBUTION FILE", "", "none"},
    /* can be one of "single", "month", "continuous" */
    {"METEOROLOGY", "MM5 PRECIPITATION DISTRIBUTION FREQUENCY", "", "single"},
    {NULL, NULL, "", NULL},
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);

  /* Assign the entries to the variables */
  if (!SScanDate(StrEnv[MM5_start].VarStr, &Start))
    ReportError(StrEnv[MM5_start].KeyName, 51);

  InitTime(Time, NULL, NULL, NULL, &Start, Time->Dt);

  if (IsEmptyStr(StrEnv[MM5_temperature].VarStr))
    ReportError(StrEnv[MM5_temperature].KeyName, 51);
  strcpy(InFiles->MM5Temp, StrEnv[MM5_temperature].VarStr);
                             
  if (IsEmptyStr(StrEnv[MM5_terrain].VarStr))
    ReportError(StrEnv[MM5_terrain].KeyName, 51);
  strcpy(InFiles->MM5Terrain, StrEnv[MM5_terrain].VarStr);
                                
  if (strncmp(StrEnv[MM5_lapse].VarStr, "none", 4)) {

    strncpy(InFiles->MM5Lapse, StrEnv[MM5_lapse].VarStr, BUFSIZE);

    CopyLCase(VarStr, StrEnv[MM5_lapse_freq].VarStr, BUFSIZE + 1);

    if (!strncmp(VarStr, "single", 6)) {
      InFiles->MM5LapseFreq = FreqSingle;
    } else if (!strncmp(VarStr, "month", 6)) {
      InFiles->MM5LapseFreq = FreqMonth;
    } else if (!strncmp(VarStr, "continuous", 10)) {
      InFiles->MM5LapseFreq = FreqContinous;
    } else {
      ReportError(StrEnv[MM5_lapse_freq].VarStr, 70);
    }
    
  } else {
    strcpy(InFiles->MM5Lapse, "");
    if (Options->TempLapse != CONSTANT) {
      ReportError(StrEnv[MM5_lapse].VarStr, 51);
    }
  }

  if (IsEmptyStr(StrEnv[MM5_lapse].VarStr))
    ReportError(StrEnv[MM5_lapse].KeyName, 51);

  if (IsEmptyStr(StrEnv[MM5_humidity].VarStr))
    ReportError(StrEnv[MM5_humidity].KeyName, 51);
  strcpy(InFiles->MM5Humidity, StrEnv[MM5_humidity].VarStr);
                                 
  if (IsEmptyStr(StrEnv[MM5_wind].VarStr))
    ReportError(StrEnv[MM5_wind].KeyName, 51);
  strcpy(InFiles->MM5Wind, StrEnv[MM5_wind].VarStr);
                             

  if (IsEmptyStr(StrEnv[MM5_shortwave].VarStr))
    ReportError(StrEnv[MM5_shortwave].KeyName, 51);
  strcpy(InFiles->MM5ShortWave, StrEnv[MM5_shortwave].VarStr);
                                  
  if (IsEmptyStr(StrEnv[MM5_longwave].VarStr))
    ReportError(StrEnv[MM5_longwave].KeyName, 51);
  strcpy(InFiles->MM5LongWave, StrEnv[MM5_longwave].VarStr);
                                 
  if (IsEmptyStr(StrEnv[MM5_precip].VarStr))
    ReportError(StrEnv[MM5_precip].KeyName, 51);
  strcpy(InFiles->MM5Precipitation, StrEnv[MM5_precip].VarStr);
                                      
  /* Use PrecipLapseFile to store the name of the MM5 precip
     distribution map file. This avoids wholesale code changes. */
  if (strncmp(StrEnv[MM5_precip_dist].VarStr, "none", 4)) {
    strcpy(InFiles->PrecipLapseFile, StrEnv[MM5_precip_dist].VarStr);

    /* Precipitation distribution can be in several forms. */

    CopyLCase(VarStr, StrEnv[MM5_precip_freq].VarStr, BUFSIZE + 1);
    
    if (!strncmp(VarStr, "single", 6)) {
      InFiles->MM5PrecipDistFreq = FreqSingle;
    } else if (!strncmp(VarStr, "month", 6)) {
      InFiles->MM5PrecipDistFreq = FreqMonth;
    } else if (!strncmp(VarStr, "continuous", 10)) {
      InFiles->MM5PrecipDistFreq = FreqContinous;
    } else {
      ReportError(StrEnv[MM5_precip_freq].VarStr, 70);
    }
  } else {
    strcpy(InFiles->PrecipLapseFile, "");
  }
                                   
  if (Options->HeatFlux == TRUE) {                               
    if (!(InFiles->MM5SoilTemp = (char **)calloc(sizeof(char *), NSoilLayers)))
      ReportError(Routine, 1);

    for (i = 0; i < NSoilLayers; i++) {
      if (!(InFiles->MM5SoilTemp[i] =
          (char *)calloc(sizeof(char), BUFSIZE + 1)))
        ReportError(Routine, 1);
      sprintf(KeyName, "MM5 SOIL TEMPERATURE FILE %d", i);
      GetInitString("METEOROLOGY", KeyName, "", VarStr,
        (unsigned long)BUFSIZE, Input);
      if (IsEmptyStr(VarStr))
        ReportError(KeyName, 51);
      strcpy(InFiles->MM5SoilTemp[i], VarStr);
    }
  }

  if (!CopyDouble(&(MM5Map->Yorig), StrEnv[MM5_ext_north].VarStr, 1))
    ReportError(StrEnv[MM5_ext_north].KeyName, 51);

  if (!CopyDouble(&(MM5Map->Xorig), StrEnv[MM5_ext_west].VarStr, 1))
    ReportError(StrEnv[MM5_ext_west].KeyName, 51);

  if (!CopyInt(&(MM5Map->NY), StrEnv[MM5_rows].VarStr, 1))
    ReportError(StrEnv[MM5_rows].KeyName, 51);

  if (!CopyInt(&(MM5Map->NX), StrEnv[MM5_cols].VarStr, 1))
    ReportError(StrEnv[MM5_cols].KeyName, 51);

  if (!CopyFloat(&(MM5Map->DY), StrEnv[MM5_dy].VarStr, 1))
    ReportError(StrEnv[MM5_dy].KeyName, 51);

  MM5Map->OffsetX = Round(((float)(MM5Map->Xorig - Map->Xorig)) /
    ((float)Map->DX));
  MM5Map->OffsetY = Round(((float)(MM5Map->Yorig - Map->Yorig)) /
    ((float)Map->DY));

  if (MM5Map->OffsetX > 0 || MM5Map->OffsetY < 0)
    ReportError("Input Options File", 31);
                
  printf("MM5 extreme north / south is %f %f \n", MM5Map->Yorig,
    MM5Map->Yorig - MM5Map->NY * MM5Map->DY);
  printf("MM5 extreme west / east is %f %f\n", MM5Map->Xorig,
    MM5Map->Xorig + MM5Map->NX * MM5Map->DY);
  printf("MM5 rows is %d \n", MM5Map->NY);
  printf("MM5 cols is %d \n", MM5Map->NX);
  printf("MM5 dy is %f \n", MM5Map->DY);
  printf("Temperature Map is %s\n", InFiles->MM5Temp);
  printf("Precip Map is %s\n", InFiles->MM5Precipitation);
  if (strlen(InFiles->PrecipLapseFile) > 0) {
    printf("Precip Distribution Map is %s\n", InFiles->PrecipLapseFile);
    switch (InFiles->MM5PrecipDistFreq) {
    case (FreqSingle):
      printf("A single precipititation distribution map is used for the entire simulation\n");
      break;
    case (FreqMonth):
      printf("Monthly precipititation distribution maps are used\n");
      break;
    case (FreqContinous):
      printf("An individual precipititation distribution maps is used for each MM5 time step.\n");
      break;
    default:
      ReportError("InitMM5", 15);
      break;
    }
  } else {
    printf("Precip is distributed evenly within MM5 cell\n");
  }                                                            
  printf("wind Map is %s\n", InFiles->MM5Wind);
  printf("shortwave Map is %s\n", InFiles->MM5ShortWave);
  printf("humidity Map is %s\n", InFiles->MM5Humidity);
                                                         
  if (strlen(InFiles->MM5Lapse) > 0) {
    printf("lapse Map is %s\n", InFiles->MM5Lapse);
  } else {
    printf("temperature lapse rate is constant\n");
  }
  printf("terrain Map is %s\n", InFiles->MM5Terrain);
  printf("MM5 offset x is %d \n", MM5Map->OffsetX);
  printf("MM5 offset y is %d \n", MM5Map->OffsetY);
  printf("dhsvm extreme north / south is %f %f \n", Map->Yorig,
    Map->Yorig - Map->NY * Map->DY);
  printf("dhsvm extreme west / east is %f %f \n", Map->Xorig,
    Map->Xorig + Map->NX * Map->DY);
  printf("fail if %d > %d\n",
    (int)((Map->NY + MM5Map->OffsetY) * Map->DY / MM5Map->DY),
    MM5Map->NY);
  printf("fail if %d > %d\n",
    (int)((Map->NX - MM5Map->OffsetX) * Map->DX / MM5Map->DY),
    MM5Map->NX);
  
  if ((int)((Map->NY + MM5Map->OffsetY) * Map->DY / MM5Map->DY) > MM5Map->NY
    || (int)((Map->NX - MM5Map->OffsetX) * Map->DX / MM5Map->DY) >
    MM5Map->NX)
    ReportError("Input Options File", 31);

}

/*******************************************************************************
  Function name: InitRadar()

  Purpose      : Read the radar information from the options file.  This
                 information is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    MAPSIZE *Map        - Information about basin area
    TIMESTRUCT *Time    - Time information
    INPUTFILES *InFiles - Filenames for various input files
    MAPSIZE *Radar      - Information about radar coverage

  Returns      : void

  Modifies     : Members of Time, InFiles and Radar

  Comments     :
*****************************************************************************/
void InitRadar(LISTPTR Input, MAPSIZE * Map, TIMESTRUCT * Time,
  INPUTFILES * InFiles, MAPSIZE * Radar)
{
  DATE Start;
  int i;
  STRINIENTRY StrEnv[] = {
    {"METEOROLOGY", "RADAR START", "", ""},
    {"METEOROLOGY", "RADAR FILE", "", ""},
    {"METEOROLOGY", "RADAR EXTREME NORTH", "", ""},
    {"METEOROLOGY", "RADAR EXTREME WEST", "", ""},
    {"METEOROLOGY", "RADAR NUMBER OF ROWS", "", ""},
    {"METEOROLOGY", "RADAR NUMBER OF COLUMNS", "", ""},
    {"METEOROLOGY", "RADAR GRID SPACING", "", ""},
    {NULL, NULL, "", NULL},
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);

  /* Assign the entries to the variables */
  if (!SScanDate(StrEnv[radar_start].VarStr, &Start))
    ReportError(StrEnv[radar_start].KeyName, 51);

  InitTime(Time, NULL, NULL, &Start, NULL, Time->Dt);

  if (IsEmptyStr(StrEnv[radar_file].VarStr))
    ReportError(StrEnv[radar_file].KeyName, 51);
  strcpy(InFiles->RadarFile, StrEnv[radar_file].VarStr);

  /**************** Determine areal extent ****************/
  strcpy(Radar->System, Map->System);

  if (!CopyDouble(&(Radar->Yorig), StrEnv[radar_north].VarStr, 1))
    ReportError(StrEnv[radar_north].KeyName, 51);

  if (!CopyDouble(&(Radar->Xorig), StrEnv[radar_west].VarStr, 1))
    ReportError(StrEnv[radar_west].KeyName, 51);

  if (!CopyInt(&(Radar->NY), StrEnv[radar_rows].VarStr, 1))
    ReportError(StrEnv[radar_rows].KeyName, 51);

  if (!CopyInt(&(Radar->NX), StrEnv[radar_cols].VarStr, 1))
    ReportError(StrEnv[radar_cols].KeyName, 51);

  if (!CopyFloat(&(Radar->DY), StrEnv[radar_grid].VarStr, 1))
    ReportError(StrEnv[radar_grid].KeyName, 51);
  Radar->DX = Radar->DY;

  Radar->DXY = sqrt(Radar->DX * Radar->DX + Radar->DY * Radar->DY);     
  Radar->X = 0;
  Radar->Y = 0;
  Radar->OffsetX = Round(((float)(Radar->Xorig - Map->Xorig)) /
    ((float)Map->DX));
  Radar->OffsetY = Round(((float)(Radar->Yorig - Map->Yorig)) /
    ((float)Map->DY));

  if (Radar->OffsetX > 0 || Radar->OffsetY < 0)
    ReportError("Input Options File", 31);
}

/*******************************************************************************
  Function name: InitWindModel()

  Purpose      : Read the wind model information from the options file.  This
                 information is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    INPUTFILES *InFiles - Filenames for various input files
    int NStats          - Number of stations
    METLOCATION *Stat   - Information about met stations

  Returns      : void

  Modifies     : Members of InFiles and Stat.  Also a global constant
                 (NWINDMAPS)

  Comments     :
*****************************************************************************/
void InitWindModel(LISTPTR Input, INPUTFILES * InFiles, int NStats,
  METLOCATION * Stat)
{
  int i;
  int WindStation;
  STRINIENTRY StrEnv[] = {
    {"METEOROLOGY", "NUMBER OF WIND MAPS", "", ""},
    {"METEOROLOGY", "WIND FILE BASENAME", "", ""},
    {"METEOROLOGY", "WIND MAP MET STATION", "", ""},
    {NULL, NULL, "", NULL},
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);

  /* Assign the entries to the variables */
  if (!CopyInt(&NWINDMAPS, StrEnv[number_of_maps].VarStr, 1))
    ReportError(StrEnv[number_of_maps].KeyName, 51);

  if (IsEmptyStr(StrEnv[wind_map_path].VarStr))
    ReportError(StrEnv[wind_map_path].KeyName, 51);
  strcpy(InFiles->WindMapPath, StrEnv[wind_map_path].VarStr);

  if (!CopyInt(&WindStation, StrEnv[wind_station].VarStr, 1))
    ReportError(StrEnv[wind_station].KeyName, 51);

  if (WindStation < 1 || WindStation > NStats)
    ReportError(StrEnv[wind_station].KeyName, 53);

  for (i = 0; i < NStats; i++)
    Stat[i].IsWindModelLocation = FALSE;
  Stat[WindStation - 1].IsWindModelLocation = TRUE;
}

/*******************************************************************************
  Function name: InitPrecipLapse()

  Purpose      : Read the file name for the precip lapse rate from the options
                 file.  This information is in the [METEOROLOGY] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    INPUTFILES *InFiles - Filenames for various input files

  Returns      : void

  Modifies     : Member of InFiles

  Comments     :
*****************************************************************************/
void InitPrecipLapse(LISTPTR Input, INPUTFILES * InFiles)
{
  int i;			/* counter */

  STRINIENTRY StrEnv[] = {
    {"METEOROLOGY", "PRECIPITATION LAPSE RATE MAP", "", ""},
    {NULL, NULL, "", NULL},
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);

  /* Assign the entries to the variables */
  if (IsEmptyStr(StrEnv[precip_lapse_rate_file].VarStr))
    ReportError(StrEnv[precip_lapse_rate_file].KeyName, 51);
  strcpy(InFiles->PrecipLapseFile, StrEnv[precip_lapse_rate_file].VarStr);
}
