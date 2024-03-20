/*
* SUMMARY:      InitTables.c - Initialize lookup tables
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Initialize lookup tables
* DESCRIP-END.
* FUNCTIONS:    InitTables()
*               InitSoilTable()
*               InitVegTable()
*               InitSnowTable()
* COMMENTS:
* $Id: InitTables.c,v3.1.2 2013/12/11 ning Exp $
*/


#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "Calendar.h"
#include "constants.h"
#include "fileio.h"
#include "getinit.h"
#define NPARAM 105 //nparam+runnumber = 104+1 =105

/*******************************************************************************/
/*				  InitTables()                                 */

/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options, 
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType,
  LAYER *Veg, double Anovapara[NPARAM])
{
  int y, x; //counter
  printf("Initializing tables\n");

  if ((Soil->NTypes = InitSoilTable(Options, SType, Input, Soil, 
    Options->Infiltration, Anovapara)) == 0)
    ReportError("Input Options File", 8);

  if ((Veg->NTypes = InitVegTable(VType, Input, Options, Veg, Anovapara)) == 0)
    ReportError("Input Options File", 8);

  InitSatVaporTable();
}

/********************************************************************************
Function Name: InitSoilTable()

Purpose      : Initialize the soil lookup table
Processes most of the following section in InFileName:
[SOILS]

Required     :
SOILTABLE **SType - Pointer to lookup table
LISTPTR Input     - Pointer to linked list with input info
LAYER *Soil       - Pointer to structure with soil layer information

Returns      : Number of soil layers

Modifies     : SoilTable and Soil

Comments     :
********************************************************************************/
int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE ** SType,
  LISTPTR Input, LAYER * Soil, int InfiltOption, double Anovapara[NPARAM])
{
  const char *Routine = "InitSoilTable";
  int i;			/* counter */
  int j;			/* counter */
  int NSoils;			/* Number of soil types */
  char KeyName[thermal_capacity + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "SOIL DESCRIPTION",
    "LATERAL CONDUCTIVITY",
    "EXPONENTIAL DECREASE",
    "DEPTH THRESHOLD",
    "MAXIMUM INFILTRATION",
    "CAPILLARY DRIVE",
    "SURFACE ALBEDO",
    "NUMBER OF SOIL LAYERS",
    "POROSITY",
    "PORE SIZE DISTRIBUTION",
    "BUBBLING PRESSURE",
    "FIELD CAPACITY",
    "WILTING POINT",
    "BULK DENSITY",
    "VERTICAL CONDUCTIVITY",
    "THERMAL CONDUCTIVITY",
    "THERMAL CAPACITY"
  };
  char SectionName[] = "SOILS";
  char VarStr[thermal_capacity + 1][BUFSIZE + 1];


  /* Get the number of different soil types */
  GetInitString(SectionName, "NUMBER OF SOIL TYPES", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(&NSoils, VarStr[0], 1))
    ReportError("NUMBER OF SOIL TYPES", 51);

  if (NSoils == 0)
    return NSoils;

  if (!(Soil->NLayers = (int *)calloc(NSoils, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*SType = (SOILTABLE *)calloc(NSoils, sizeof(SOILTABLE))))
    ReportError((char *)Routine, 1);

  /********** Read information and allocate memory for each soil type *********/

  Soil->MaxLayers = 0;

  for (i = 0; i < NSoils; i++) {
      
      //the first soil(CLAY).
      if (i == 0) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= thermal_capacity; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[soil_description]))
              ReportError(KeyName[soil_description], 51);

          strcpy((*SType)[i].Desc, VarStr[soil_description]);
          (*SType)[i].Index = i;
          printf("(*SType)[%d].Desc:%s\n", i, (*SType)[i].Desc);

          if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
              ReportError(KeyName[lateral_ks], 51);
          (*SType)[i].KsLat = Anovapara[0]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLat:%lf\n", i, (*SType)[i].KsLat);

          if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
              ReportError(KeyName[exponent], 51);
          (*SType)[i].KsLatExp = Anovapara[1]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLatExp:%lf\n", i, (*SType)[i].KsLatExp);

          if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
              ReportError(KeyName[depth_thresh], 51);

          if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
              ReportError(KeyName[max_infiltration], 51);
          (*SType)[i].MaxInfiltrationRate = Anovapara[2]; //replace parameter in configfile by anova
          printf("(*SType)[%d].MaxInfiltrationRate:%lf\n", i, (*SType)[i].MaxInfiltrationRate);

          if (InfiltOption == DYNAMIC) {
              if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
                  ReportError(KeyName[capillary_drive], 51);
          }
          else (*SType)[i].G_Infilt = NOT_APPLICABLE;

          if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
              ReportError(KeyName[soil_albedo], 51);
          (*SType)[i].Albedo = Anovapara[3]; //replace parameter in configfile by anova
          printf("(*SType)[%d].Albedo:%lf\n", i, (*SType)[i].Albedo);

          if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
              ReportError(KeyName[number_of_layers], 51);
          Soil->NLayers[i] = (*SType)[i].NLayers;

          if (Soil->NLayers[i] > Soil->MaxLayers)
              Soil->MaxLayers = Soil->NLayers[i];

          /* allocate memory for the soil layers */
          if (!((*SType)[i].Porosity = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].PoreDist = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Press = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].FCap = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].WP = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Dens = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ks = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhDry = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhSol = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ch = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          sprintf(VarStr[porosity], "%lf %lf %lf", Anovapara[4], 0.93 * Anovapara[4], 0.9 * Anovapara[4]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity],
              (*SType)[i].NLayers))
              ReportError(KeyName[porosity], 51);
              printf("(*SType)[%d].Porosity:%lf\n", i, *(*SType)[i].Porosity);

          sprintf(VarStr[pore_size], "%lf %lf %lf", Anovapara[5], Anovapara[5], Anovapara[5]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
              (*SType)[i].NLayers))
              ReportError(KeyName[pore_size], 51);
          printf("(*SType)[%d].PoreDist:%lf\n", i, *(*SType)[i].PoreDist);

          sprintf(VarStr[bubbling_pressure], "%lf %lf %lf", Anovapara[6], Anovapara[6], Anovapara[6]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
              (*SType)[i].NLayers))
              ReportError(KeyName[bubbling_pressure], 51);
          printf("(*SType)[%d].Press:%lf\n", i, *(*SType)[i].Press);

          sprintf(VarStr[field_capacity], "%lf %lf %lf", Anovapara[7], 0.93 * Anovapara[7], 0.9 * Anovapara[7]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
              (*SType)[i].NLayers))
              ReportError(KeyName[field_capacity], 51);
          printf("(*SType)[%d].FCap:%lf\n", i, *(*SType)[i].FCap);

          sprintf(VarStr[wilting_point], "%lf %lf %lf", Anovapara[8], 0.93 * Anovapara[8], 0.9 * Anovapara[8]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
              ReportError(KeyName[wilting_point], 51);
          printf("(*SType)[%d].WP:%lf\n", i, *(*SType)[i].WP);

          sprintf(VarStr[bulk_density], "%lf %lf %lf", Anovapara[9], Anovapara[9], Anovapara[9]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
              ReportError(KeyName[bulk_density], 51);
          printf("(*SType)[%d].Dens:%lf\n", i, *(*SType)[i].Dens);

          sprintf(VarStr[vertical_ks], "%lf %lf %lf", Anovapara[10], Anovapara[10], Anovapara[10]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
              ReportError(KeyName[vertical_ks], 51);
          printf("(*SType)[%d].Ks:%lf\n", i, *(*SType)[i].Ks);

          sprintf(VarStr[solids_thermal], "%lf %lf %lf", Anovapara[11], Anovapara[11], Anovapara[11]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
              ReportError(KeyName[solids_thermal], 51);
          printf("(*SType)[%d].KhSol:%lf\n", i, *(*SType)[i].KhSol);

          sprintf(VarStr[thermal_capacity], "%lf %lf %lf", Anovapara[12], Anovapara[12], Anovapara[12]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
              ReportError(KeyName[thermal_capacity], 51);
          printf("(*SType)[%d].Ch:%lf\n", i, *(*SType)[i].Ch);
      }

      //the second soil(SILTY LOAM)
      if (i == 1) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= thermal_capacity; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[soil_description]))
              ReportError(KeyName[soil_description], 51);

          strcpy((*SType)[i].Desc, VarStr[soil_description]);
          (*SType)[i].Index = i;
          printf("(*SType)[%d].Desc:%s\n", i, (*SType)[i].Desc);

          if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
              ReportError(KeyName[lateral_ks], 51);
          (*SType)[i].KsLat = Anovapara[13]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLat:%lf\n", i, (*SType)[i].KsLat);

          if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
              ReportError(KeyName[exponent], 51);
          (*SType)[i].KsLatExp = Anovapara[14]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLatExp:%lf\n", i, (*SType)[i].KsLatExp);

          if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
              ReportError(KeyName[depth_thresh], 51);

          if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
              ReportError(KeyName[max_infiltration], 51);
          (*SType)[i].MaxInfiltrationRate = Anovapara[15]; //replace parameter in configfile by anova
          printf("(*SType)[%d].MaxInfiltrationRate:%lf\n", i, (*SType)[i].MaxInfiltrationRate);

          if (InfiltOption == DYNAMIC) {
              if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
                  ReportError(KeyName[capillary_drive], 51);
          }
          else (*SType)[i].G_Infilt = NOT_APPLICABLE;

          if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
              ReportError(KeyName[soil_albedo], 51);
          (*SType)[i].Albedo = Anovapara[16]; //replace parameter in configfile by anova
          printf("(*SType)[%d].Albedo:%lf\n", i, (*SType)[i].Albedo);

          if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
              ReportError(KeyName[number_of_layers], 51);
          Soil->NLayers[i] = (*SType)[i].NLayers;

          if (Soil->NLayers[i] > Soil->MaxLayers)
              Soil->MaxLayers = Soil->NLayers[i];

          /* allocate memory for the soil layers */
          if (!((*SType)[i].Porosity = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].PoreDist = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Press = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].FCap = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].WP = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Dens = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ks = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhDry = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhSol = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ch = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          sprintf(VarStr[porosity], "%lf %lf %lf", Anovapara[17], 0.93 * Anovapara[17], 0.9 * Anovapara[17]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
              ReportError(KeyName[porosity], 51);
          printf("(*SType)[%d].Porosity:%lf\n", i, *(*SType)[i].Porosity);

          sprintf(VarStr[pore_size], "%lf %lf %lf", Anovapara[18], Anovapara[18], Anovapara[18]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
              (*SType)[i].NLayers))
              ReportError(KeyName[pore_size], 51);
          printf("(*SType)[%d].PoreDist:%lf\n", i, *(*SType)[i].PoreDist);

          sprintf(VarStr[bubbling_pressure], "%lf %lf %lf", Anovapara[19], Anovapara[19], Anovapara[19]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
              (*SType)[i].NLayers))
              ReportError(KeyName[bubbling_pressure], 51);
          printf("(*SType)[%d].Press:%lf\n", i, *(*SType)[i].Press);

          sprintf(VarStr[field_capacity], "%lf %lf %lf", Anovapara[20], 0.93 * Anovapara[20], 0.9 * Anovapara[20]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
              (*SType)[i].NLayers))
              ReportError(KeyName[field_capacity], 51);
          printf("(*SType)[%d].FCap:%lf\n", i, *(*SType)[i].FCap);

          sprintf(VarStr[wilting_point], "%lf %lf %lf", Anovapara[21], 0.93 * Anovapara[21], 0.9 * Anovapara[21]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
              ReportError(KeyName[wilting_point], 51);
          printf("(*SType)[%d].WP:%lf\n", i, *(*SType)[i].WP);

          sprintf(VarStr[bulk_density], "%lf %lf %lf", Anovapara[22], Anovapara[22], Anovapara[22]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
              ReportError(KeyName[bulk_density], 51);
          printf("(*SType)[%d].Dens:%lf\n", i, *(*SType)[i].Dens);

          sprintf(VarStr[vertical_ks], "%lf %lf %lf", Anovapara[23], Anovapara[23], Anovapara[23]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
              ReportError(KeyName[vertical_ks], 51);
          printf("(*SType)[%d].Ks:%lf\n", i, *(*SType)[i].Ks);

          sprintf(VarStr[solids_thermal], "%lf %lf %lf", Anovapara[24], Anovapara[24], Anovapara[24]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
              ReportError(KeyName[solids_thermal], 51);
          printf("(*SType)[%d].KhSol:%lf\n", i, *(*SType)[i].KhSol);

          sprintf(VarStr[thermal_capacity], "%lf %lf %lf", Anovapara[25], Anovapara[25], Anovapara[25]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
              ReportError(KeyName[thermal_capacity], 51);
          printf("(*SType)[%d].Ch:%lf\n", i, *(*SType)[i].Ch);
      }

      //the third soil(LOAM).
      if (i == 2) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= thermal_capacity; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[soil_description]))
              ReportError(KeyName[soil_description], 51);

          strcpy((*SType)[i].Desc, VarStr[soil_description]);
          (*SType)[i].Index = i;
          printf("(*SType)[%d].Desc:%s\n", i, (*SType)[i].Desc);

          if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
              ReportError(KeyName[lateral_ks], 51);
          (*SType)[i].KsLat = Anovapara[26]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLat:%lf\n", i, (*SType)[i].KsLat);

          if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
              ReportError(KeyName[exponent], 51);
          (*SType)[i].KsLatExp = Anovapara[27]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLatExp:%lf\n", i, (*SType)[i].KsLatExp);

          if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
              ReportError(KeyName[depth_thresh], 51);

          if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
              ReportError(KeyName[max_infiltration], 51);
          (*SType)[i].MaxInfiltrationRate = Anovapara[28]; //replace parameter in configfile by anova
          printf("(*SType)[%d].MaxInfiltrationRate:%lf\n", i, (*SType)[i].MaxInfiltrationRate);

          if (InfiltOption == DYNAMIC) {
              if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
                  ReportError(KeyName[capillary_drive], 51);
          }
          else (*SType)[i].G_Infilt = NOT_APPLICABLE;

          if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
              ReportError(KeyName[soil_albedo], 51);
          (*SType)[i].Albedo = Anovapara[29]; //replace parameter in configfile by anova
          printf("(*SType)[%d].Albedo:%lf\n", i, (*SType)[i].Albedo);

          if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
              ReportError(KeyName[number_of_layers], 51);
          Soil->NLayers[i] = (*SType)[i].NLayers;

          if (Soil->NLayers[i] > Soil->MaxLayers)
              Soil->MaxLayers = Soil->NLayers[i];

          /* allocate memory for the soil layers */
          if (!((*SType)[i].Porosity = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].PoreDist = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Press = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].FCap = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].WP = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Dens = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ks = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhDry = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhSol = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ch = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          sprintf(VarStr[porosity], "%lf %lf %lf", Anovapara[30], 0.93 * Anovapara[30], 0.9 * Anovapara[30]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
              ReportError(KeyName[porosity], 51);
          printf("(*SType)[%d].Porosity:%lf\n", i, *(*SType)[i].Porosity);

          sprintf(VarStr[pore_size], "%lf %lf %lf", Anovapara[31], Anovapara[31], Anovapara[31]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
              (*SType)[i].NLayers))
              ReportError(KeyName[pore_size], 51);
          printf("(*SType)[%d].PoreDist:%lf\n", i, *(*SType)[i].PoreDist);

          sprintf(VarStr[bubbling_pressure], "%lf %lf %lf", Anovapara[32], Anovapara[32], Anovapara[32]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
              (*SType)[i].NLayers))
              ReportError(KeyName[bubbling_pressure], 51);
          printf("(*SType)[%d].Press:%lf\n", i, *(*SType)[i].Press);

          sprintf(VarStr[field_capacity], "%lf %lf %lf", Anovapara[33], 0.93 * Anovapara[33], 0.9 * Anovapara[33]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
              (*SType)[i].NLayers))
              ReportError(KeyName[field_capacity], 51);
          printf("(*SType)[%d].FCap:%lf\n", i, *(*SType)[i].FCap);

          sprintf(VarStr[wilting_point], "%lf %lf %lf", Anovapara[34], 0.93 * Anovapara[34], 0.9 * Anovapara[34]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
              ReportError(KeyName[wilting_point], 51);
          printf("(*SType)[%d].WP:%lf\n", i, *(*SType)[i].WP);

          sprintf(VarStr[bulk_density], "%lf %lf %lf", Anovapara[35], Anovapara[35], Anovapara[35]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
              ReportError(KeyName[bulk_density], 51);
          printf("(*SType)[%d].Dens:%lf\n", i, *(*SType)[i].Dens);

          sprintf(VarStr[vertical_ks], "%lf %lf %lf", Anovapara[36], Anovapara[36], Anovapara[36]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
              ReportError(KeyName[vertical_ks], 51);
          printf("(*SType)[%d].Ks:%lf\n", i, *(*SType)[i].Ks);

          sprintf(VarStr[solids_thermal], "%lf %lf %lf", Anovapara[37], Anovapara[37], Anovapara[37]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
              ReportError(KeyName[solids_thermal], 51);
          printf("(*SType)[%d].KhSol:%lf\n", i, *(*SType)[i].KhSol);

          sprintf(VarStr[thermal_capacity], "%lf %lf %lf", Anovapara[38], Anovapara[38], Anovapara[38]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
              ReportError(KeyName[thermal_capacity], 51);
          printf("(*SType)[%d].Ch:%lf\n", i, *(*SType)[i].Ch);
      }
      
      //the fourth soil(SANDY CLAY LOAM)
      if (i == 3) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= thermal_capacity; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[soil_description]))
              ReportError(KeyName[soil_description], 51);

          strcpy((*SType)[i].Desc, VarStr[soil_description]);
          (*SType)[i].Index = i;
          printf("(*SType)[%d].Desc:%s\n", i, (*SType)[i].Desc);

          if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
              ReportError(KeyName[lateral_ks], 51);
          (*SType)[i].KsLat = Anovapara[39]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLat:%lf\n", i, (*SType)[i].KsLat);

          if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
              ReportError(KeyName[exponent], 51);
          (*SType)[i].KsLatExp = Anovapara[40]; //replace parameter in configfile by anova
          printf("(*SType)[%d].KsLatExp:%lf\n", i, (*SType)[i].KsLatExp);

          if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
              ReportError(KeyName[depth_thresh], 51);

          if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
              ReportError(KeyName[max_infiltration], 51);
          (*SType)[i].MaxInfiltrationRate = Anovapara[41]; //replace parameter in configfile by anova
          printf("(*SType)[%d].MaxInfiltrationRate:%lf\n", i, (*SType)[i].MaxInfiltrationRate);

          if (InfiltOption == DYNAMIC) {
              if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
                  ReportError(KeyName[capillary_drive], 51);
          }
          else (*SType)[i].G_Infilt = NOT_APPLICABLE;

          if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
              ReportError(KeyName[soil_albedo], 51);
          (*SType)[i].Albedo = Anovapara[42]; //replace parameter in configfile by anova
          printf("(*SType)[%d].Albedo:%lf\n", i, (*SType)[i].Albedo);

          if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
              ReportError(KeyName[number_of_layers], 51);
          Soil->NLayers[i] = (*SType)[i].NLayers;

          if (Soil->NLayers[i] > Soil->MaxLayers)
              Soil->MaxLayers = Soil->NLayers[i];

          /* allocate memory for the soil layers */
          if (!((*SType)[i].Porosity = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].PoreDist = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Press = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].FCap = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].WP = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Dens = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ks = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhDry = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].KhSol = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
          if (!((*SType)[i].Ch = (float*)calloc((*SType)[i].NLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          sprintf(VarStr[porosity], "%lf %lf %lf", Anovapara[43], 0.93 * Anovapara[43], 0.9 * Anovapara[43]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
              ReportError(KeyName[porosity], 51);
          printf("(*SType)[%d].Porosity:%lf\n", i, *(*SType)[i].Porosity);

          sprintf(VarStr[pore_size], "%lf %lf %lf", Anovapara[44], Anovapara[44], Anovapara[44]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
              (*SType)[i].NLayers))
              ReportError(KeyName[pore_size], 51);
          printf("(*SType)[%d].PoreDist:%lf\n", i, *(*SType)[i].PoreDist);

          sprintf(VarStr[bubbling_pressure], "%lf %lf %lf", Anovapara[45], Anovapara[45], Anovapara[45]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
              (*SType)[i].NLayers))
              ReportError(KeyName[bubbling_pressure], 51);
          printf("(*SType)[%d].Press:%lf\n", i, *(*SType)[i].Press);

          sprintf(VarStr[field_capacity], "%lf %lf %lf", Anovapara[46], 0.93 * Anovapara[46], 0.9 * Anovapara[46]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
              (*SType)[i].NLayers))
              ReportError(KeyName[field_capacity], 51);
          printf("(*SType)[%d].FCap:%lf\n", i, *(*SType)[i].FCap);

          sprintf(VarStr[wilting_point], "%lf %lf %lf", Anovapara[47], 0.93 * Anovapara[47], 0.9 * Anovapara[47]);  //read para will be calibrated.
          if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
              ReportError(KeyName[wilting_point], 51);
          printf("(*SType)[%d].WP:%lf\n", i, *(*SType)[i].WP);

          sprintf(VarStr[bulk_density], "%lf %lf %lf", Anovapara[48], Anovapara[48], Anovapara[48]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
              ReportError(KeyName[bulk_density], 51);
          printf("(*SType)[%d].Dens:%lf\n", i, *(*SType)[i].Dens);

          sprintf(VarStr[vertical_ks], "%lf %lf %lf", Anovapara[49], Anovapara[49], Anovapara[49]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
              ReportError(KeyName[vertical_ks], 51);
          printf("(*SType)[%d].Ks:%lf\n", i, *(*SType)[i].Ks);

          sprintf(VarStr[solids_thermal], "%lf %lf %lf", Anovapara[50], Anovapara[50], Anovapara[50]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
              ReportError(KeyName[solids_thermal], 51);
          printf("(*SType)[%d].KhSol:%lf\n", i, *(*SType)[i].KhSol);

          sprintf(VarStr[thermal_capacity], "%lf %lf %lf", Anovapara[51], Anovapara[51], Anovapara[51]); //read para will be sensitivity analysised.
          if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
              ReportError(KeyName[thermal_capacity], 51);
          printf("(*SType)[%d].Ch:%lf\n", i, *(*SType)[i].Ch);
      }

       else {
       /* Read the key-entry pairs from the input file */
       for (j = 0; j <= thermal_capacity; j++) {
           sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
           GetInitString(SectionName, KeyName[j], "", VarStr[j],
               (unsigned long)BUFSIZE, Input);
       }

       /* Assign the entries to the appropriate variables */
       if (IsEmptyStr(VarStr[soil_description]))
           ReportError(KeyName[soil_description], 51);

       strcpy((*SType)[i].Desc, VarStr[soil_description]);
       (*SType)[i].Index = i;

       if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
           ReportError(KeyName[lateral_ks], 51);

       if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
           ReportError(KeyName[exponent], 51);

       if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
           ReportError(KeyName[depth_thresh], 51);

       if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
           ReportError(KeyName[max_infiltration], 51);

       if (InfiltOption == DYNAMIC) {
           if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
               ReportError(KeyName[capillary_drive], 51);
       }
       else (*SType)[i].G_Infilt = NOT_APPLICABLE;

       if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
           ReportError(KeyName[soil_albedo], 51);

       if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
           ReportError(KeyName[number_of_layers], 51);
       Soil->NLayers[i] = (*SType)[i].NLayers;

       if (Soil->NLayers[i] > Soil->MaxLayers)
           Soil->MaxLayers = Soil->NLayers[i];

       /* allocate memory for the soil layers */
       if (!((*SType)[i].Porosity = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].PoreDist = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].Press = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].FCap = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].WP = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].Dens = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].Ks = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].KhDry = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].KhSol = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);
       if (!((*SType)[i].Ch = (float*)calloc((*SType)[i].NLayers,
           sizeof(float))))
           ReportError((char*)Routine, 1);

       if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
           ReportError(KeyName[porosity], 51);

       if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
           (*SType)[i].NLayers))
           ReportError(KeyName[pore_size], 51);

       if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
           (*SType)[i].NLayers))
           ReportError(KeyName[bubbling_pressure], 51);

       if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
           (*SType)[i].NLayers))
           ReportError(KeyName[field_capacity], 51);

       if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
           ReportError(KeyName[wilting_point], 51);

       if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
           ReportError(KeyName[bulk_density], 51);

       if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
           ReportError(KeyName[vertical_ks], 51);

       if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal], (*SType)[i].NLayers))
           ReportError(KeyName[solids_thermal], 51);

       if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity], (*SType)[i].NLayers))
           ReportError(KeyName[thermal_capacity], 51);
    }
  }

  for (i = 0; i < NSoils; i++)
    for (j = 0; j < (*SType)[i].NLayers; j++) {
      (*SType)[i].KhDry[j] = CalcKhDry((*SType)[i].Dens[j]);
      if (((*SType)[i].Porosity[j] < (*SType)[i].FCap[j])
        || ((*SType)[i].Porosity[j] < (*SType)[i].WP[j])
        || ((*SType)[i].FCap[j] < (*SType)[i].WP[j]))
        ReportError((*SType)[i].Desc, 11);
    }

  return NSoils;
}

/********************************************************************************
Function Name: InitVegTable()

Purpose      : Initialize the vegetation lookup table
Processes most of the following section in the input file:
[VEGETATION]

Required     :
VEGTABLE **VType - Pointer to lookup table
LISTPTR Input    - Pointer to linked list with input info
LAYER *Veg       - Pointer to structure with veg layer information

Returns      : Number of vegetation types

Modifies     : VegTable and Veg

Comments     :
********************************************************************************/
int InitVegTable(VEGTABLE **VType, LISTPTR Input, OPTIONSTRUCT *Options, LAYER *Veg, double Anovapara[NPARAM])
{
  const char *Routine = "InitVegTable";
  int i;			/* Counter */
  int j;			/* Counter */
  int k;            /* counter */
  float impervious;	/* flag to check whether impervious layers are specified */

  int NVegs;		/* Number of vegetation types */

  char KeyName[understory_monalb + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "VEGETATION DESCRIPTION",
    "OVERSTORY PRESENT",
    "UNDERSTORY PRESENT",
    "FRACTIONAL COVERAGE",
    "HEMI FRACT COVERAGE",
    "TRUNK SPACE",
    "AERODYNAMIC ATTENUATION",
    "RADIATION ATTENUATION",
    "DIFFUSE RADIATION ATTENUATION",
    "CLUMPING FACTOR",
    "LEAF ANGLE A",
    "LEAF ANGLE B",
    "SCATTERING PARAMETER",
    "MAX SNOW INT CAPACITY",
    "MASS RELEASE DRIP RATIO",
    "SNOW INTERCEPTION EFF",
    "IMPERVIOUS FRACTION",
    "DETENTION FRACTION",
    "DETENTION DECAY",
    "HEIGHT",
    "MAXIMUM RESISTANCE",
    "MINIMUM RESISTANCE",
    "MOISTURE THRESHOLD",
    "VAPOR PRESSURE DEFICIT",
    "RPC",
    "NUMBER OF ROOT ZONES",
    "ROOT ZONE DEPTHS",
    "OVERSTORY ROOT FRACTION",
    "UNDERSTORY ROOT FRACTION",
    "MONTHLY LIGHT EXTINCTION",
    "CANOPY VIEW ADJ FACTOR",
    "OVERSTORY MONTHLY LAI",
    "UNDERSTORY MONTHLY LAI",
    "OVERSTORY MONTHLY ALB",
    "UNDERSTORY MONTHLY ALB"
  };
  char SectionName[] = "VEGETATION";
  char VarStr[understory_monalb + 1][BUFSIZE + 1];
  float maxLAI;

  /* Get the number of different vegetation types */
  GetInitString(SectionName, "NUMBER OF VEGETATION TYPES", "", VarStr[0],
    (unsigned long)BUFSIZE, Input);
  if (!CopyInt(&NVegs, VarStr[0], 1))
    ReportError("NUMBER OF VEGETATION TYPES", 51);

  if (NVegs == 0)
    return NVegs;

  if (!(Veg->NLayers = (int *)calloc(NVegs, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*VType = (VEGTABLE *)calloc(NVegs, sizeof(VEGTABLE))))
    ReportError((char *)Routine, 1);

   /******* Read information and allocate memory for each vegetation type ******/

  Veg->MaxLayers = 0;
  impervious = 0.0;
  for (i = 0; i < NVegs; i++) {
      
      //the first vegetation(RAINFED CROPLAND)
      if (i == 0) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= understory_monalb; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[veg_description]))
              ReportError(KeyName[veg_description], 51);
          strcpy((*VType)[i].Desc, VarStr[veg_description]);
          MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                                  uppercase and removed spaces so
                                                  it is easier to compare */
          if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
              (*VType)[i].Index = GLACIER;
          }
          else
              (*VType)[i].Index = i;

          (*VType)[i].NVegLayers = 0;

          if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
              (*VType)[i].OverStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
              (*VType)[i].OverStory = FALSE;
          else
              ReportError(KeyName[overstory], 51);

          if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
              (*VType)[i].UnderStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
              (*VType)[i].UnderStory = FALSE;
          else
              ReportError(KeyName[understory], 51);

          Veg->NLayers[i] = (*VType)[i].NVegLayers;
          if ((*VType)[i].NVegLayers > Veg->MaxLayers)
              Veg->MaxLayers = (*VType)[i].NVegLayers;

          if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
              ReportError(KeyName[number_of_root_zones], 51);

          if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
              ReportError(KeyName[imperv_frac], 51);
          impervious += (*VType)[i].ImpervFrac;

          if ((*VType)[i].ImpervFrac > 0) {
              if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
                  ReportError(KeyName[detention_frac], 51);
              if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
                  ReportError(KeyName[detention_decay], 51);
          }
          else {
              (*VType)[i].DetentionFrac = 0.;
              (*VType)[i].DetentionDecay = 0.;
          }

          /* allocate memory for the vegetation layers */
          if (!((*VType)[i].Fract = (float*)calloc((*VType)[i].NVegLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          if (Options->CanopyRadAtt == VARIABLE) {
              if (!((*VType)[i].HemiFract = (float*)calloc((*VType)[i].NVegLayers,
                  sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          else {
              (*VType)[i].HemiFract = NULL;
          }

          if (!((*VType)[i].Height = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMax = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMin = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MoistThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].VpdThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Rpc = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Albedo = (float*)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MaxInt = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))

              ReportError((char*)Routine, 1);
          if (!((*VType)[i].LAI = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RootFract = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);

          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].RootFract[j] =
                  (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          if (!((*VType)[i].RootDepth = (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].LAIMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].LAIMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          if (!((*VType)[i].AlbedoMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].AlbedoMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          /* assign the entries to the appropriate variables */
          /* allocation of zero memory is not supported on some
          compilers */
          if ((*VType)[i].OverStory == TRUE) {
              if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
                  ReportError(KeyName[fraction], 51);
              if (Options->CanopyRadAtt == VARIABLE) {
                  if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
                      ReportError(KeyName[hemifraction], 51);
                  if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
                      ReportError(KeyName[clumping_factor], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
                      ReportError(KeyName[leaf_angle_a], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
                      ReportError(KeyName[leaf_angle_b], 51);
                  if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
                      ReportError(KeyName[scat], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
              }
              else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
                  if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
                      ReportError(KeyName[beam_attn], 51);
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
                  (*VType)[i].Taud = NOT_APPLICABLE;
              }
              else if (Options->ImprovRadiation == TRUE) {
                  if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
                      ReportError(KeyName[diff_attn], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
              }

              if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
                  ReportError(KeyName[trunk_space], 51);

              if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
                  ReportError(KeyName[aerodynamic_att], 51);

              if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
                  ReportError(KeyName[snow_int_cap], 51);

              if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
                  ReportError(KeyName[mass_drip_ratio], 51);

              if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
                  ReportError(KeyName[snow_int_eff], 51);

              if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[overstory_fraction], 51);

              if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
                  ReportError(KeyName[overstory_monlai], 51);

              maxLAI = -9999;
              for (k = 0; k < 12; k++) {
                  if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
                      maxLAI = (*VType)[i].LAIMonthly[0][k];
              }

              if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
                  ReportError(KeyName[overstory_monalb], 51);

              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[1] = 1.0;
                  if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);

                  if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);

                  if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
              }
          }
          else { // e.g.OverStory = FALSE, UnderStory = TRUE.
              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[0] = 1.0;
                  sprintf(VarStr[understory_fraction], "%lf %lf %lf", Anovapara[52], 1 - Anovapara[52], 0.00); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);
                  printf("(*VType)[%d].RootFract[0]:%lf\n", i, *(*VType)[i].RootFract[0]);

                  sprintf(VarStr[understory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 1.02 * Anovapara[53], 1.01 * Anovapara[53],
                      1.02 * Anovapara[53], 1.12 * Anovapara[53], 1.34 * Anovapara[53], 1.79 * Anovapara[53], 1.72 * Anovapara[53],
                      1.36 * Anovapara[53], 1.76 * Anovapara[53], 1.05 * Anovapara[53], Anovapara[53], Anovapara[53]);
                  //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);
                  printf("(*VType)[%d].LAIMonthly[0]:%lf\n", i, *(*VType)[i].LAIMonthly[0]);

                  sprintf(VarStr[understory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[54], Anovapara[54], Anovapara[54],
                      Anovapara[54], Anovapara[54], Anovapara[54], Anovapara[54], Anovapara[54], Anovapara[54], Anovapara[54], Anovapara[54],
                      Anovapara[54]);
                  if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
                         printf("(*VType)[%d].AlbedoMonthly[0]:%lf\n", i, *(*VType)[i].AlbedoMonthly[0]);
              }
              (*VType)[i].Trunk = NOT_APPLICABLE;
              (*VType)[i].Cn = NOT_APPLICABLE;
              (*VType)[i].Atten = NOT_APPLICABLE;
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
          }

          sprintf(VarStr[height], "%lf", Anovapara[55]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
              ReportError(KeyName[height], 51);
          printf("(*VType)[%d].Height:%lf\n", i, *(*VType)[i].Height);

          sprintf(VarStr[max_resistance], "%lf", Anovapara[56]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[max_resistance], 51);
          printf("(*VType)[%d].RsMax:%lf\n", i, *(*VType)[i].RsMax);

          sprintf(VarStr[min_resistance], "%lf", Anovapara[57]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[min_resistance], 51);
          printf("(*VType)[%d].RsMin:%lf\n", i, *(*VType)[i].RsMin);

          sprintf(VarStr[moisture_threshold], "%lf", Anovapara[58]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[moisture_threshold], 51);
          printf("(*VType)[%d].MoistThres:%lf\n", i, *(*VType)[i].MoistThres);

          if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
              ReportError(KeyName[vpd], 51);

          if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
              ReportError(KeyName[rpc], 51);

          sprintf(VarStr[root_zone_depth], "%lf %lf %lf", Anovapara[59], Anovapara[59], Anovapara[59]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
              (*VType)[i].NSoilLayers))
              ReportError(KeyName[root_zone_depth], 51);
          printf("(*VType)[%d].RootDepth:%lf\n", i, *(*VType)[i].RootDepth);

          /* Calculate the wind speed profiles and the aerodynamical resistances
          for each layer.  The values are normalized for a reference height wind
          speed of 1 m/s, and are adjusted each timestep using actual reference
          height wind speeds */
          CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
              (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
              (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
              &((*VType)[i].RaSnow));
      }

      //the second vegetation(IRRIGATED CROPLAND)
      if (i == 2) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= understory_monalb; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[veg_description]))
              ReportError(KeyName[veg_description], 51);
          strcpy((*VType)[i].Desc, VarStr[veg_description]);
          MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                                  uppercase and removed spaces so
                                                  it is easier to compare */
          if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
              (*VType)[i].Index = GLACIER;
          }
          else
              (*VType)[i].Index = i;

          (*VType)[i].NVegLayers = 0;

          if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
              (*VType)[i].OverStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
              (*VType)[i].OverStory = FALSE;
          else
              ReportError(KeyName[overstory], 51);

          if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
              (*VType)[i].UnderStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
              (*VType)[i].UnderStory = FALSE;
          else
              ReportError(KeyName[understory], 51);

          Veg->NLayers[i] = (*VType)[i].NVegLayers;
          if ((*VType)[i].NVegLayers > Veg->MaxLayers)
              Veg->MaxLayers = (*VType)[i].NVegLayers;

          if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
              ReportError(KeyName[number_of_root_zones], 51);

          if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
              ReportError(KeyName[imperv_frac], 51);
          impervious += (*VType)[i].ImpervFrac;

          if ((*VType)[i].ImpervFrac > 0) {
              if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
                  ReportError(KeyName[detention_frac], 51);
              if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
                  ReportError(KeyName[detention_decay], 51);
          }
          else {
              (*VType)[i].DetentionFrac = 0.;
              (*VType)[i].DetentionDecay = 0.;
          }

          /* allocate memory for the vegetation layers */
          if (!((*VType)[i].Fract = (float*)calloc((*VType)[i].NVegLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          if (Options->CanopyRadAtt == VARIABLE) {
              if (!((*VType)[i].HemiFract = (float*)calloc((*VType)[i].NVegLayers,
                  sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          else {
              (*VType)[i].HemiFract = NULL;
          }

          if (!((*VType)[i].Height = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMax = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMin = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MoistThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].VpdThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Rpc = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Albedo = (float*)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MaxInt = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))

              ReportError((char*)Routine, 1);
          if (!((*VType)[i].LAI = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RootFract = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);

          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].RootFract[j] =
                  (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          if (!((*VType)[i].RootDepth = (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].LAIMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].LAIMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          if (!((*VType)[i].AlbedoMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].AlbedoMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          /* assign the entries to the appropriate variables */
          /* allocation of zero memory is not supported on some
          compilers */
          if ((*VType)[i].OverStory == TRUE) {
              if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
                  ReportError(KeyName[fraction], 51);
              if (Options->CanopyRadAtt == VARIABLE) {
                  if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
                      ReportError(KeyName[hemifraction], 51);
                  if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
                      ReportError(KeyName[clumping_factor], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
                      ReportError(KeyName[leaf_angle_a], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
                      ReportError(KeyName[leaf_angle_b], 51);
                  if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
                      ReportError(KeyName[scat], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
              }
              else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
                  if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
                      ReportError(KeyName[beam_attn], 51);
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
                  (*VType)[i].Taud = NOT_APPLICABLE;
              }
              else if (Options->ImprovRadiation == TRUE) {
                  if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
                      ReportError(KeyName[diff_attn], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
              }

              if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
                  ReportError(KeyName[trunk_space], 51);

              if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
                  ReportError(KeyName[aerodynamic_att], 51);

              if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
                  ReportError(KeyName[snow_int_cap], 51);

              if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
                  ReportError(KeyName[mass_drip_ratio], 51);

              if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
                  ReportError(KeyName[snow_int_eff], 51);

              if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[overstory_fraction], 51);

              if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
                  ReportError(KeyName[overstory_monlai], 51);

              maxLAI = -9999;
              for (k = 0; k < 12; k++) {
                  if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
                      maxLAI = (*VType)[i].LAIMonthly[0][k];
              }

              if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
                  ReportError(KeyName[overstory_monalb], 51);

              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[1] = 1.0;
                  if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);

                  if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);

                  if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
              }
          }
          else { // e.g.OverStory = FALSE, UnderStory = TRUE.
              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[0] = 1.0;
                  sprintf(VarStr[understory_fraction], "%lf %lf %lf", Anovapara[60], 1 - Anovapara[60], 0.00); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);
                  printf("(*VType)[%d].RootFract[0]:%lf\n", i, *(*VType)[i].RootFract[0]);

                  sprintf(VarStr[understory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 1.02 * Anovapara[61], 1.01 * Anovapara[61],
                      1.02 * Anovapara[61], 1.12 * Anovapara[61], 1.34 * Anovapara[61], 1.79 * Anovapara[61], 1.72 * Anovapara[61],
                      1.36 * Anovapara[61], 1.76 * Anovapara[61], 1.05 * Anovapara[61], Anovapara[61], Anovapara[61]);
                  //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);
                  printf("(*VType)[%d].LAIMonthly[0]:%lf\n", i, *(*VType)[i].LAIMonthly[0]);

                  sprintf(VarStr[understory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[62], Anovapara[62], Anovapara[62],
                      Anovapara[62], Anovapara[62], Anovapara[62], Anovapara[62], Anovapara[62], Anovapara[62], Anovapara[62], Anovapara[62],
                      Anovapara[62]);
                  if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
                  printf("(*VType)[%d].AlbedoMonthly[0]:%lf\n", i, *(*VType)[i].AlbedoMonthly[0]);
              }
              (*VType)[i].Trunk = NOT_APPLICABLE;
              (*VType)[i].Cn = NOT_APPLICABLE;
              (*VType)[i].Atten = NOT_APPLICABLE;
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
          }

          sprintf(VarStr[height], "%lf", Anovapara[63]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
              ReportError(KeyName[height], 51);
          printf("(*VType)[%d].Height:%lf\n", i, *(*VType)[i].Height);

          sprintf(VarStr[max_resistance], "%lf", Anovapara[64]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[max_resistance], 51);
          printf("(*VType)[%d].RsMax:%lf\n", i, *(*VType)[i].RsMax);

          sprintf(VarStr[min_resistance], "%lf", Anovapara[65]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[min_resistance], 51);
          printf("(*VType)[%d].RsMin:%lf\n", i, *(*VType)[i].RsMin);

          sprintf(VarStr[moisture_threshold], "%lf", Anovapara[66]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[moisture_threshold], 51);
          printf("(*VType)[%d].MoistThres:%lf\n", i, *(*VType)[i].MoistThres);

          if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
              ReportError(KeyName[vpd], 51);

          if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
              ReportError(KeyName[rpc], 51);

          sprintf(VarStr[root_zone_depth], "%lf %lf %lf", Anovapara[67], Anovapara[67], Anovapara[67]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
              (*VType)[i].NSoilLayers))
              ReportError(KeyName[root_zone_depth], 51);
          printf("(*VType)[%d].RootDepth:%lf\n", i, *(*VType)[i].RootDepth);

          /* Calculate the wind speed profiles and the aerodynamical resistances
          for each layer.  The values are normalized for a reference height wind
          speed of 1 m/s, and are adjusted each timestep using actual reference
          height wind speeds */
          CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
              (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
              (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
              &((*VType)[i].RaSnow));
      }

      //the third vegetation(OPEN EVERGREEN BROADLEAVED FOREST)
      if (i == 4) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= understory_monalb; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[veg_description]))
              ReportError(KeyName[veg_description], 51);
          strcpy((*VType)[i].Desc, VarStr[veg_description]);
          MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                                  uppercase and removed spaces so
                                                  it is easier to compare */
          if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
              (*VType)[i].Index = GLACIER;
          }
          else
              (*VType)[i].Index = i;

          (*VType)[i].NVegLayers = 0;

          if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
              (*VType)[i].OverStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
              (*VType)[i].OverStory = FALSE;
          else
              ReportError(KeyName[overstory], 51);

          if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
              (*VType)[i].UnderStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
              (*VType)[i].UnderStory = FALSE;
          else
              ReportError(KeyName[understory], 51);

          Veg->NLayers[i] = (*VType)[i].NVegLayers;
          if ((*VType)[i].NVegLayers > Veg->MaxLayers)
              Veg->MaxLayers = (*VType)[i].NVegLayers;

          if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
              ReportError(KeyName[number_of_root_zones], 51);

          if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
              ReportError(KeyName[imperv_frac], 51);
          impervious += (*VType)[i].ImpervFrac;

          if ((*VType)[i].ImpervFrac > 0) {
              if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
                  ReportError(KeyName[detention_frac], 51);
              if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
                  ReportError(KeyName[detention_decay], 51);
          }
          else {
              (*VType)[i].DetentionFrac = 0.;
              (*VType)[i].DetentionDecay = 0.;
          }

          /* allocate memory for the vegetation layers */
          if (!((*VType)[i].Fract = (float*)calloc((*VType)[i].NVegLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          if (Options->CanopyRadAtt == VARIABLE) {
              if (!((*VType)[i].HemiFract = (float*)calloc((*VType)[i].NVegLayers,
                  sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          else {
              (*VType)[i].HemiFract = NULL;
          }

          if (!((*VType)[i].Height = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMax = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMin = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MoistThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].VpdThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Rpc = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Albedo = (float*)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MaxInt = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))

              ReportError((char*)Routine, 1);
          if (!((*VType)[i].LAI = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RootFract = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);

          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].RootFract[j] =
                  (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          if (!((*VType)[i].RootDepth = (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].LAIMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].LAIMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          if (!((*VType)[i].AlbedoMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].AlbedoMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          /* assign the entries to the appropriate variables */
          /* allocation of zero memory is not supported on some
          compilers */
          if ((*VType)[i].OverStory == TRUE) {
              if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
                  ReportError(KeyName[fraction], 51);
              (*VType)[i].Fract[0] = Anovapara[68]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Fract[0]:%lf\n", i, (*VType)[i].Fract[0]);

              if (Options->CanopyRadAtt == VARIABLE) {
                  if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
                      ReportError(KeyName[hemifraction], 51);
                  if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
                      ReportError(KeyName[clumping_factor], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
                      ReportError(KeyName[leaf_angle_a], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
                      ReportError(KeyName[leaf_angle_b], 51);
                  if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
                      ReportError(KeyName[scat], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
              }

              else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
                  if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
                      ReportError(KeyName[beam_attn], 51);
                  (*VType)[i].Atten = Anovapara[69]; //replace parameter in configfile by anova
                  printf("(*VType)[%d].Atten:%lf\n", i, (*VType)[i].Atten);

                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
                  (*VType)[i].Taud = NOT_APPLICABLE;
              }
              else if (Options->ImprovRadiation == TRUE) {
                  if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
                      ReportError(KeyName[diff_attn], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
              }

              if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
                  ReportError(KeyName[trunk_space], 51);
              (*VType)[i].Trunk = Anovapara[70]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Trunk:%lf\n", i, (*VType)[i].Trunk);

              if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
                  ReportError(KeyName[aerodynamic_att], 51);
              (*VType)[i].Cn = Anovapara[71]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Cn:%lf\n", i, (*VType)[i].Cn);


              if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
                  ReportError(KeyName[snow_int_cap], 51);

              if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
                  ReportError(KeyName[mass_drip_ratio], 51);

              if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
                  ReportError(KeyName[snow_int_eff], 51);

              sprintf(VarStr[overstory_fraction], "%lf %lf %lf", Anovapara[72], 2 * Anovapara[72], 1 - (3 * Anovapara[72])); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[overstory_fraction], 51);
              printf("(*VType)[%d].RootFract[0]:%lf\n", i, *(*VType)[i].RootFract[0]);

              sprintf(VarStr[overstory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 1.13 * Anovapara[73], 1.18 * Anovapara[73],
                  1.27 * Anovapara[73], 1.3 * Anovapara[73], 1.34 * Anovapara[73], 1.4 * Anovapara[73], 1.35 * Anovapara[73], 1.32 * Anovapara[73],
                  1.3 * Anovapara[73], 1.27 * Anovapara[73], 1.2 * Anovapara[73], Anovapara[73]); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
                  ReportError(KeyName[overstory_monlai], 51);
              printf("(*VType)[%d].LAIMonthly[0]:%lf\n", i, *(*VType)[i].LAIMonthly[0]);

              maxLAI = -9999;
              for (k = 0; k < 12; k++) {
                  if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
                      maxLAI = (*VType)[i].LAIMonthly[0][k];
              }

              sprintf(VarStr[overstory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[74], Anovapara[74], Anovapara[74],
                  Anovapara[74], Anovapara[74], Anovapara[74], Anovapara[74], Anovapara[74], Anovapara[74], Anovapara[74], Anovapara[74],
                  Anovapara[74]); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
                  ReportError(KeyName[overstory_monalb], 51);
              printf("(*VType)[%d].AlbedoMonthly[0]:%lf\n", i, *(*VType)[i].AlbedoMonthly[0]);

              if ((*VType)[i].UnderStory == TRUE) { //UnderStory = TRUE.
                  (*VType)[i].Fract[1] = 1.0;
                  sprintf(VarStr[understory_fraction], "%lf %lf %lf", Anovapara[75], 1 - Anovapara[75], 0.00); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);
                  printf("(*VType)[%d].RootFract[1]:%lf\n", i, *(*VType)[i].RootFract[1]);

                  sprintf(VarStr[understory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[76], 1.28 * Anovapara[76], 1.67 * Anovapara[76],
                      1.67 * Anovapara[76], 1.67 * Anovapara[76], 3.14 * Anovapara[76], 3 * Anovapara[76], 1.78 * Anovapara[76], Anovapara[76], Anovapara[76],
                      Anovapara[76], Anovapara[76]); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);
                  printf("(*VType)[%d].LAIMonthly[1]:%lf\n", i, *(*VType)[i].LAIMonthly[1]);

                  sprintf(VarStr[understory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[77], Anovapara[77], Anovapara[77],
                      Anovapara[77], Anovapara[77], Anovapara[77], Anovapara[77], Anovapara[77], Anovapara[77], Anovapara[77], Anovapara[77],
                      Anovapara[77]); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
                  printf("(*VType)[%d].AlbedoMonthly[1]:%lf\n", i, *(*VType)[i].AlbedoMonthly[1]);
              }
          }
          else {
              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[0] = 1.0;
                  if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);

                  if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);

                  if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
              }
              (*VType)[i].Trunk = NOT_APPLICABLE;
              (*VType)[i].Cn = NOT_APPLICABLE;
              (*VType)[i].Atten = NOT_APPLICABLE;
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
          }

          sprintf(VarStr[height], "%lf %lf", Anovapara[78], Anovapara[79]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
              ReportError(KeyName[height], 51);
          printf("(*VType)[%d].Height.overstory:%lf\n", i, *(*VType)[i].Height);
          printf("(*VType)[%d].Height.understory:%lf\n", i, *((*VType)[i].Height + 1));

          sprintf(VarStr[max_resistance], "%lf %lf", Anovapara[80], Anovapara[81]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[max_resistance], 51);
          //printf("(*VType)[%d].RsMax:%lf\n", i, *(*VType)[i].RsMax);
          printf("(*VType)[%d].RsMax.overstory:%lf\n", i, *(*VType)[i].RsMax);
          printf("(*VType)[%d].RsMax.understory:%lf\n", i, *((*VType)[i].RsMax + 1));

          sprintf(VarStr[min_resistance], "%lf %lf", Anovapara[82], Anovapara[83]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[min_resistance], 51);
          printf("(*VType)[%d].RsMin.overstory:%lf\n", i, *(*VType)[i].RsMin);
          printf("(*VType)[%d].RsMin.understory:%lf\n", i, *((*VType)[i].RsMin + 1));

          sprintf(VarStr[moisture_threshold], "%lf %lf", Anovapara[84], Anovapara[84]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[moisture_threshold], 51);
          printf("(*VType)[%d].MoistThres:%lf\n", i, *(*VType)[i].MoistThres);

          if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
              ReportError(KeyName[vpd], 51);

          if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
              ReportError(KeyName[rpc], 51);

          sprintf(VarStr[root_zone_depth], "%lf %lf %lf", Anovapara[85], Anovapara[85], Anovapara[85]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
              (*VType)[i].NSoilLayers))
              ReportError(KeyName[root_zone_depth], 51);
          printf("(*VType)[%d].RootDepth:%lf\n", i, *(*VType)[i].RootDepth);

          /* Calculate the wind speed profiles and the aerodynamical resistances
          for each layer.  The values are normalized for a reference height wind
          speed of 1 m/s, and are adjusted each timestep using actual reference
          height wind speeds */
          CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
              (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
              (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
              &((*VType)[i].RaSnow));
      }

      //the fourth vegetation (Closed evergreen needle-leaved forest)
      if (i == 8) {
          /* Read the key-entry pairs from the input file */
          for (j = 0; j <= understory_monalb; j++) {
              sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
              GetInitString(SectionName, KeyName[j], "", VarStr[j],
                  (unsigned long)BUFSIZE, Input);
          }

          /* Assign the entries to the appropriate variables */
          if (IsEmptyStr(VarStr[veg_description]))
              ReportError(KeyName[veg_description], 51);
          strcpy((*VType)[i].Desc, VarStr[veg_description]);
          MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                                  uppercase and removed spaces so
                                                  it is easier to compare */
          if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
              (*VType)[i].Index = GLACIER;
          }
          else
              (*VType)[i].Index = i;

          (*VType)[i].NVegLayers = 0;

          if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
              (*VType)[i].OverStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
              (*VType)[i].OverStory = FALSE;
          else
              ReportError(KeyName[overstory], 51);

          if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
              (*VType)[i].UnderStory = TRUE;
              ((*VType)[i].NVegLayers)++;
          }
          else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
              (*VType)[i].UnderStory = FALSE;
          else
              ReportError(KeyName[understory], 51);

          Veg->NLayers[i] = (*VType)[i].NVegLayers;
          if ((*VType)[i].NVegLayers > Veg->MaxLayers)
              Veg->MaxLayers = (*VType)[i].NVegLayers;

          if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
              ReportError(KeyName[number_of_root_zones], 51);

          if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
              ReportError(KeyName[imperv_frac], 51);
          impervious += (*VType)[i].ImpervFrac;

          if ((*VType)[i].ImpervFrac > 0) {
              if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
                  ReportError(KeyName[detention_frac], 51);
              if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
                  ReportError(KeyName[detention_decay], 51);
          }
          else {
              (*VType)[i].DetentionFrac = 0.;
              (*VType)[i].DetentionDecay = 0.;
          }

          /* allocate memory for the vegetation layers */
          if (!((*VType)[i].Fract = (float*)calloc((*VType)[i].NVegLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);

          if (Options->CanopyRadAtt == VARIABLE) {
              if (!((*VType)[i].HemiFract = (float*)calloc((*VType)[i].NVegLayers,
                  sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          else {
              (*VType)[i].HemiFract = NULL;
          }

          if (!((*VType)[i].Height = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMax = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RsMin = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MoistThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].VpdThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Rpc = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].Albedo = (float*)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].MaxInt = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))

              ReportError((char*)Routine, 1);
          if (!((*VType)[i].LAI = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].RootFract = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);

          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].RootFract[j] =
                  (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }
          if (!((*VType)[i].RootDepth = (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
              ReportError((char*)Routine, 1);

          if (!((*VType)[i].LAIMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].LAIMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          if (!((*VType)[i].AlbedoMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
              ReportError((char*)Routine, 1);
          for (j = 0; j < (*VType)[i].NVegLayers; j++) {
              if (!((*VType)[i].AlbedoMonthly[j] = (float*)calloc(12, sizeof(float))))
                  ReportError((char*)Routine, 1);
          }

          /* assign the entries to the appropriate variables */
          /* allocation of zero memory is not supported on some
          compilers */
          if ((*VType)[i].OverStory == TRUE) {
              if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
                  ReportError(KeyName[fraction], 51);
              (*VType)[i].Fract[0] = Anovapara[86]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Fract[0]:%lf\n", i, (*VType)[i].Fract[0]);

              if (Options->CanopyRadAtt == VARIABLE) {
                  if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
                      ReportError(KeyName[hemifraction], 51);
                  if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
                      ReportError(KeyName[clumping_factor], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
                      ReportError(KeyName[leaf_angle_a], 51);
                  if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
                      ReportError(KeyName[leaf_angle_b], 51);
                  if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
                      ReportError(KeyName[scat], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
              }

              else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
                  if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
                      ReportError(KeyName[beam_attn], 51);
                  (*VType)[i].Atten = Anovapara[87]; //replace parameter in configfile by anova
                  printf("(*VType)[%d].Atten:%lf\n", i, (*VType)[i].Atten);

                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
                  (*VType)[i].Taud = NOT_APPLICABLE;
              }
              else if (Options->ImprovRadiation == TRUE) {
                  if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
                      ReportError(KeyName[diff_attn], 51);
                  (*VType)[i].Atten = NOT_APPLICABLE;
                  (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
                  (*VType)[i].Scat = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleA = NOT_APPLICABLE;
                  (*VType)[i].LeafAngleB = NOT_APPLICABLE;
              }

              if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
                  ReportError(KeyName[trunk_space], 51);
              (*VType)[i].Trunk = Anovapara[88]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Trunk:%lf\n", i, (*VType)[i].Trunk);

              if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
                  ReportError(KeyName[aerodynamic_att], 51);
              (*VType)[i].Cn = Anovapara[89]; //replace parameter in configfile by anova
              printf("(*VType)[%d].Cn:%lf\n", i, (*VType)[i].Cn);

              if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
                  ReportError(KeyName[snow_int_cap], 51);

              if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
                  ReportError(KeyName[mass_drip_ratio], 51);

              if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
                  ReportError(KeyName[snow_int_eff], 51);

              sprintf(VarStr[overstory_fraction], "%lf %lf %lf", Anovapara[90], 2 * Anovapara[90], 1 - (3 * Anovapara[90])); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[overstory_fraction], 51);
              printf("(*VType)[%d].RootFract[0]:%lf\n", i, *(*VType)[i].RootFract[0]);

              sprintf(VarStr[overstory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 1.13 * Anovapara[91], 1.18 * Anovapara[91],
                  1.27 * Anovapara[91], 1.3 * Anovapara[91], 1.34 * Anovapara[91], 1.4 * Anovapara[91], 1.35 * Anovapara[91], 1.32 * Anovapara[91],
                  1.3 * Anovapara[91], 1.27 * Anovapara[91], 1.2 * Anovapara[91], Anovapara[91]); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
                  ReportError(KeyName[overstory_monlai], 51);
              printf("(*VType)[%d].LAIMonthly[0]:%lf\n", i, *(*VType)[i].LAIMonthly[0]);

              maxLAI = -9999;
              for (k = 0; k < 12; k++) {
                  if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
                      maxLAI = (*VType)[i].LAIMonthly[0][k];
              }

              sprintf(VarStr[overstory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[92], Anovapara[92], Anovapara[92],
                  Anovapara[92], Anovapara[92], Anovapara[92], Anovapara[92], Anovapara[92], Anovapara[92], Anovapara[92], Anovapara[92],
                  Anovapara[92]); //read para will be sensitivity analysised.
              if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
                  ReportError(KeyName[overstory_monalb], 51);
              printf("(*VType)[%d].AlbedoMonthly[0]:%lf\n", i, *(*VType)[i].AlbedoMonthly[0]);

              if ((*VType)[i].UnderStory == TRUE) { //UnderStory = TRUE.
                  (*VType)[i].Fract[1] = 1.0;
                  sprintf(VarStr[understory_fraction], "%lf %lf %lf", Anovapara[93], 1 - Anovapara[93], 0.00); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);
                  printf("(*VType)[%d].RootFract[1]:%lf\n", i, *(*VType)[i].RootFract[1]);

                  sprintf(VarStr[understory_monlai], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[94], 1.28 * Anovapara[94], 1.67 * Anovapara[94],
                      1.67 * Anovapara[94], 1.67 * Anovapara[94], 3.14 * Anovapara[94], 3 * Anovapara[94], 1.78 * Anovapara[94], Anovapara[94], Anovapara[94],
                      Anovapara[94], Anovapara[94]); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);
                  printf("(*VType)[%d].LAIMonthly[1]:%lf\n", i, *(*VType)[i].LAIMonthly[1]);

                  sprintf(VarStr[understory_monalb], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Anovapara[95], Anovapara[95], Anovapara[95],
                      Anovapara[95], Anovapara[95], Anovapara[95], Anovapara[95], Anovapara[95], Anovapara[95], Anovapara[95], Anovapara[95],
                      Anovapara[95]); //read para will be sensitivity analysised.
                  if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
                  printf("(*VType)[%d].AlbedoMonthly[1]:%lf\n", i, *(*VType)[i].AlbedoMonthly[1]);
              }
          }
          else {
              if ((*VType)[i].UnderStory == TRUE) {
                  (*VType)[i].Fract[0] = 1.0;
                  if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
                      (*VType)[i].NSoilLayers))
                      ReportError(KeyName[understory_fraction], 51);

                  if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
                      ReportError(KeyName[understory_monlai], 51);

                  if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
                      ReportError(KeyName[understory_monalb], 51);
              }
              (*VType)[i].Trunk = NOT_APPLICABLE;
              (*VType)[i].Cn = NOT_APPLICABLE;
              (*VType)[i].Atten = NOT_APPLICABLE;
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
          }

          sprintf(VarStr[height], "%lf %lf", Anovapara[96], Anovapara[97]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
              ReportError(KeyName[height], 51);
          printf("(*VType)[%d].Height.overstory:%lf\n", i, *(*VType)[i].Height);
          printf("(*VType)[%d].Height.understory:%lf\n", i, *((*VType)[i].Height + 1));

          sprintf(VarStr[max_resistance], "%lf %lf", Anovapara[98], Anovapara[99]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[max_resistance], 51);
         // printf("(*VType)[%d].RsMax:%lf\n", i, *(*VType)[i].RsMax);
          printf("(*VType)[%d].RsMax.overstory:%lf\n", i, *(*VType)[i].RsMax);
          printf("(*VType)[%d].RsMax.understory:%lf\n", i, *((*VType)[i].RsMax + 1));

          sprintf(VarStr[min_resistance], "%lf %lf", Anovapara[100], Anovapara[101]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[min_resistance], 51);
          printf("(*VType)[%d].RsMin.overstory:%lf\n", i, *(*VType)[i].RsMin);
          printf("(*VType)[%d].RsMin.understory:%lf\n", i, *((*VType)[i].RsMin + 1));

          sprintf(VarStr[moisture_threshold], "%lf %lf", Anovapara[102], Anovapara[102]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
              (*VType)[i].NVegLayers))
              ReportError(KeyName[moisture_threshold], 51);
          printf("(*VType)[%d].MoistThres:%lf\n", i, *(*VType)[i].MoistThres);

          if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
              ReportError(KeyName[vpd], 51);

          if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
              ReportError(KeyName[rpc], 51);

          sprintf(VarStr[root_zone_depth], "%lf %lf %lf", Anovapara[103], Anovapara[103], Anovapara[103]); //read para will be sensitivity analysised.
          if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
              (*VType)[i].NSoilLayers))
              ReportError(KeyName[root_zone_depth], 51);
          printf("(*VType)[%d].RootDepth:%lf\n", i, *(*VType)[i].RootDepth);

          /* Calculate the wind speed profiles and the aerodynamical resistances
          for each layer.  The values are normalized for a reference height wind
          speed of 1 m/s, and are adjusted each timestep using actual reference
          height wind speeds */
          CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
              (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
              (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
              &((*VType)[i].RaSnow));
      }
      
      else {
      /* Read the key-entry pairs from the input file */
      for (j = 0; j <= understory_monalb; j++) {
          sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
          GetInitString(SectionName, KeyName[j], "", VarStr[j],
              (unsigned long)BUFSIZE, Input);
      }

      /* Assign the entries to the appropriate variables */
      if (IsEmptyStr(VarStr[veg_description]))
          ReportError(KeyName[veg_description], 51);
      strcpy((*VType)[i].Desc, VarStr[veg_description]);
      MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
                                              uppercase and removed spaces so
                                              it is easier to compare */
      if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
          (*VType)[i].Index = GLACIER;
      }
      else
          (*VType)[i].Index = i;

      (*VType)[i].NVegLayers = 0;

      if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
          (*VType)[i].OverStory = TRUE;
          ((*VType)[i].NVegLayers)++;
      }
      else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
          (*VType)[i].OverStory = FALSE;
      else
          ReportError(KeyName[overstory], 51);

      if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
          (*VType)[i].UnderStory = TRUE;
          ((*VType)[i].NVegLayers)++;
      }
      else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
          (*VType)[i].UnderStory = FALSE;
      else
          ReportError(KeyName[understory], 51);

      Veg->NLayers[i] = (*VType)[i].NVegLayers;
      if ((*VType)[i].NVegLayers > Veg->MaxLayers)
          Veg->MaxLayers = (*VType)[i].NVegLayers;

      if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
          ReportError(KeyName[number_of_root_zones], 51);

      if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
          ReportError(KeyName[imperv_frac], 51);
      impervious += (*VType)[i].ImpervFrac;

      if ((*VType)[i].ImpervFrac > 0) {
          if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
              ReportError(KeyName[detention_frac], 51);
          if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
              ReportError(KeyName[detention_decay], 51);
      }
      else {
          (*VType)[i].DetentionFrac = 0.;
          (*VType)[i].DetentionDecay = 0.;
      }

      /* allocate memory for the vegetation layers */
      if (!((*VType)[i].Fract = (float*)calloc((*VType)[i].NVegLayers,
          sizeof(float))))
          ReportError((char*)Routine, 1);

      if (Options->CanopyRadAtt == VARIABLE) {
          if (!((*VType)[i].HemiFract = (float*)calloc((*VType)[i].NVegLayers,
              sizeof(float))))
              ReportError((char*)Routine, 1);
      }
      else {
          (*VType)[i].HemiFract = NULL;
      }

      if (!((*VType)[i].Height = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].RsMax = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].RsMin = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].MoistThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].VpdThres = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].Rpc = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].Albedo = (float*)calloc(((*VType)[i].NVegLayers + 1), sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].MaxInt = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))

          ReportError((char*)Routine, 1);
      if (!((*VType)[i].LAI = (float*)calloc((*VType)[i].NVegLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].RootFract = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
          ReportError((char*)Routine, 1);

      for (j = 0; j < (*VType)[i].NVegLayers; j++) {
          if (!((*VType)[i].RootFract[j] =
              (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
              ReportError((char*)Routine, 1);
      }
      if (!((*VType)[i].RootDepth = (float*)calloc((*VType)[i].NSoilLayers, sizeof(float))))
          ReportError((char*)Routine, 1);

      if (!((*VType)[i].LAIMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
          ReportError((char*)Routine, 1);
      for (j = 0; j < (*VType)[i].NVegLayers; j++) {
          if (!((*VType)[i].LAIMonthly[j] = (float*)calloc(12, sizeof(float))))
              ReportError((char*)Routine, 1);
      }

      if (!((*VType)[i].AlbedoMonthly = (float**)calloc((*VType)[i].NVegLayers, sizeof(float*))))
          ReportError((char*)Routine, 1);
      for (j = 0; j < (*VType)[i].NVegLayers; j++) {
          if (!((*VType)[i].AlbedoMonthly[j] = (float*)calloc(12, sizeof(float))))
              ReportError((char*)Routine, 1);
      }

      /* assign the entries to the appropriate variables */
      /* allocation of zero memory is not supported on some
      compilers */
      if ((*VType)[i].OverStory == TRUE) {
          if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
              ReportError(KeyName[fraction], 51);
          if (Options->CanopyRadAtt == VARIABLE) {
              if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
                  ReportError(KeyName[hemifraction], 51);
              if (!CopyFloat(&((*VType)[i].ClumpingFactor), VarStr[clumping_factor], 1))
                  ReportError(KeyName[clumping_factor], 51);
              if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
                  ReportError(KeyName[leaf_angle_a], 51);
              if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
                  ReportError(KeyName[leaf_angle_b], 51);
              if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
                  ReportError(KeyName[scat], 51);
              (*VType)[i].Atten = NOT_APPLICABLE;
          }
          else if (Options->CanopyRadAtt == FIXED && Options->ImprovRadiation == FALSE) {
              if (!CopyFloat(&((*VType)[i].Atten), VarStr[beam_attn], 1))
                  ReportError(KeyName[beam_attn], 51);
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
              (*VType)[i].Scat = NOT_APPLICABLE;
              (*VType)[i].LeafAngleA = NOT_APPLICABLE;
              (*VType)[i].LeafAngleB = NOT_APPLICABLE;
              (*VType)[i].Taud = NOT_APPLICABLE;
          }
          else if (Options->ImprovRadiation == TRUE) {
              if (!CopyFloat(&((*VType)[i].Taud), VarStr[diff_attn], 1))
                  ReportError(KeyName[diff_attn], 51);
              (*VType)[i].Atten = NOT_APPLICABLE;
              (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
              (*VType)[i].Scat = NOT_APPLICABLE;
              (*VType)[i].LeafAngleA = NOT_APPLICABLE;
              (*VType)[i].LeafAngleB = NOT_APPLICABLE;
          }

          if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
              ReportError(KeyName[trunk_space], 51);

          if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
              ReportError(KeyName[aerodynamic_att], 51);

          if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
              ReportError(KeyName[snow_int_cap], 51);

          if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
              ReportError(KeyName[mass_drip_ratio], 51);

          if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
              ReportError(KeyName[snow_int_eff], 51);

          if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
              (*VType)[i].NSoilLayers))
              ReportError(KeyName[overstory_fraction], 51);

          if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
              ReportError(KeyName[overstory_monlai], 51);

          maxLAI = -9999;
          for (k = 0; k < 12; k++) {
              if ((*VType)[i].LAIMonthly[0][k] > maxLAI)
                  maxLAI = (*VType)[i].LAIMonthly[0][k];
          }

          if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb], 12))
              ReportError(KeyName[overstory_monalb], 51);

          if ((*VType)[i].UnderStory == TRUE) {
              (*VType)[i].Fract[1] = 1.0;
              if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[understory_fraction], 51);

              if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai], 12))
                  ReportError(KeyName[understory_monlai], 51);

              if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb], 12))
                  ReportError(KeyName[understory_monalb], 51);
          }
      }
      else {
          if ((*VType)[i].UnderStory == TRUE) {
              (*VType)[i].Fract[0] = 1.0;
              if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
                  (*VType)[i].NSoilLayers))
                  ReportError(KeyName[understory_fraction], 51);

              if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai], 12))
                  ReportError(KeyName[understory_monlai], 51);

              if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb], 12))
                  ReportError(KeyName[understory_monalb], 51);
          }
          (*VType)[i].Trunk = NOT_APPLICABLE;
          (*VType)[i].Cn = NOT_APPLICABLE;
          (*VType)[i].Atten = NOT_APPLICABLE;
          (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
      }

      if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
          ReportError(KeyName[height], 51);

      if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
          (*VType)[i].NVegLayers))
          ReportError(KeyName[max_resistance], 51);

      if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
          (*VType)[i].NVegLayers))
          ReportError(KeyName[min_resistance], 51);

      if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
          (*VType)[i].NVegLayers))
          ReportError(KeyName[moisture_threshold], 51);

      if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
          ReportError(KeyName[vpd], 51);

      if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
          ReportError(KeyName[rpc], 51);

      if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
          (*VType)[i].NSoilLayers))
          ReportError(KeyName[root_zone_depth], 51);

      /* Calculate the wind speed profiles and the aerodynamical resistances
      for each layer.  The values are normalized for a reference height wind
      speed of 1 m/s, and are adjusted each timestep using actual reference
      height wind speeds */
      CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
          (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
          (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
          &((*VType)[i].RaSnow));

          /* Run the improved radiation scheme in which the tree height, solar altitude and fractional coverage
          are all taken into account into the radiation calculation */
          if (Options->ImprovRadiation == TRUE) {
              if ((*VType)[i].OverStory == TRUE) {
                  if (!CopyFloat((*VType)[i].MonthlyExtnCoeff, VarStr[monextn], 12))
                      ReportError(KeyName[monextn], 51);
                  if (!CopyFloat(&((*VType)[i].VfAdjust), VarStr[vf_adj], 1))
                      ReportError(KeyName[vf_adj], 51);
                  (*VType)[i].Vf = (*VType)[i].Fract[0] * (*VType)[i].VfAdjust;
              }
              else {
                  if ((*VType)[i].UnderStory == TRUE) {
                      for (k = 0; k < 12; k++)
                          (*VType)[i].MonthlyExtnCoeff[k] = 0;
                      (*VType)[i].VfAdjust = 1.0;
                      /* assuming 100% coverage if understory=TRUE & overstory=FALSE*/
                      (*VType)[i].Vf = (*VType)[i].Fract[0] * (*VType)[i].VfAdjust;
                  }
              }
          }
  } 
  }
  if (impervious) {
    GetInitString(SectionName, "IMPERVIOUS SURFACE ROUTING FILE", "", VarStr[0],
      (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(VarStr[0]))
      ReportError("IMPERVIOUS SURFACE ROUTING FILE", 51);
    strcpy(Options->ImperviousFilePath, VarStr[veg_description]);
  }

  return NVegs;
}

