/*
 * SUMMARY:      InitMetMaps.c - Initialize meteorological maps
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize meteorological maps
 * DESCRIP-END.
 * FUNCTIONS:    InitMetMaps()
 *               InitEvapMap()
 *               InitPrecipMap()
 *               InitRadarMap()
 *               InitRadMap()
 * COMMENTS:
 * $Id: InitMetMaps.c,v 1.6 2006/10/03 22:50:22 nathalie Exp $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "rad.h"
#include "sizeofnt.h"
#include "varid.h"


 /*****************************************************************************
   InitMetMaps()
 *****************************************************************************/
void InitMetMaps(LISTPTR Input, int NDaySteps, MAPSIZE *Map, MAPSIZE *Radar,
  OPTIONSTRUCT *Options, char *WindPath, char *PrecipLapseFile,
  float ***PrecipLapseMap, float ***PrismMap,
  unsigned char ****ShadowMap, float ***SkyViewMap,
  EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap, float ***PptMultiplierMap,
  RADARPIX ***RadarMap, PIXRAD ***RadMap,
  SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
  LAYER *Veg, TOPOPIX **TopoMap, float ****MM5Input,
  float ****WindModel)
{
  int y, x;

  printf("Initializing meteorological maps\n");

  InitEvapMap(Map, EvapMap, SoilMap, Soil, VegMap, Veg, TopoMap);
  InitPrecipMap(Map, PrecipMap, VegMap, Veg, TopoMap);
  InitPptMultiplierMap(Options, Map, PptMultiplierMap);                                                            

  if (Options->MM5 == TRUE) {
    InitMM5Maps(Soil->MaxLayers, Map->NY, Map->NX, MM5Input, RadMap, Options);
    /* If called for, use the precip lapse map for MM5 precip
       distribution, avoiding lots of function interfaces changes */
    if (strlen(PrecipLapseFile) > 0) {
      InitPrecipLapseMap(PrecipLapseFile, Map, PrecipLapseMap);
    }
    if (Options->Shading == TRUE)
      InitShadeMap(Options, NDaySteps, Map, ShadowMap, SkyViewMap);
  }
  else {
    if (Options->PrecipType == RADAR)
      InitRadarMap(Radar, RadarMap);
    if (Options->PrecipLapse == MAP)
      InitPrecipLapseMap(PrecipLapseFile, Map, PrecipLapseMap);
    if (Options->Prism == TRUE)
      InitPrismMap(Map->NY, Map->NX, PrismMap);
    if (Options->Shading == TRUE)
      InitShadeMap(Options, NDaySteps, Map, ShadowMap, SkyViewMap);

    if (!((*SkyViewMap) = (float **)calloc(Map->NY, sizeof(float *))))
      ReportError("InitMetMaps()", 1);
    for (y = 0; y < Map->NY; y++) {
      if (!((*SkyViewMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
        ReportError("InitMetMaps()", 1);
    }
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        (*SkyViewMap)[y][x] = 1.0;
      }
    }
    if (Options->WindSource == MODEL)
      InitWindModelMaps(WindPath, Map, WindModel);

    InitRadMap(Map, RadMap);
  }
  if (Options->MM5 == TRUE && Options->QPF == TRUE && Options->Prism == TRUE)
    InitPrismMap(Map->NY, Map->NX, PrismMap);
}


/*****************************************************************************
  InitEvapMap()
*****************************************************************************/
void InitEvapMap(MAPSIZE *Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap,
  LAYER *Soil, VEGPIX **VegMap, LAYER *Veg,
  TOPOPIX **TopoMap)
{
  const char *Routine = "InitEvapMap";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */

  if (DEBUG)
    printf("Initializing evaporation map\n");

  if (!(*EvapMap = (EVAPPIX **)calloc(Map->NY, sizeof(EVAPPIX *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*EvapMap)[y] = (EVAPPIX *)calloc(Map->NX, sizeof(EVAPPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
        assert(VegMap[y][x].Veg > 0 && SoilMap[y][x].Soil > 0);

        if (!((*EvapMap)[y][x].EPot =
          (float *)calloc(NVeg + 1, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].EAct =
          (float *)calloc(NVeg + 1, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].EInt = (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*EvapMap)[y][x].ESoil =
          (float **)calloc(NVeg, sizeof(float *))))
          ReportError((char *)Routine, 1);

        for (i = 0; i < NVeg; i++) {
          if (!((*EvapMap)[y][x].ESoil[i] =
            (float *)calloc(NSoil, sizeof(float))))
            ReportError((char *)Routine, 1);
        }
      }
    }
  }
}

/*****************************************************************************
  InitPrecipMap()
*****************************************************************************/
void InitPrecipMap(MAPSIZE * Map, PRECIPPIX *** PrecipMap, VEGPIX ** VegMap,
  LAYER * Veg, TOPOPIX ** TopoMap)
{
  const char *Routine = "InitPrecipMap";
  int x;			/* counter */
  int y;			/* counter */
  int NVeg;			/* Number of veg layers at current pixel */

  if (DEBUG)
    printf("Initializing precipitation map\n");

  if (!(*PrecipMap = (PRECIPPIX **)calloc(Map->NY, sizeof(PRECIPPIX *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*PrecipMap)[y] = (PRECIPPIX *)calloc(Map->NX, sizeof(PRECIPPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        if (!((*PrecipMap)[y][x].IntRain =
          (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
        if (!((*PrecipMap)[y][x].IntSnow =
          (float *)calloc(NVeg, sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        (*PrecipMap)[y][x].PrecipStart = TRUE;
    }
  }
}

/*******************************************************************************
  InitMM5Maps()
*******************************************************************************/
void InitMM5Maps(int NSoilLayers, int NY, int NX, float ****MM5Input,
  PIXRAD ***RadMap, OPTIONSTRUCT *Options)
{
  char *Routine = "InitMM5Maps";
  int NTotalMaps = NSoilLayers + N_MM5_MAPS;
  int n;
  int y;

  if (Options->HeatFlux == FALSE)
    NTotalMaps -= NSoilLayers;

  if (!((*MM5Input) = (float ***)calloc(NTotalMaps, sizeof(float **))))
    ReportError(Routine, 1);

  for (n = 0; n < NTotalMaps; n++) {
    if (!((*MM5Input)[n] = (float **)calloc(NY, sizeof(float *))))
      ReportError(Routine, 1);
    for (y = 0; y < NY; y++) {
      if (!((*MM5Input)[n][y] = (float *)calloc(NX, sizeof(float))))
        ReportError(Routine, 1);
    }
  }

  /* Initiate radiation map */
  if (!(*RadMap = (PIXRAD **)calloc(NY, sizeof(PIXRAD *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < NY; y++) {
    if (!((*RadMap)[y] = (PIXRAD *)calloc(NX, sizeof(PIXRAD))))
      ReportError((char *)Routine, 1);
  }
}

/*******************************************************************************
  InitWindModelMaps()
*******************************************************************************/
void InitWindModelMaps(char *WindPath, MAPSIZE *Map, float ****WindModel)
{
  char *Routine = "InitWindModelMaps";
  char InFileName[NAMESIZE + 1];
  char Str[NAMESIZE + 1];
  int NumberType;
  int n;
  int x;
  int y;
  float *Array = NULL;

  if (!((*WindModel) = (float ***)calloc(NWINDMAPS, sizeof(float **))))
    ReportError(Routine, 1);

  for (n = 0; n < NWINDMAPS; n++) {
    if (!((*WindModel)[n] = (float **)calloc(Map->NY, sizeof(float *))))
      ReportError(Routine, 1);
    for (y = 0; y < Map->NY; y++) {
      if (!((*WindModel)[n][y] = (float *)calloc(Map->NX, sizeof(float))))
        ReportError(Routine, 1);
    }
  }

  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);
  NumberType = NC_FLOAT;

  /* Read the wind model maps */
  for (n = 0; n < NWINDMAPS; n++) {
    sprintf(Str, "%02d", n + 1);
    sprintf(InFileName, "%s%s%s", WindPath, Str, fileext);
    Read2DMatrix(InFileName, Array, NumberType, Map, 0, "", 0);
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        (*WindModel)[n][y][x] = Array[y * Map->NX + x];
      }
    }
  }
  free(Array);
}

/*******************************************************************************
  InitRadarMap()
*******************************************************************************/
void InitRadarMap(MAPSIZE *Radar, RADARPIX ***RadarMap)
{
  const char *Routine = "InitRadarMap";
  int y;			/* counter */

  if (DEBUG)
    printf("Initializing radar precipitation map\n");

  if (!(*RadarMap = (RADARPIX **)calloc(Radar->NY, sizeof(RADARPIX *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Radar->NY; y++) {
    if (!((*RadarMap)[y] = (RADARPIX *)calloc(Radar->NX, sizeof(RADARPIX))))
      ReportError((char *)Routine, 1);
  }
}

/******************************************************************************
  InitRadMap()
******************************************************************************/
void InitRadMap(MAPSIZE *Map, PIXRAD ***RadMap)
{
  const char *Routine = "InitRadMap";
  int y;			/* counter */

  if (DEBUG)
    printf("Initializing radiation map\n");

  if (!(*RadMap = (PIXRAD **)calloc(Map->NY, sizeof(PIXRAD *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*RadMap)[y] = (PIXRAD *)calloc(Map->NX, sizeof(PIXRAD))))
      ReportError((char *)Routine, 1);
  }
}

/******************************************************************************/
/*			       InitPrecipLapseMap                             */
/******************************************************************************/
void InitPrecipLapseMap(char *PrecipLapseFile, MAPSIZE *Map, float ***PrecipLapseMap)
{
  const char *Routine = "InitPrecipLapseMap";
  int NumberType;
  int x;			/* counter */
  int y;			/* counter */
  float *Array = NULL;

  if (!((*PrecipLapseMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*PrecipLapseMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }

  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);
  NumberType = NC_FLOAT;

  Read2DMatrix(PrecipLapseFile, Array, NumberType, Map, 0, "", 0);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*PrecipLapseMap)[y][x] = Array[y * Map->NX + x];
    }
  }

  free(Array);
}

/******************************************************************************/
/*			       InitPRISMMap                             */
/******************************************************************************/

void InitPrismMap(int NY, int NX, float ***PrismMap)
{
  const char *Routine = "InitPRISMMap";
  int x;			/* counter */
  int y;			/* counter */

  if (!((*PrismMap) = (float **)calloc(NY, sizeof(float *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < NY; y++) {
    if (!((*PrismMap)[y] = (float *)calloc(NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      (*PrismMap)[y][x] = 1.0;
    }
  }

}

/******************************************************************************/
/*				  InitShadeMap                                */
/******************************************************************************/
void InitShadeMap(OPTIONSTRUCT * Options, int NDaySteps, MAPSIZE *Map,
  unsigned char ****ShadowMap, float ***SkyViewMap)
{
  const char *Routine = "InitShadeMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int x;			/* counter */
  int y;			/* counter */
  int n;
  int NumberType;
  float *Array = NULL;

  if (!((*ShadowMap) =
    (unsigned char ***)calloc(NDaySteps, sizeof(unsigned char **))))
    ReportError((char *)Routine, 1);
  for (n = 0; n < NDaySteps; n++) {
    if (!((*ShadowMap)[n] =
      (unsigned char **)calloc(Map->NY, sizeof(unsigned char *))))
      ReportError((char *)Routine, 1);
    for (y = 0; y < Map->NY; y++) {
      if (!((*ShadowMap)[n][y] =
        (unsigned char *)calloc(Map->NX, sizeof(unsigned char))))
        ReportError((char *)Routine, 1);
    }
  }

  if (!((*SkyViewMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SkyViewMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SkyViewMap)[y][x] = 1.0;
    }
  }

  GetVarName(305, 0, VarName);
  GetVarNumberType(305, &NumberType);
  if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *)Routine, 1);
  Read2DMatrix(Options->SkyViewDataPath, Array, NumberType, Map, 0,
    VarName, 0);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      (*SkyViewMap)[y][x] = Array[y * Map->NX + x];
    }
  }

  free(Array);
}

/*****************************************************************************
InitPptMultiplierMap()
*****************************************************************************/
void InitPptMultiplierMap(OPTIONSTRUCT * Options, MAPSIZE *Map, float ***PptMultiplierMap)
{
  const char *Routine = "InitPptMultiplierMap";
  char VarName[BUFSIZE + 1];
  int i, j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  float *Array;


  /* Get the canopy gap map filename from the [VEGETATION] section */
  if (!((*PptMultiplierMap) = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*PptMultiplierMap)[y] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
  }
  if (PRECIP_MULTIPLIER > NA) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*PptMultiplierMap)[y][x] = PRECIP_MULTIPLIER;
      }
    }
  }
  else if (!IsEmptyStr(Options->PrecipMultiplierMapPath)) {
    /* Read the map path */
    GetVarName(100, 0, VarName);
    GetVarNumberType(100, &NumberType);
    if (!(Array = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    Read2DMatrix(Options->PrecipMultiplierMapPath, Array, NumberType, Map, 0, VarName, 0);

    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*PptMultiplierMap)[y][x] = Array[i];
      }
    }
    free(Array);
  }
  else {
    printf("No valid input of precipitation multiplier ...\n");
    exit(88);
  }

}                                                               