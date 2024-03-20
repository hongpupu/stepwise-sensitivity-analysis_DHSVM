/*
* SUMMARY:      EvapoTranspiration.c - Calculate evapotranspiration
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Calculate evapotranspiration
* DESCRIP-END.
* FUNCTIONS:    EvapoTranspiration()
* COMMENTS:
* $Id: EvapoTranspiration.c,v 1.5 2007/03/02 22:02:01 lancuo Exp $
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "functions.h"

/*****************************************************************************
EvapoTranspiration()
*****************************************************************************/
void EvapoTranspiration(int Layer, int ImpvRad, int Dt, float F, PIXMET *Met,
  float NetRad, float Rp, VEGTABLE *VType, SOILTABLE *SType,
  float MoistureFlux, float *Moist, float *SoilTemp, float *Int,
  float *EPot, float *EInt, float **ESoil, float *EAct, float *ETot,
  float *Adjust, float Ra, VEGPIX *LocalVeg)
{
  float *Rc;			/* canopy resistance associated with
                        conditions in each soil layer (s/m) */
  float DryArea;		/* relative dry leaf area  */
  float DryEvapTime;	/* amount of time remaining during a timestep
                        after the interception storage is depleted (sec) */
  //float F;			    /* Fractional coverage by vegetation layer */
  float SoilMoisture;	/* Amount of water in each soil layer (m) */
  float WetArea;		/* relative leaf area wetted by interception storage */
  float WetEvapRate;	/* evaporation rate from wetted fraction per unit ground area (m/s) */
  float WetEvapTime;	/* amount of time needed to evaporate the amount of water
                        in interception storage (sec) */
  int i;			    /* counter */

  /* make F an input variable */
  //F = LocalVeg->Fract[Layer];
  
  /* removing this line as Net Radiation already accounted for F in radiation calculation */
  NetRad /= F;

  /* Convert the water amounts related to partial canopy cover to a pixel depth
  as if the entire pixel is covered. These depths will be converted back later on. */
  *Int /= F;
  MoistureFlux /= F;
  LocalVeg->MaxInt[Layer] /= F;

  /* allocate memory for the canopy resistance array */
  if (!(Rc = (float *)calloc(VType->NSoilLayers, sizeof(float))))
    ReportError("EvapoTranspiration()", 1);

  /* Calculate the evaporation rate in m/s */
  EPot[Layer] = (Met->Slope * NetRad + Met->AirDens * CP * Met->Vpd / Ra) /
    (WATER_DENSITY * Met->Lv * (Met->Slope + Met->Gamma));

  /* The potential evaporation rate accounts for the amount of moisture that
  the atmosphere can absorb. If we do not account for the amount of
  evaporation from overlying evaporation, we can end up with the situation
  that all vegetation layers and the soil layer transpire/evaporate at the
  potential rate, resulting in an overprediction of the actual evaporation
  rate. Thus we subtract the amount of evaporation that has already
  been calculated for overlying layers from the potential evaporation.
  Another mechanism that could be used to account for this would be to
  decrease the vapor pressure deficit while going down through the canopy
  (not implemented here) */
  EPot[Layer] -= MoistureFlux / Dt;

  if (EPot[Layer] < 0)
    EPot[Layer] = 0;

  /* WetArea = pow(*Int/VType->MaxInt[Layer], (double) 2.0/3.0); */

  WetArea = cbrt(*Int / LocalVeg->MaxInt[Layer]);
  WetArea = MIN(WetArea, 1);
  WetArea = WetArea * WetArea;
  DryArea = 1 - WetArea;

  /* calculate the amount of water that can evaporate from the interception
  storage.  Given this evaporation rate, calculate the amount of time it
  will take to evaporate the entire amount of intercepted water.  If this
  time period is shorter than the length of the time interval, the
  previously wetted leaves can transpire during the remaining part of the
  time interval */

  /* WORK IN PROGRESS:  the amount of interception storage can be replenished
  during (low intensity) rainstorms */

  WetEvapRate = WetArea * EPot[Layer];
  if (WetEvapRate > 0) {
    WetEvapTime = *Int / WetEvapRate;
    if (WetEvapTime > Dt) {
      WetEvapTime = Dt;
    }
  }
  else if (*Int > 0)
    WetEvapTime = Dt;
  else
    WetEvapTime = 0;

  if (WetEvapRate > 0) {
    if (WetEvapTime < Dt) {
      EInt[Layer] = *Int;
      *Int = 0.0;
      DryEvapTime = Dt - WetEvapTime;
    }
    else {
      EInt[Layer] = Dt * WetEvapRate;
      *Int -= EInt[Layer];
      DryEvapTime = 0.0;
    }
  }
  else {
    EInt[Layer] = 0.0;
    if (*Int > 0)
      DryEvapTime = 0.0;
    else
      DryEvapTime = Dt;
  }

  /* Correct the evaporation from interception and the interception storage for
  the fractional overstory coverage */
  EInt[Layer] *= F;
  *ETot += EInt[Layer];
  *Int *= F;
  LocalVeg->MaxInt[Layer] *= F;

  /* calculate the canopy conductances associated with the conditions in
  each of the soil layers */
  for (i = 0; i < VType->NSoilLayers; i++)
    Rc[i] = CanopyResistance(LocalVeg->LAI[Layer], VType->RsMin[Layer],
      VType->RsMax[Layer], VType->Rpc[Layer],
      VType->VpdThres[Layer], VType->MoistThres[Layer],
      SType->WP[i], SoilTemp[i], Moist[i], Met->Vpd, Rp);

  /* calculate the transpiration rate for the current vegetation layer,
  and adjust the soil moisture content in each of the soil layers */
  for (i = 0; i < VType->NSoilLayers; i++) {
    ESoil[Layer][i] = (Met->Slope + Met->Gamma) /
      (Met->Slope + Met->Gamma * (1 + Rc[i] / Ra)) * VType->RootFract[Layer][i] *
      EPot[Layer] * Adjust[i];

    /* calculate the amounts of water transpirated during each timestep based
    on the evaporation and transpiration rates.  While there is still water
    in interception storage only the area that is not covered by intercep-
    ted water will transpire.  When all of the interception storage has
    disappeared all leaves will contribute to the transpiration */
    ESoil[Layer][i] *= WetEvapTime * (1 - WetArea) + DryEvapTime;

    SoilMoisture = Moist[i] * VType->RootDepth[i] * Adjust[i];
    if (SoilMoisture < ESoil[Layer][i])
      ESoil[Layer][i] = SoilMoisture;

    /* correct the evaporation for the fractional overstory coverage and update
    the soil moisture */
    ESoil[Layer][i] *= F;
    SoilMoisture -= ESoil[Layer][i];

    Moist[i] = SoilMoisture / (VType->RootDepth[i] * Adjust[i]);

  }

  for (i = 0, EAct[Layer] = 0; i < VType->NSoilLayers; i++) {
    *ETot += ESoil[Layer][i];
    EAct[Layer] += ESoil[Layer][i];
   // if (EAct[Layer] < 0)
   //   printf("veg[%d]+soil[%d]: EAct=%f; ESoil=%f\n", Layer, i, EAct[Layer], ESoil[Layer][i]);
  }


  /* clean up */
  free(Rc);
}

