#ifndef _DEFINITIONS_HH_
#define _DEFINITIONS_HH_

#include"libs.hh"

const UInt_t    NEVENTS             = 50; //Number of particles to propagate
const Double_t  AMU                 = 0.9314940038; //GeV 90 from AME 2016
const Int_t     BEAM_Z              = 50; //primary Z
const Int_t     BEAM_A              = 124;//Primary A
const Double_t  BEAM_MASS_GEV_C2    = 123.905279 * AMU;// Primary mass in GeV/c2 - 124Sn
const Double_t  BEAM_DTHETA_DEG     = 1e-9;//max spread of theta angle in deg (+-)
const Double_t  BEAM_DTHETA_DEG     = 0.02*180./3.14;//max spread of theta angle in deg (+-)
const Double_t  BEAM_RADIUS_CM      = 3;//radius of the beamspot in cm 
const Double_t  BEAM_EKIN_GEV_U     = 0.9;//per unit in GeV
const Double_t  SCALE_MIN           = 0.15;// scaling for subtraction: p_min = Ptot-Ptot*scale_min
const Double_t  SCALE_MAX           = 0.15; // scaling for aadition:   p_max = Ptot+Ptot*scale_max
//const Double_t  SCALE_MIN           = 1e-19;// scaling for subtraction: p_min = Ptot-Ptot*scale_min
//const Double_t  SCALE_MAX           = 1e-19; // scaling for aadition:   p_max = Ptot+Ptot*scale_max
const Double_t  GLAD_CURRENT        = -2983.0; 

//Plotting ranges in cm
const  Double_t Ymin = -100;
const  Double_t Ymax =  100;
const  Double_t Xmin = -1050;
const  Double_t Xmax =  1050;
const  Double_t Zmin = -300;
const  Double_t Zmax =  1800;
const  UInt_t NbinsX = 1000;
const  UInt_t NbinsY = 1000;
const  UInt_t NbinsZ = 1000;

const Double_t Target_to_GLAD_flange = 176.65;
const Double_t Target_to_MWPC = 164.85;
const Double_t Target_to_TP = Target_to_GLAD_flange + 164.;
const Double_t TP_to_F10 = 381.;
const Double_t TP_to_F11 = 1347.2;
const Double_t TP_to_F12 = TP_to_F11 + 20;
const Double_t TP_to_TOFD = TP_to_F12 + 14;

#endif

