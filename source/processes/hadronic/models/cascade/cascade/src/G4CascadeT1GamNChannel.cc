//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4CascadeT1GamNChannel.cc 67796 2013-03-08 06:18:39Z mkelsey $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 06 20 18 - N. Toro added new processes for pion-nucleon backscatter
//   would be more elegant to modify kin. distribution for
//   existing pi-nucleon processes to populate this phase space but
//   tricky because I don't really understand how these distributions
//   were derived, and comparing Geant4 distribution to data
//   from 10.1103/PhysRevLett.94.012003 suggests that the current
//   distributions underpopulate backscatter region already for 2.45 GeV
//  
// 06 20 18 - N. Toro reversed the ordering of pion and nucleon in 2-body
//   final states so that kinematic distribution is consistent.

#include "G4CascadeGamPChannel.hh"
#include "G4CascadeGamNChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

/* * * * * * ** * * * * * * ** * * * * * * ** * * * * * * ** * * * * * * ** * * * * * * ** 
Note that Bertini obtains cross-sections by interpolating data at the following energies:

static const G4double bins[30] = {
  0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
  0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
  2.4,  3.2,  4.2,   5.6,   7.5,   10.0,  13.0,  18.0,  24.0, 32.0 };
}

This code up-scales Kaon production processes at energies of 2.4 to 10 GeV; photons with energies between 1.8 and 2.4 GeV, or between 10 and 13 GeV, will receive a partial (reduced) upscaling.

note also that 2-body processes have zero input cross-section above 7.5 GeV -- NOT ACCEPTABLE
MODEL FOR 8 GEV STUDIES. -NT
* * * * * * ** * * * * * * ** * * * * * * ** * * * * * * ** * * * * * * ** * * * * * * **/

namespace {


  // Can safely bias 2- and 3-body simultaneously by ~60
  // or 2-body only by ~70.  At higher biases, interpolation 
  // artifacts lead to sizable energy-dependence in actual bias.
  // For 4 GeV beam, could do 2-body biasing up to ~100.
  static const G4double upK2f = 45.; 
  static const G4double upK3f = 45.; 
                                  
  static const G4double sigmaK2[6] = {0.001, 0.0007, 0.0004, 0.0003, 0., 0.};
  static const G4double sigmaK3[6] = {0.0004, 0.001, 0.001, 0.0008, 0.0007, 0.0006};
  static const G4double totalXS[6] = {0.138,  0.1296, 0.1275, 0.124,  0.122,  0.12};
  static G4double upK2(int i);
  static G4double upK3(int i);
  static G4double biasedTotal(int i);

  // gamma p : Outgoing particle types of a given multiplicity
  static const G4int gamp2bfs[8][2] =
    {{gam,pro}, {pi0, pro}, {pip, neu}, {pro,pi0}, {neu,pip}, {kpl, lam}, {kpl, s0}, {k0, sp}};

  static const G4int gamp3bfs[6][3] =
  {{pro,pip,pim}, {pro,pi0,pi0}, {neu,pip,pi0},
   {pro,kpl,kmi},  {pro,k0,k0b},  {neu,kpl,k0b}};

  static const G4int gamp4bfs[4][4] =
  {{pro,pip,pim,pi0}, {pro,pi0,pi0,pi0}, {neu,pip,pip,pim}, {neu,pip,pi0,pi0}};

  static const G4int gamp5bfs[5][5] =
  {{pro,pip,pip,pim,pim}, {pro,pip,pim,pi0,pi0}, {pro,pi0,pi0,pi0,pi0},
   {neu,pip,pip,pim,pi0}, {neu,pip,pi0,pi0,pi0}};

  static const G4int gamp6bfs[6][6] =
  {{pro,pip,pip,pim,pim,pi0}, {pro,pip,pim,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0},
   {neu,pip,pip,pip,pim,pim}, {neu,pip,pip,pim,pi0,pi0}, {neu,pip,pi0,pi0,pi0,pi0}};

  static const G4int gamp7bfs[7][7] =
  {{pro,pip,pip,pip,pim,pim,pim}, {pro,pip,pip,pim,pim,pi0,pi0},
   {pro,pip,pim,pi0,pi0,pi0,pi0}, {pro,pi0,pi0,pi0,pi0,pi0,pi0},
   {neu,pip,pip,pim,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pi0},
   {neu,pip,pi0,pi0,pi0,pi0,pi0}};

  static const G4int gamp8bfs[7][8] =
  {{pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pi0},
   {pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pi0,pi0},
   {neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {neu,pip,pip,pip,pip,pim,pim,pim},
   {neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}};

  static const G4int gamp9bfs[7][9] =
  {{pro,pip,pip,pip,pim,pim,pim,pi0,pi0}, {pro,pip,pip,pim,pim,pi0,pi0,pi0,pi0},
   {pro,pip,pip,pip,pip,pim,pim,pim,pim}, {pro,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
   {neu,pip,pip,pip,pim,pim,pi0,pi0,pi0}, {neu,pip,pip,pip,pip,pim,pim,pim,pi0},
   {neu,pip,pip,pim,pi0,pi0,pi0,pi0,pi0}};
}

namespace {
  // gamma n : Outgoing particle types of a given multiplicity
  static const G4int gamn2bfs[8][2] =
    {{gam,neu}, {pi0,neu}, {pim, pro}, {neu,pi0}, {pro,pim}, {k0, lam}, {k0, s0}, {kpl, sm}};

  static const G4int gamn3bfs[6][3] =
  {{neu,pip,pim}, {neu,pi0,pi0}, {pro,pim,pi0},
   {neu,k0,k0b},  {neu,kpl,kmi},   {pro,k0,kmi}};

  static const G4int gamn4bfs[4][4] =
  {{neu,pip,pim,pi0}, {neu,pi0,pi0,pi0}, {pro,pip,pim,pim}, {pro,pim,pi0,pi0}};

  static const G4int gamn5bfs[5][5] =
  {{neu,pip,pip,pim,pim}, {neu,pip,pim,pi0,pi0}, {neu,pi0,pi0,pi0,pi0},
   {pro,pip,pim,pim,pi0}, {pro,pim,pi0,pi0,pi0}};

  static const G4int gamn6bfs[6][6] =
  {{neu,pip,pip,pim,pim,pi0}, {neu,pip,pim,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0},
   {pro,pip,pip,pim,pim,pim}, {pro,pip,pim,pim,pi0,pi0}, {pro,pim,pi0,pi0,pi0,pi0}};

  static const G4int gamn7bfs[7][7] =
  {{neu,pip,pip,pip,pim,pim,pim}, {neu,pip,pip,pim,pim,pi0,pi0},
   {neu,pip,pim,pi0,pi0,pi0,pi0}, {neu,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,pip,pim,pim,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pim,pi0},
   {pro,pim,pi0,pi0,pi0,pi0,pi0}};

  static const G4int gamn8bfs[7][8] =
  {{neu,pip,pip,pim,pim,pi0,pi0,pi0}, {neu,pip,pip,pip,pim,pim,pim,pi0},
   {neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pip,pip,pim,pim,pim,pi0,pi0},
   {pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pim},
   {pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}};

  static const G4int gamn9bfs[7][9] =
  {{neu,pip,pip,pip,pim,pim,pim,pi0,pi0}, {neu,pip,pip,pim,pim,pi0,pi0,pi0,pi0},
   {neu,pip,pip,pip,pip,pim,pim,pim,pim}, {neu,pip,pim,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,pip,pip,pim,pim,pim,pi0,pi0,pi0}, {pro,pip,pip,pip,pim,pim,pim,pim,pi0},
   {pro,pip,pim,pim,pi0,pi0,pi0,pi0,pi0}};
}


namespace {
  // Total gamma p cross section as a function of kinetic energy
  static const G4double gampTotXSec[30] =
   {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 
    0.0002, 0.0723, 0.2327, 0.54,   0.24,   0.19,   0.28,   0.223, 0.18, 0.146, 
    biasedTotal(0), biasedTotal(1), biasedTotal(2),
    biasedTotal(3), biasedTotal(4), biasedTotal(5),
    /*0.138,  0.1296, 0.1275, 0.124,  0.122,  0.12, */  
    0.1185, 0.117,  0.115,  0.115};

  static const G4double gampCrossSections[50][30] = {
  //
  // multiplicity 2 (6 channels)
  //
  // gamma p elastic ( = gamma n elastic )
   {0.0002,  0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
    0.0002,  0.0002, 0.0003, 0.0007, 0.0027, 0.0014, 0.0010, 0.0012, 0.0005,
    0.0003,  0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0001, 0.0001, 0.0001,
    0.0001,  0.0001, 0.0001},
 
  // pi0 p ( = pi0 n ) with new data above 500 MeV
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
    0.0,    0.01,   0.08,   0.3,    0.1143, 0.03,  0.029,  0.0229, 0.0106, 0.005,
    0.0018, 0.0009, 0.0004, 0.0002, 0.0,    0.0,   0.0,    0.0,    0.0,    0.0}, 
 
  // pi+ n ( = pi- p ) with new data above 800 MeV
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,   0.0,
    0.0,    0.062,  0.152,  0.2373, 0.1143, 0.084, 0.0983, 0.0419, 0.015, 0.006,
    0.0035, 0.0022, 0.0014, 0.0007, 0.0003, 0.0,   0.0,    0.0,    0.0,   0.0},

  // p pi0 ( = n pi0 ) backscatter -- from integral of Anderson 1969 data.  Start at 2.4 GeV 
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
   1.5e-4, 6.4e-5, 2.85e-5, 1.2e-5, 5.0e-6 , 2.1e-6, 9.6e-7,  3.6e-7,  1.5e-7, 6.4e-8}, 
 
  // n pi+ ( = p pi- ) backscatter -- from integral of Anderson 1969 data.  Start at 2.4 GeV
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
   1.5e-4, 6.4e-5, 2.85e-5, 1.2e-5, 5.0e-6 , 2.1e-6, 9.6e-7,  3.6e-7,  1.5e-7, 6.4e-8}, 

  // L K+ ( = L K0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,   0.0,  0.0,  0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,  0.0,  0.0012, 0.0019, 0.0007,
    upK2(0)*0.0004, upK2(1)*0.0003, upK2(2)*0.0002, upK2(3)*0.0001, upK2(4)*0.0,  upK2(5)*0.0,  0.0,  0.0,    0.0,    0.0}, 

  // S0 K+ ( = S0 K0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0, 0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0, 0.0012, 0.0006,
    upK2(0)*0.0003, upK2(1)*0.0002, upK2(2)*0.0001, upK2(3)*0.0001, upK2(4)*0.0, upK2(5)*0.0,  0.0,  0.0, 0.0,    0.0},
 
  // S+ K0  ( = S- K+ )
   {0.0,    0.0,    0.0,    0.0,    0.0,   0.0,  0.0,  0.0, 0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,  0.0,  0.0, 0.0014, 0.0006,
    upK2(0)*0.0003, upK2(1)*0.0002, upK2(2)*0.0001, upK2(3)*0.0001, upK2(4)*0.0, upK2(5)*0.0,  0.0,  0.0, 0.0,    0.0},
  //
  // multiplicity 3 (6 channels)
  //
  // p pi+ pi- ( = n pi+ pi- )
   {0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,   0.0,   0.0,
    0.0,   0.0,   0.0,   0.0,   0.004, 0.045, 0.072, 0.0728, 0.060, 0.041,
    0.033, 0.026, 0.022, 0.018, 0.016, 0.015, 0.014, 0.0135, 0.013, 0.012},

  // p pi0 pi0 ( = n pi0 pi0 )
   {0.0,   0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,   0.0,
    0.0,   0.0,    0.0,    0.0,    0.002, 0.009,  0.0299, 0.0265, 0.025, 0.016,
    0.011, 0.0089, 0.0076, 0.0071, 0.006, 0.0065, 0.0066, 0.007,  0.006, 0.006},

  // n pi+ pi0 ( = p pi- pi0 ) 
   {0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,     0.0,
    0.0,   0.0,   0.0,    0.0,    0.004,  0.021,  0.0457, 0.0498, 0.0418, 0.029,
    0.017, 0.013, 0.0081, 0.0043, 0.0017, 0.0006, 0.0002, 0.0,    0.0,     0.0},

  // p K+ K- ( = n K0 K0b )
   {0.0,    0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    upK3(0)*0.0004, upK3(1)*0.001, upK3(2)*0.001, upK3(3)*0.0008, upK3(4)*0.0007, upK3(5)*0.0006, 0.0005, 0.0003, 0.0002, 0.0001},   

  // p K0 K0b  (= n K+ K- )
   {0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    upK3(0)*0.0004, upK3(1)*0.001, upK3(2)*0.0009, upK3(3)*0.0006, upK3(4)*0.0007, upK3(5)*0.0006, 0.0005, 0.0003, 0.0002, 0.0001},  

  // n K+ K0b  (= p K0 K- )
   {0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    upK3(0)*0.0004, upK3(1)*0.001, upK3(2)*0.0009, upK3(3)*0.0006, upK3(4)*0.0007, upK3(5)*0.0006, 0.0005, 0.0003, 0.0002, 0.0001},
  //
  // multiplicity 4 (4 channels)
  //
  // p pi+ pi- pi0 ( = n pi+ pi- pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0021, 0.0041, 0.0132,  0.024,
    0.025,  0.0225, 0.0196, 0.0162, 0.0125, 0.0088, 0.0059, 0.003,  0.0013, 0.0004},

  // p pi0 pi0 pi0 ( = n pi0 pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0004, 0.0007, 0.0015, 0.0025,
    0.004,  0.0038, 0.0033, 0.0027, 0.0021, 0.0015, 0.0009, 0.0005, 0.0002, 0.0001},

  // n pi+ pi+ pi- ( = p pi+ pi- pi- )
   {0.0,   0.0,   0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
    0.0,   0.0,   0.0,    0.0,    0.0,    0.0,   0.0007, 0.0013, 0.0038,  0.0088,
    0.011, 0.010, 0.0081, 0.0064, 0.0046, 0.003, 0.0018, 0.0008, 0.0003, 0.0001},

  // n pi+ pi0 pi0 (= p pi- pi0 pi0 )
   {0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0007, 0.0013, 0.0029,  0.006,
    0.011, 0.008, 0.006, 0.0053, 0.0037, 0.0027, 0.0017, 0.0008, 0.0003, 0.0001},
  //
  // multiplicity 5 (5 channels)
  //
  // p pi+ pi+ pi- pi- ( = n pi+ pi+ pi- pi- )
   {0.0,    0.0,   0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,   0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0004, 0.0012,
    0.0032, 0.005, 0.0055, 0.0052, 0.005, 0.0048, 0.0046, 0.0042, 0.0038, 0.0033},

  // p pi+ pi- pi0 pi0 ( = n pi+ pi- pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0005, 0.0016,
    0.0041, 0.0064, 0.0076, 0.0081, 0.0075, 0.008, 0.0086, 0.0084, 0.0076, 0.0066},   

  // p pi0 pi0 pi0 pi0 ( = n pi0 pi0 pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0001,
    0.0005, 0.0005, 0.0006, 0.0007, 0.0007, 0.0007, 0.0008, 0.0007, 0.0006, 0.0005},   

  // n pi+ pi+ pi- pi0 ( = p pi+ pi- pi- pi0 )
   {0.0,   0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,   0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0005, 0.0016,
    0.006, 0.0064, 0.0076, 0.0081, 0.009, 0.0093, 0.0086, 0.0084, 0.0076, 0.0066},    

  // n pi+ pi0 pi0 pi0 ( = p pi- pi0 pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0006,
    0.0013, 0.0022, 0.0026, 0.0027, 0.0025, 0.0027, 0.0028, 0.0028, 0.0026, 0.0022},
  //
  // multiplicity 6 (6 channels)
  //
  // p pi+ pi+ pi- pi- pi0 ( = n pi+ pi+ pi- pi- pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0002,
    0.0014, 0.0038, 0.0062, 0.0083, 0.0090, 0.0086, 0.0082, 0.0075, 0.0068, 0.0059},

  // p pi+ pi- pi0 pi0 pi0 ( = n pi+ pi- pi0 pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0009, 0.0017, 0.0029, 0.005, 0.0045, 0.0049, 0.0051, 0.0050, 0.0046, 0.0040},

  // p pi0 pi0 pi0 pi0 pi0 ( = n pi0 pi0 pi0 pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0002, 0.0002, 0.0002, 0.0003, 0.0003, 0.0002, 0.0002, 0.0002},

  // n pi+ pi+ pi+ pi- pi- ( = p pi+ pi+ pi- pi- pi- )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0001,
    0.0003, 0.0009, 0.0026, 0.0032, 0.0035, 0.0033, 0.0031, 0.0028, 0.0025, 0.0022},

  // n pi+ pi+ pi- pi0 pi0 ( = p pi+ pi- pi- pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0002,
    0.0005, 0.0017, 0.0054, 0.0073, 0.0079, 0.0085, 0.0088, 0.0084, 0.0075, 0.0066},

  // n pi+ pi0 pi0 pi0 pi0 ( = p pi- pi0 pi0 pi0 pi0 )
   {0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0001, 0.0003, 0.0009, 0.0012, 0.0013, 0.0013, 0.0014, 0.0014, 0.0012, 0.0011},
  //
  // multiplicity 7 (7 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- ( = n pi+ pi+ pi+ pi- pi- pi- )
   {0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0004, 0.0006, 0.0009, 0.0009, 0.0009, 0.0009, 0.0008, 0.0008},

  // p pi+ pi+ pi- pi- pi0 pi0 ( = n pi+ pi+ pi- pi- pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,   0.0,   0.0,    0.0,   0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,   0.0,   0.0,    0.0,   0.0,    0.0,
    0.0,  0.0003, 0.0012, 0.0021, 0.004, 0.004, 0.0047, 0.005, 0.0045, 0.004},

  // p pi+ pi- pi0 pi0 pi0 pi0 ( = n pi+ pi- pi0 pi0 pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0004, 0.0007, 0.0013, 0.0014, 0.0016, 0.0017, 0.0015, 0.0013},

  // p pi0 pi0 pi0 pi0 pi0 pi0 ( = n pi0 pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0,    0.0,    0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0,    0.0,    0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0001, 0.0001, 0.0, 0.0},

  // n pi+ pi+ pi- pi0 pi0 pi0 ( = p pi+ pi- pi- pi0 pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0003, 0.0008, 0.0014, 0.002, 0.0023, 0.0025, 0.0027, 0.0025, 0.0025},

  // n pi+ pi+ pi+ pi- pi- pi0 ( = p pi+ pi+ pi- pi- pi- pi0 )
   {0.0,  0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0003, 0.0008, 0.0014, 0.002, 0.0023, 0.0025, 0.0027, 0.0025, 0.0025},

  // n pi+ pi0 pi0 pi0 pi0 pi0 ( = p pi- pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0001, 0.0002, 0.0002, 0.0003, 0.0003, 0.0003, 0.0002, 0.0002},
  //
  // multiplicity 8 (7 channels)
  //
  // p pi+ pi+ pi- pi- pi0 pi0 pi0 ( = n pi+ pi+ pi- pi- pi0 pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0003, 0.0006, 0.0015, 0.0026, 0.0033, 0.0034, 0.0034, 0.0033},

  // p pi+ pi+ pi+ pi- pi- pi- pi0 ( = n pi+ pi+ pi+ pi- pi- pi- pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0002, 0.0005, 0.0011, 0.0017, 0.0023, 0.0022, 0.0022, 0.0022},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 ( = n pi+ pi- pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0001, 0.0002, 0.0004, 0.0005, 0.0007, 0.0007, 0.0007, 0.0007},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 ( = n pi0 pi0 pi0 pi0 pi0 pi0 pi0 )
  // negligible 

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 ( = p pi+ pi+ pi- pi- pi- pi0 pi0 )
   {0.0,  0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0,    0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,  0.0001, 0.0003, 0.0006, 0.002, 0.0026, 0.0033, 0.0034, 0.0034, 0.0033},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 ( = p pi+ pi- pi- pi0 pi0 pi0 pi0 )
   {0.0, 0.0, 0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,    0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0001, 0.0003, 0.001, 0.0013, 0.0017, 0.0017, 0.0017, 0.0016},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- ( = p pi+ pi+ pi+ pi- pi- pi- pi- )
   {0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0001, 0.0002, 0.0004, 0.0006, 0.0005, 0.0005, 0.0005},

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 ( = p pi- pi0 pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0, 0.0,  0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001},
  //
  // multiplicity 9 (7 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- pi0 pi0 ( = n pi+ pi+ pi+ pi- pi- pi- pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0002, 0.0005, 0.0013, 0.0019, 0.0023, 0.0026, 0.0026, 0.0024},

  // p pi+ pi+ pi- pi- pi0 pi0 pi0 pi0 ( = n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0002, 0.0003, 0.0011, 0.0014, 0.0018, 0.0020, 0.0020, 0.0018},

  // p pi+ pi+ pi+ pi+ pi- pi- pi- pi- ( = n pi+ pi+ pi+ pi+ pi- pi- pi- pi- )
   {0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0001, 0.0002, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0 ( = n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0, 0.0,  0.0001, 0.0001, 0.0002, 0.0003, 0.0002, 0.0002, 0.0002},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 ( = n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 )
  // negligible

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 pi0 ( = p pi+ pi+ pi- pi- pi- pi0 pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0002, 0.0005, 0.0013, 0.0018, 0.0023, 0.0026, 0.0026, 0.0024},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- pi0 ( = p pi+ pi+ pi+ pi- pi- pi- pi- pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0001, 0.0002, 0.0007, 0.0009, 0.0011, 0.0013, 0.0013, 0.0012},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0 ( = p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0 )
   {0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0, 0.0,  0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.0009, 0.0009, 0.0008}};

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 ( = p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0 )
  // negligible

}

namespace
{
  static double upK2(int i)
  {
    double denom = totalXS[i] - upK2f*sigmaK2[i] - upK3f*sigmaK3[i];
    double num = totalXS[i]  - sigmaK2[i] - sigmaK3[i];
    return upK2f * num/denom;
  }

  static double upK3(int i)
  {
    double denom = totalXS[i] - upK2f*sigmaK2[i] - upK3f*sigmaK3[i];
    double num = totalXS[i]  - sigmaK2[i] - sigmaK3[i];
    return  upK3f * num/denom;
  }
  
  static double biasedTotal(int i)
  {
    double num = totalXS[i]  - sigmaK2[i] - sigmaK3[i];
    return num + upK2(i)*sigmaK2[i] + upK3(i)* sigmaK3[i];
  }


}


// Initialize both |T Tz> = |1/2 +-1/2> channels using gamma-p cross section
// tables
const G4CascadeGamPChannelData::data_t
G4CascadeGamPChannelData::data(gamp2bfs, gamp3bfs, gamp4bfs, gamp5bfs,
                               gamp6bfs, gamp7bfs, gamp8bfs, gamp9bfs,
                               gampCrossSections, gampTotXSec, gam*pro,
                               "GammaProton");

const G4CascadeGamNChannelData::data_t
G4CascadeGamNChannelData::data(gamn2bfs, gamn3bfs, gamn4bfs, gamn5bfs,
                               gamn6bfs, gamn7bfs, gamn8bfs, gamn9bfs,
                               gampCrossSections, gampTotXSec, gam*neu,
                               "GammaNeutron");
