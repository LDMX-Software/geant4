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
//
// Author:  Natalia Toro (SLAC)
// Date:    20 June 2018
//
// Description: CM angular distribution of the nucleon in gamma N -> N pi
//              backscatter. Anderson et al. (1969), PRL 23, 721, found a
//              near-universal k^3 dsigma/du for u > -1.8 GeV^2 at 4-18 GeV
//              (pi+ n, similar for p pi0). It is well fit by
//              k^3 dsigma/du = (0.4/((u + 0.07)^2 + 0.15) + 0.35) 1e-3 mb GeV,
//              and the table below integrates that fit.

#include "G4GamP2NucPiBackscatterAngDst.hh"

namespace
{
static const G4double eBins[7] = {2.450, 2.900, 4.0, 6.0, 9.0, 12.0, 18.0};

static const G4double angleBins[14] = {-1.000, 0.000, 0.174, 0.342, 0.500, 0.643, 0.766,
                                       0.866,  0.900, 0.940, 0.960, 0.985, 0.993, 1.000};

static const G4double integralTable[7][14] = {
  {0.0000, 0.0000, 0.0327, 0.1104, 0.2047, 0.3291, 0.5006, 0.7067, 0.7856, 0.8769, 0.9202,
   0.9713, 0.9868, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0405, 0.1346, 0.2522, 0.4142, 0.6325, 0.7265, 0.8413, 0.8972,
   0.9631, 0.9830, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0000, 0.0035, 0.1215, 0.2690, 0.4824, 0.5939, 0.7532, 0.8384,
   0.9421, 0.9735, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1104, 0.3065, 0.4164, 0.6070, 0.7320,
   0.9022, 0.9553, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1516, 0.2540, 0.4423, 0.5911,
   0.8413, 0.9274, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0438, 0.1475, 0.3286, 0.4792,
   0.7812, 0.8991, 1.0000},
  {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1785, 0.3225,
   0.6691, 0.8424, 1.0000}};
}  // namespace

// Constructor passes arrays to templated base class

G4GamP2NucPiBackscatterAngDst::G4GamP2NucPiBackscatterAngDst(G4int verbose)
  : G4NumIntTwoBodyAngDst<7, 14>("G4GamP2NucPiBackscatterAngDst", eBins, angleBins,
                                 integralTable, 1.5, verbose)
{}
