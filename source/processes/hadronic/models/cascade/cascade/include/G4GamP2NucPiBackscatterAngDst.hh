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
// Description: numerically integrated CM angular distribution of the
//              nucleon in the gamma N -> N pi backscatter (u-channel)
//              reactions, from the Anderson et al. (1969) data

#ifndef G4GAMP2NUCPIBACKSCATTERANGDST_HH
#define G4GAMP2NUCPIBACKSCATTERANGDST_HH

#include "G4NumIntTwoBodyAngDst.hh"

class G4GamP2NucPiBackscatterAngDst : public G4NumIntTwoBodyAngDst<7, 14>
{
  public:
    G4GamP2NucPiBackscatterAngDst(G4int verbose = 0);
    ~G4GamP2NucPiBackscatterAngDst() override = default;
};

#endif
