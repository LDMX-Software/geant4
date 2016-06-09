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
// The code was written by :
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef CML2InputDataH
#define CML2InputDataH

#include "G4UImessenger.hh"
#include "ML2MainMessenger.hh"

class CML2MainMessenger;

class CML2CInputData
{
public:
	CML2CInputData(void);
	~CML2CInputData(void);

	inline void setbOnlyVisio(G4bool val){this->bOnlyVisio=val;};
	inline void setPhaseSpaceCentre(G4ThreeVector val){this->inputData.generalData.centrePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setPhaseSpaceHalfSize(G4ThreeVector val){this->inputData.generalData.halfSizePhaseSpace.set(val.getX(), val.getY(), val.getZ());};
	inline void setbSavePhaseSPace(G4bool val){this->inputData.generalData.bSavePhaseSpace=val;};
	inline void setbStopAtPhaseSpace(G4bool val){this->inputData.generalData.bStopAtPhaseSpace=val;};
	inline void setPhaseSpaceOutFile(G4String val){this->inputData.generalData.PhaseSpaceOutFile=val;};

	inline void setbSaveROG(G4bool val){this->inputData.generalData.bSaveROG=val;};
	inline void setROGOutFile(G4String val){this->inputData.generalData.ROGOutFile=val;};

	inline void setMinNumberOfEvents(G4double val){this->inputData.generalData.minNumberOfEvents=val;};

	inline void setBCompareExp(G4bool val){this->inputData.generalData.bCompareExp=val;};
	inline void setFileExperimentalData(G4String val){this->inputData.generalData.fileExperimentalData=val;};
	inline void setNBeams(G4int val){this->inputData.generalData.nBeam=val;};
	inline void setNMaxParticlesInRamPlanePhaseSpace(G4int val){this->inputData.generalData.nMaxParticlesInRamPlanePhaseSpace=val;};

	inline void setSaving_in_Selected_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_Selected_Voxels_every_events=val;};
	inline void setSaving_in_ROG_Voxels_every_events(G4int val){this->inputData.generalData.saving_in_ROG_Voxels_every_events=val;};
	inline void setMax_N_particles_in_PhSp_File(G4int val){this->inputData.generalData.max_N_particles_in_PhSp_File=val;};

	G4bool bOnlyVisio;
	SInputData inputData;
private:
	CML2MainMessenger *ML2MainMessenger;
};


#endif
