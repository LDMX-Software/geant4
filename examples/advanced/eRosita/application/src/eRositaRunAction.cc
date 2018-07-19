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
// $Id: eRositaRunAction.cc 107396 2017-11-10 08:28:08Z gcosmo $
// 

#include "eRositaRunAction.hh"
#include "AnalysisManager.hh"

#include "G4Run.hh"


eRositaRunAction::eRositaRunAction()
{}


eRositaRunAction::~eRositaRunAction()
{}



void eRositaRunAction::BeginOfRunAction(const G4Run* run)
{
  timerRun.Start();	

  G4cout << "--- Run " << run->GetRunID() << " start." << G4endl;

  AnalysisManager::Instance();
}



void eRositaRunAction::EndOfRunAction(const G4Run* run)
{ 

  AnalysisManager::Instance()->Destroy();

  G4cout << "--- Run " << run->GetRunID() 
	 << run -> GetRunID()
         << " ends (Number of events = "
         << run -> GetNumberOfEvent() << ")." 
	 << G4endl;

  timerRun.Stop();
  std::cout << "  "  << timerRun << std::endl;
}





