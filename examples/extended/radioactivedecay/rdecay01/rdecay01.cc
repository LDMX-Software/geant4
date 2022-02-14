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
/// \file radioactivedecay/rdecay01/rdecay01.cc
/// \brief Main program of the radioactivedecay/rdecay01 example
//
//
<<<<<<< HEAD
// $Id: rdecay01.cc 94458 2015-11-18 14:36:25Z gcosmo $
// 
=======
//
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

<<<<<<< HEAD
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
=======

int main(int argc,char** argv) {
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {
 
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  //construct the run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();  
  if (argc==3) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
  }  

  // set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);
      
  // set user action classes
  // 
  runManager->SetUserInitialization(new ActionInitialization);
    
  //Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  if (argc!=1)   // batch mode  
    { 
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);  
    }
    
  else           // define visualization and UI terminal for interactive mode 
    { 
#ifdef G4UI_USE
     G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
#ifdef G4VIS_USE
     UI->ApplyCommand("/control/execute vis.mac");          
#endif
     ui->SessionStart();
     delete ui;
#endif
    }

#ifdef G4VIS_USE
  delete visManager;
#endif     

  // job termination
  //
  delete runManager;

  return 0;
}

<<<<<<< HEAD
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
