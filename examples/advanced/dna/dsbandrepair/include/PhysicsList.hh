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
//
/// \file PhysicsList.hh
/// \brief Definition of the PhysicsList class

#ifndef PhysicsList_h
#define PhysicsList_h 1
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <memory>

#include "G4VPhysicsConstructor.hh"
#include "G4GenericMessenger.hh"
#include "G4EmDNAChemistry_option2.hh"
#include "G4EmDNAChemistry_option3_Extended.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;
    
    void ConstructParticle() override;
    void ConstructProcess() override;

    void RegisterPhysicsList(const G4String& name);
    void RegisterChemListConstructor(const G4String& name);
    void SetChemListName(const G4String& cname) {fChemListName = cname;}
private:
    std::unique_ptr<G4VPhysicsConstructor>    fDNAPhysicsList{nullptr};
    std::unique_ptr<G4VPhysicsConstructor> fEmDNAChemistryList{nullptr};
    std::unique_ptr<G4GenericMessenger> fMessenger;
    G4String                                  fPhysDNAName{""};
    G4String                                  fChemListName{""};
    void DefineCommands();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif