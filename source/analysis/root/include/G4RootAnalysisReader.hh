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

// The main manager for Root analysis reader.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4RootAnalysisReader_h
#define G4RootAnalysisReader_h 1

#include "G4ToolsAnalysisReader.hh"
#include "globals.hh"

#include "tools/histo/h1d" 
#include "tools/histo/h2d" 
#include "tools/histo/h3d" 
#include "tools/histo/p1d" 
#include "tools/histo/p2d" 
#include "tools/rroot/ntuple"

class G4RootRFileManager;
class G4H1ToolsManager;
class G4H2ToolsManager;
class G4H3ToolsManager;
class G4P1ToolsManager;
class G4P2ToolsManager;
class G4RootRNtupleManager;
  
class G4RootAnalysisReader : public G4ToolsAnalysisReader
{
  public:
    explicit G4RootAnalysisReader(G4bool isMaster = true);
    virtual ~G4RootAnalysisReader();
    
    // static methods
    static G4RootAnalysisReader* Instance();
    
    // Access methods
    tools::rroot::ntuple* GetNtuple() const;
    tools::rroot::ntuple* GetNtuple(G4int ntupleId) const;
    using G4VAnalysisReader::GetNtuple;
    
  protected:
    // virtual methods from base class
    virtual G4int  ReadH1Impl(const G4String& h1Name,  const G4String& fileName,
                              G4bool isUserFileName) final;
    virtual G4int  ReadH2Impl(const G4String& h2Name,  const G4String& fileName,
                              G4bool isUserFileName) final;
    virtual G4int  ReadH3Impl(const G4String& h3Name,  const G4String& fileName,
                              G4bool isUserFileName) final;
    virtual G4int  ReadP1Impl(const G4String& p1Name,  const G4String& fileName,
                              G4bool isUserFileName) final;
    virtual G4int  ReadP2Impl(const G4String& p2Name,  const G4String& fileName,
                              G4bool isUserFileName) final;
    virtual G4int  ReadNtupleImpl(const G4String& ntupleName,  const G4String& fileName,
                              G4bool isUserFileName) final;
 
  private:
    // static data members
    static G4RootAnalysisReader* fgMasterInstance;
    static G4ThreadLocal G4RootAnalysisReader* fgInstance;    

    // methods
    //    
    tools::rroot::buffer* GetBuffer(
                             const G4String& fileName,
                             const G4String& name,
                             const G4String& inFunction);
    G4bool Reset();

    // data members
    G4RootRNtupleManager* fNtupleManager;
    G4RootRFileManager*   fFileManager;
};

#include "G4RootAnalysisReader.icc"

#endif

