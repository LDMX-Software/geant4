// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01AnalysisManager.hh,v 1.3 2000/10/31 13:10:00 barrand Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// 

// Example Analysis Manager implementing virtual function
//   RegisterAnalysisSystems.  Exploits C-pre-processor variables
//   G4ANALYSIS_USE_JAS, etc., which are set by the GNUmakefiles if
//   environment variables of the same name are set.

// So all you have to do is set environment variables and compile and
//   instantiate this in your main().

#ifndef AnaEx01AnalysisManager_h
#define AnaEx01AnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "G4AnalysisManager.hh"

class G4Run;
class G4Event;
class G4Step;

class IHistogram1D;
//class ITuple;

class AnaEx01AnalysisManager: public G4AnalysisManager {
public:
  AnaEx01AnalysisManager(const G4String&);
public:
  virtual void BeginOfRun(const G4Run*); 
  virtual void EndOfRun(const G4Run*); 
  virtual void BeginOfEvent(const G4Event*); 
  virtual void EndOfEvent(const G4Event*); 
  virtual void Step(const G4Step*);
private:
  G4int fCalorimeterCollID;                
  IHistogram1D* fEAbs;
  IHistogram1D* fLAbs;
  IHistogram1D* fEGap;
  IHistogram1D* fLGap;
  //ITuple* fTuple;
};

#endif

#endif