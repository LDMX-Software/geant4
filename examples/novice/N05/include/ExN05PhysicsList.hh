// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PhysicsList.hh,v 1.3 1999/04/23 04:43:12 kurasige Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// 
#ifndef ExN05PhysicsList_h
#define ExN05PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExN05PhysicsList: public G4VUserPhysicsList
{
public:
  ExN05PhysicsList();
  virtual ~ExN05PhysicsList();
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  // 
  virtual void SetCuts();
  
protected:
  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  virtual void ConstructMesons();
  virtual void ConstructBarions();
  virtual void ConstructIons();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();

  virtual void ConstructGeneral();
  virtual void ConstructEM();

  virtual void AddTransportation();
};

#endif


