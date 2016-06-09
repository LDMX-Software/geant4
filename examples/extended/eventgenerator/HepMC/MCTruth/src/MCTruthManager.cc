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
// $Id: MCTruthManager.cc,v 1.2 2006/12/13 15:42:36 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthManager class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
//
// --------------------------------------------------------------

#include "MCTruthManager.hh"

static MCTruthManager* instance = 0;

MCTruthManager::MCTruthManager() : config(0)
{}

MCTruthManager::~MCTruthManager() 
{} 

MCTruthManager* MCTruthManager::GetInstance()
{
  if( instance == 0 )
  {
    instance = new MCTruthManager();
  }
  return instance;
}

void MCTruthManager::NewEvent()
{
  // first delete the old event
  delete event;
  // and now instaciate a new one
  event = new HepMC::GenEvent();
}

void MCTruthManager::AddParticle(G4LorentzVector& momentum,
                                 G4LorentzVector& prodpos, 
                                 G4LorentzVector& endpos,
                                 G4int pdg_id, G4int partID, G4int motherID,
                                 G4bool directParent)
{
  // we create a new particle with barcode = partID
  HepMC::GenParticle* particle = new HepMC::GenParticle(momentum, pdg_id);
  particle->suggest_barcode(partID);
  // we initialize the 'segmentations' map
  // for the moment particle is not 'segmented' 
  segmentations[partID] = 1;

  // we create the GenVertex corresponding to the end point of the track
  HepMC::GenVertex* endvertex = new HepMC::GenVertex(endpos);
  
  // barcode of the endvertex = - barcode of the track
  endvertex->suggest_barcode(-partID);
  endvertex->add_particle_in(particle);
  event->add_vertex(endvertex);
  
  if(motherID) // not a primary
  {
    // here we could try to improve speed by searching only through particles which 
    // belong to the given primary tree
    HepMC::GenParticle* mother = event->barcode_to_particle(motherID);
    //
    if(mother)
    {
      // we first check whether the mother's end vertex corresponds to the particle's
      // production vertex
      HepMC::GenVertex* motherendvtx = mother->end_vertex();
      G4LorentzVector motherendpos = motherendvtx->position();
      
      if( motherendpos.x() == prodpos.x() &&
          motherendpos.y() == prodpos.y() &&
          motherendpos.z() == prodpos.z() ) // if yes, we attach the particle
      {
        motherendvtx->add_particle_out(particle);
      }
      else // if not, we check whether the mother is biological or adopted
      {            
        if(!directParent) // adopted
        {  
          G4bool found = false;

          // first check if any of the dummy particles
          // has the end vertex at the right place
          //
          for(HepMC::GenVertex::particles_out_const_iterator 
                it=motherendvtx->particles_out_const_begin();
              it!=motherendvtx->particles_out_const_end();it++)
          {
            if((*it)->pdg_id()==-999999)
            {
              G4LorentzVector dummypos = (*it)->end_vertex()->position();
              
              if( dummypos.x() == prodpos.x() &&
                  dummypos.y() == prodpos.y() &&
                  dummypos.z() == prodpos.z() ) 
              {
                (*it)->end_vertex()->add_particle_out(particle);
                found = true;
                break;
              }
            }
          }

          // and if not, create a dummy particle connecting
          // to the end vertex of the mother
          //
          if(!found)
          {
            HepMC::GenVertex* childvtx = new HepMC::GenVertex(prodpos);
            childvtx->add_particle_out(particle);

            // the dummy vertex gets the barcode -500000
            // minus the daughter particle barcode
            //
            childvtx->suggest_barcode(-500000-partID);
            event->add_vertex(childvtx);
            
            HepMC::GenParticle* dummypart =
               new HepMC::GenParticle(G4LorentzVector(),-999999);

            // the dummy particle gets the barcode 500000
            // plus the daughter particle barcode
            //
            dummypart->suggest_barcode(500000+partID);
            childvtx->add_particle_in(dummypart);
            motherendvtx->add_particle_out(dummypart);
          }
        }
        else // biological
        {
          // in case mother was already 'split' we need to look for
          // the right 'segment' to add the new daugther.
          // We use Time coordinate to locate the place for the new vertex

          G4int number_of_segments = segmentations[motherID];
          G4int segment = 0;

          // we loop through the segments
          //         
          while ( !((mother->end_vertex()->position().t()>prodpos.t()) && 
                    (mother->production_vertex()->position().t()<prodpos.t())) )
          {
            segment++;
            if (segment == number_of_segments) 
              G4cerr << "Problem!!!! Time coordinates incompatible!" << G4endl;
            
            mother = event->barcode_to_particle(segment*10000000 + motherID);
          }
          
          // now, we 'split' the appropriate 'segment' of the mother particle
          // into two particles and create a new vertex
          //
          HepMC::GenVertex* childvtx = new HepMC::GenVertex(prodpos);
          childvtx->add_particle_out(particle);
          event->add_vertex(childvtx);

          // we first detach the mother from its original vertex
          //
          HepMC::GenVertex* orig_mother_end_vtx = mother->end_vertex();
          orig_mother_end_vtx->remove_particle(mother);

          // and attach it to the new vertex
          //
          childvtx->add_particle_in(mother);

          // now we create a new particle representing the mother after
          // interaction the barcode of the new particle is 10000000 + the
          // original barcode
          //
          HepMC::GenParticle* mothertwo = new HepMC::GenParticle(*mother);
          mothertwo->suggest_barcode(segmentations[motherID]*10000000
                                    + mother->barcode());

          // we also reset the barcodes of the vertices
          //
          orig_mother_end_vtx->suggest_barcode(-segmentations[motherID]
                                               *10000000 - mother->barcode());
          childvtx->suggest_barcode(-mother->barcode());

          // we attach it to the new vertex where interaction took place
          //
          childvtx->add_particle_out(mothertwo);

          // and we attach it to the original endvertex
          //
          orig_mother_end_vtx->add_particle_in(mothertwo);

          // and finally ... the increase the 'segmentation counter'
          //
          segmentations[motherID] = segmentations[motherID]+1;
        }
      }
    }
    else 
      // mother GenParticle is not there for some reason...
      // if this happens, we need to revise the philosophy... 
      // a solution would be to create HepMC particles
      // at the begining of each track
    {
      G4cerr << "barcode " <<  motherID << " mother not there! "<<  G4endl;
    }
  }
  else // primary
  {
    HepMC::GenVertex* primaryvtx = new HepMC::GenVertex(prodpos);
    primaryvtx->add_particle_out(particle);
    event->add_vertex(primaryvtx);

    // add id to the list of primaries
    //
    primarybarcodes.push_back(partID);
  } 
}

void MCTruthManager::PrintEvent()
{
  event->print();

  // looping over primaries and print the decay tree for each of them
  //
  for(std::vector<int>::const_iterator primarybar=primarybarcodes.begin();
      primarybar!=primarybarcodes.end();primarybar++)
  {
    printTree(event->barcode_to_particle(*primarybar), " | ");
  }
}

void MCTruthManager::printTree(HepMC::GenParticle* particle, G4String offset)
{
  G4cout << offset << "---  barcode: " << particle->barcode() << " pdg: " 
         << particle->pdg_id() << " energy: " << particle->momentum().e() 
         << " production vertex: "<< particle->production_vertex()->position()
         << G4endl;

  for(HepMC::GenVertex::particles_out_const_iterator 
        it=particle->end_vertex()->particles_out_const_begin();
        it!=particle->end_vertex()->particles_out_const_end();
      it++)
  {
    G4String deltaoffset = "";

    if( std::fmod((*it)->barcode(),10000000) != std::fmod(particle->barcode(),10000000) )
      {
        deltaoffset = " | ";
      }

    printTree((*it), offset + deltaoffset);
  } 
}