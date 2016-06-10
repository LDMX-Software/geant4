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
// G4CascadeCoalescence:  Factory model for final-state interactions to
//   produce light ions from cascade nucleons.  The algorithm implemented
//   here is descirbed in Section 2.3 of the LAQGSM documentation (p. 11-12)
//   [http://lib-www.lanl.gov/la-pubs/00818645.pdf].
//
// 20110917  Michael Kelsey
// 20110920  M. Kelsey -- Use environment variables to set momentum cuts for tuning,
//	     replace polymorphic argument lists with use of "ClusterCandidate"

#ifndef G4CASCADE_COALESCENCE_HH
#define G4CASCADE_COALESCENCE_HH

#include "globals.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzVector.hh"
#include <vector>
#include <set>

class G4CollisionOutput;
class G4InuclElementaryParticle;


class G4CascadeCoalescence {
public:
  G4CascadeCoalescence(G4int verbose=0);
  virtual ~G4CascadeCoalescence();

  // Final state particle list is modified directly
  void FindClusters(G4CollisionOutput& finalState);

  void setVerboseLevel(G4int verbose) { verboseLevel = verbose; }

private:
  typedef std::vector<size_t> ClusterCandidate;	// Indices of constituents

  G4int verboseLevel;				// Control diagnostic messages

  static const G4double dpMaxDoublet;		// Relative momenta for clusters
  static const G4double dpMaxTriplet;
  static const G4double dpMaxAlpha;

  std::vector<ClusterCandidate> allClusters;	// List of candidates found
  std::set<size_t> triedClusters;		// Hashes of combinatorics
  std::set<size_t> usedNucleons;		// List of converted nucleons

  G4CollisionOutput* thisFinalState;		// Pointers to current event
  const std::vector<G4InuclElementaryParticle>* thisHadrons;

  ClusterCandidate thisCluster;			// Reusable buffer for attempts
  G4InuclNuclei thisLightIon;			// Reusable construction buffer

  // Processing stages -- search, construct, cleanup
  void selectCandidates();
  void createNuclei();
  void removeNucleons();

  // Do combinatorics of given nucleons to make candidates
  void tryClusters(size_t idx1, size_t idx2);
  void tryClusters(size_t idx1, size_t idx2, size_t idx3);
  void tryClusters(size_t idx1, size_t idx2, size_t idx3, size_t idx4);

  // Create cluster candidate with listed indices
  void fillCluster(size_t idx1, size_t idx2);
  void fillCluster(size_t idx1, size_t idx2, size_t idx3);
  void fillCluster(size_t idx1, size_t idx2, size_t idx3, size_t idx4);

  // Convert cluster to hash index (for combinatoric reduction)
  size_t clusterHash(const ClusterCandidate& clus) const;

  // Check if candidate cluster has already been evaluated
  bool clusterTried(const ClusterCandidate& clus) const {
    return triedClusters.find(clusterHash(clus)) != triedClusters.end();
  }

  // Check if indexed nucleon is already in a cluster
  bool nucleonUsed(size_t idx) const {
    return usedNucleons.find(idx) != usedNucleons.end();
  }

  // Evaluate conditions for cluster to form light ion
  bool allNucleons(const ClusterCandidate& clus) const;
  bool goodCluster(const ClusterCandidate& clus) const;
  G4int clusterType(const ClusterCandidate& aCluster) const;

  // Extract hadron from final state list
  const G4InuclElementaryParticle& getHadron(size_t idx) const {
    return (*thisHadrons)[idx];
  }

  // Convert candidate nucleon set into output nucleus (true == success)
  bool makeLightIon(const ClusterCandidate& aCluster);

  // Kinematics for cluster evaluations
  G4LorentzVector getClusterMomentum(const ClusterCandidate& aCluster) const;
  mutable G4LorentzVector pCluster;	// Reusable buffer to reduce churn

  G4double maxDeltaP(const ClusterCandidate& aCluster) const;

  // Report cluster arguments for validation
  void reportArgs(const G4String& name, const ClusterCandidate& clus) const;
  void reportResult(const G4String& name, const G4InuclNuclei& nucl) const;
};

#endif	/* G4CASCADE_COALESCENCE_HH */