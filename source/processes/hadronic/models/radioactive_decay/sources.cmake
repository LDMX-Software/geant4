# - G4hadronic_radioactivedecay module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_radioactivedecay
  PUBLIC_HEADERS
    G4AlphaDecay.hh
    G4BatemanParameters.hh
    G4BetaDecayCorrections.hh
    G4BetaDecayType.hh
    G4BetaMinusDecay.hh
    G4BetaPlusDecay.hh
    G4ECDecay.hh
    G4ITDecay.hh
    G4SFDecay.hh
    G4NeutronDecay.hh
    G4NuclearDecay.hh
    G4NucleusLimits.hh
    G4ProtonDecay.hh
    G4RadioactiveDecay.hh
    G4RadioactiveDecayMessenger.hh
    G4Radioactivation.hh
    G4RadioactivationMessenger.hh
    G4RadioactiveDecayMode.hh
    G4RadioactiveDecayRatesToDaughter.hh
    G4RadioactiveDecayChainsFromParent.hh
    G4RadioactivityTable.hh
    G4TritonDecay.hh
    G4UIcmdWithNucleusLimits.hh
    G4UserLimitsForRD.hh
  SOURCES
    G4AlphaDecay.cc
    G4BatemanParameters.cc
    G4BetaDecayCorrections.cc
    G4BetaDecayType.cc
    G4BetaMinusDecay.cc
    G4BetaPlusDecay.cc
    G4ECDecay.cc
    G4ITDecay.cc
    G4SFDecay.cc
    G4NeutronDecay.cc
    G4NuclearDecay.cc
    G4NucleusLimits.cc
    G4ProtonDecay.cc
    G4RadioactiveDecay.cc
    G4RadioactiveDecayMessenger.cc
    G4Radioactivation.cc
    G4RadioactivationMessenger.cc
    G4RadioactiveDecayMode.cc
    G4RadioactiveDecayRatesToDaughter.cc
    G4RadioactiveDecayChainsFromParent.cc
    G4RadioactivityTable.cc
    G4TritonDecay.cc
    G4UIcmdWithNucleusLimits.cc)

geant4_module_link_libraries(G4hadronic_radioactivedecay
  PUBLIC
    G4globman
    G4hadronic_util
    G4heprandom
    G4intercoms
    G4partman
    G4procman
    G4track
  PRIVATE
    G4baryons
    G4bosons
    G4emlowenergy
    G4emutils
    G4geometrymng
    G4had_fission
    G4hadronic_deex_management
    G4hadronic_deex_photon_evaporation
    G4hadronic_mgt
    G4hadronic_util
    G4hepgeometry
    G4ions
    G4leptons
    G4materials)
