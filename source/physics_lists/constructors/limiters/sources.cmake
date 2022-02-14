# - G4phys_ctor_limiters module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_limiters
  PUBLIC_HEADERS
    G4ImportanceBiasing.hh
    G4WeightWindowBiasing.hh
    G4GenericBiasingPhysics.hh
    G4FastSimulationPhysics.hh
    G4MaxTimeCuts.hh
    G4MinEkineCuts.hh
    G4NeutronTrackingCut.hh
    G4SpecialCuts.hh
    G4StepLimiterPhysics.hh
    G4ParallelWorldPhysics.hh
  SOURCES
    G4ImportanceBiasing.cc
    G4WeightWindowBiasing.cc
    G4GenericBiasingPhysics.cc
    G4FastSimulationPhysics.cc
    G4MaxTimeCuts.cc
    G4MinEkineCuts.cc
    G4NeutronTrackingCut.cc
    G4SpecialCuts.cc
    G4StepLimiterPhysics.cc
    G4ParallelWorldPhysics.cc)

geant4_module_link_libraries(G4phys_ctor_limiters
  PUBLIC
    G4biasing_imp
    G4geombias
    G4globman
    G4procman
    G4run
  PRIVATE
    G4baryons
    G4biasing_gen
    G4emutils
    G4hadronic_mgt
    G4navigation
    G4parameterisation
    G4partman
    G4phys_builders
    G4phys_ctor_factory
    G4scoring
    G4track
    G4transportation)
