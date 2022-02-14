# - G4muons module build definition

# Define the Geant4 Module.
geant4_add_module(G4muons
  PUBLIC_HEADERS
    G4EnergyLossForExtrapolator.hh
    G4ErrorEnergyLoss.hh
    G4ModifiedMephi.hh
    G4MuBetheBlochModel.hh
    G4MuBremsstrahlung.hh
    G4MuBremsstrahlungModel.hh
    G4MuIonisation.hh
    G4MuMultipleScattering.hh
    G4MuPairProduction.hh
    G4MuPairProductionModel.hh
    G4TablesForExtrapolator.hh
    G4ePairProduction.hh
  SOURCES
    G4EnergyLossForExtrapolator.cc
    G4ErrorEnergyLoss.cc
    G4ModifiedMephi.cc
    G4MuBetheBlochModel.cc
    G4MuBremsstrahlung.cc
    G4MuBremsstrahlungModel.cc
    G4MuIonisation.cc
    G4MuMultipleScattering.cc
    G4MuPairProduction.cc
    G4MuPairProductionModel.cc
    G4TablesForExtrapolator.cc
    G4ePairProduction.cc)

geant4_module_link_libraries(G4muons
  PUBLIC
    G4emutils
    G4globman
    G4materials
    G4procman
  PRIVATE
    G4baryons
    G4bosons
    G4cuts
    G4emstandard
    G4heprandom
    G4leptons
    G4partman
    G4track)
