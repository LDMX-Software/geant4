#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4expat
# Package: Geant4.src.G4externals.G4expat
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Created on: 02/06/2011
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

# - Ensure clients of the static expat library use the XML_STATIC
# preprocessor def
if(TARGET G4expat-static)
  target_compile_definitions(G4expat-static PUBLIC XML_STATIC)
endif()

