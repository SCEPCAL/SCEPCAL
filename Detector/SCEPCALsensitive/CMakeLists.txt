project(SCEPCALsensitive)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

add_library(SCEPCALsensitive SHARED ${sources} ${headers})

target_include_directories(SCEPCALsensitive PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(SCEPCALsensitive PROPERTIES PUBLIC_HEADER "${headers}")

target_link_libraries(
  SCEPCALsensitive
  SCEPCALsegmentation
  edm4dr
  edm4dr::edm4drDict
  DD4hep::DDCore
  DD4hep::DDG4
)

dd4hep_generate_rootmap(SCEPCALsensitive)

install(TARGETS SCEPCALsensitive EXPORT DetectorTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

#---Geant4 Testing-----------------------------------------------------------------
#

#---------------------------  Plugin library for the simulation framework  ---------
# dd4hep_add_dictionary(SCEPCALsensitive
  # SOURCES ${DD4hep_DIR}/include/ROOT/Warnings.h include/ScepcalHit.h
  # LINKDEF ${DD4hep_DIR}/include/ROOT/LinkDef.h
  # OUTPUT  ${LIBRARY_OUTPUT_PATH}
  # USES    DD4hep::DDCore DD4hep::DDG4 Geant4::Interface
  # )

#----  Example of a client library with user defined plugins  --------------------
# dd4hep_add_plugin( SCEPCALsensitive
#   GENERATED G__SCEPCALsensitive.cxx
#   SOURCES   src/*.cpp
#   USES      DD4hep::DDCore DD4hep::DDG4 Geant4::Interface ROOT::Core ROOT::Geom ROOT::GenVector ROOT::RIO
#   )
# # install(TARGETS SCEPCALsensitive DESTINATION lib)


