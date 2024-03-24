project(SCEPCALsimG4Full)

file(GLOB libsrc
  ${PROJECT_SOURCE_DIR}/src/lib/*.cpp
)

file(GLOB libheaders
  ${PROJECT_SOURCE_DIR}/include/*.h
)

gaudi_add_library(SCEPCALsimG4Full
  SOURCES ${libsrc}
  LINK
  edm4dr
  ${Geant4_LIBRARIES}
  SCEPCALsegmentation
)

set_target_properties(SCEPCALsimG4Full PROPERTIES PUBLIC_HEADER "${libheaders}")

file(GLOB compsrc
  ${PROJECT_SOURCE_DIR}/src/components/*.cpp
)

file(GLOB compheaders
  ${PROJECT_SOURCE_DIR}/src/components/*.h
)

gaudi_add_module(SCEPCALsimG4FullPlugins
  SOURCES ${compsrc}
  LINK
  Gaudi::GaudiKernel
  k4FWCore::k4FWCore
  ${Geant4_LIBRARIES}
  SCEPCALsegmentation
  SCEPCALsimG4Full
)

set_target_properties(SCEPCALsimG4FullPlugins PROPERTIES PUBLIC_HEADER "${compheaders}")

install(TARGETS SCEPCALsimG4Full EXPORT SimTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

install(TARGETS SCEPCALsimG4FullPlugins EXPORT SimTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
