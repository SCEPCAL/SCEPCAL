project(ddSCEPCAL)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

dd4hep_set_compiler_flags()
dd4hep_use_python_executable()

dd4hep_add_plugin(ddSCEPCAL SOURCES src/*.cpp USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
  SCEPCALsegmentation
)

target_include_directories(ddSCEPCAL PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(ddSCEPCAL PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS ddSCEPCAL
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

install(DIRECTORY compact DESTINATION ${CMAKE_INSTALL_DATADIR})

dd4hep_configure_scripts( ddSCEPCAL DEFAULT_SETUP WITH_TESTS )
