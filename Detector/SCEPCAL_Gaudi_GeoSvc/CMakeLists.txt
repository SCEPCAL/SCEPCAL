project(SCEPCAL_Gaudi_GeoSvc)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

add_library(SCEPCAL_Gaudi_GeoSvc SHARED ${sources} ${headers})

target_include_directories(SCEPCAL_Gaudi_GeoSvc PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${DD4hep_DIR}/include
)

target_link_libraries(
  SCEPCAL_Gaudi_GeoSvc
  DD4hep::DDCore
  DD4hep::DDG4
)

install(TARGETS SCEPCAL_Gaudi_GeoSvc
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
