cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(BetheHeitler)


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

# LightOutput is needed by this simulation
#if( DEFINED ${LightOutput_dir} )
#	message( "LightOutput_dir is defined" )
#else( DEFINED ${LightOutput_dir} )
#	set( LightOutput_dir "${PROJECT_SOURCE_DIR}/../LightOutput-1.1" )
#endif( DEFINED ${LightOutput_dir} )

# Include directories need for this compile
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories("${PROJECT_SOURCE_DIR}/Fakelucid/include")
#include_directories("${LightOutput_dir}/include")

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#message("LightOuptut_dir = ${LightOutput_dir}" )

# get the extra libraries needed
#find_library(G4FakeLucid_lib G4FakeLucid ${PROJECT_SOURCE_DIR}/Fakelucid )
#find_library(LightOutput_lib LightOutput "${LightOutput_dir}" )

#message(${G4FakeLucid_lib} )
#message(${LightOutput_lib} )
# Add the executable, and link it to the Geant4 libraries
#add_executable(fluxnolucid fluxnolucid.cc ${sources} ${headers}) 
add_executable(BetheHeitler BetheHeitler.cc ${sources}) 
target_link_libraries(BetheHeitler ${Geant4_LIBRARIES} )
#target_link_libraries(fluxnolucid ${Geant4_LIBRARIES} ${G4FakeLucid_lib}
#${LightOutput_lib})
 
set(TOCOPY
  gui.mac
  vis.mac
  )

foreach(_script ${TOCOPY})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#add_custom_target(fluxNL depends fluxnolucid)
