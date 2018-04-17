cmake_minimum_required( VERSION 2.8.6 )

project( labelorientcalcs )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( SRCFILES ../itkCentralProcessing/vectorfileread.cxx ../itkCentralProcessing/itkhelpers.cxx ../itkCentralProcessing/itkprocessing.cxx ../itkCentralProcessing/itklineprocessing.cxx ../itkCentralProcessing/itkanglehists.cxx ../itkCentralProcessing/itkmathops.cxx)
#add_definitions(-DmyMESSAGE='"Hello World!"')
#add_definitions(-DmyMESSAGE)
#remove_definitions(-DmyMESSAGE)

#remove_definitions(-DWITH_LINEPROFILES)
#remove_definitions(-DFLOATIMAGET)
#remove_definitions(-DWITH_AERS)
#set(EXECNAME1 "labelorientcalcs")
#add_executable( ${EXECNAME1}  Code.cxx ${SRCFILES} )
#target_link_libraries( ${EXECNAME1} ${ITK_LIBRARIES} )

add_definitions(-DFLOATIMAGET)
add_definitions(-DWITH_LINEPROFILES)
set(EXECNAME2 "labelorientcalcs_withlineprofiles_floatt")
add_executable( ${EXECNAME2}  Code.cxx ${SRCFILES} )
target_link_libraries( ${EXECNAME2} ${ITK_LIBRARIES} )


#remove_definitions(-DWITH_LINEPROFILES)
#add_definitions(-DFLOATIMAGET)
#add_definitions(-DWITH_AERS)
#set(EXECNAME3 "labelorientcalcs_withaers_floatt")
#add_executable( ${EXECNAME3}  Code.cxx ${SRCFILES} )
#target_link_libraries( ${EXECNAME3} ${ITK_LIBRARIES} )

