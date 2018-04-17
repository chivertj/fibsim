cmake_minimum_required( VERSION 2.8.6 )

project( simulate3Dfibres )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

add_executable( simulate3Dfibres Code.cxx ../itkCentralProcessing/vectorfileread.cxx ../itkCentralProcessing/itkhelpers.cxx ../itkCentralProcessing/itkprocessing.cxx ../itkCentralProcessing/itklineprocessing.cxx ../itkCentralProcessing/itkfibresimulator.cxx ../itkCentralProcessing/itkanglehists.cxx ../itkCentralProcessing/itkmathops.cxx)
target_link_libraries( simulate3Dfibres ${ITK_LIBRARIES} )

