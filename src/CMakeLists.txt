PROJECT(ProstateErrorMapRendering)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

## KWWidgets
find_package(KWWidgets REQUIRED)
if(NOT KWWidgets_FOUND)
  message(FATAL_ERROR "Need KWWidgets to build")
endif(NOT KWWidgets_FOUND)
include(${KWWidgets_USE_FILE})

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ADD_EXECUTABLE(ProstateErrorMapRendering 
  ProstateErrorMapRendering.cxx 
  vtkKWProstateErrorMapRenderingWidget.cxx
  vtkKWProstateErrorMapRenderingWidget.h 
  Random.cxx 
  Random.h 
  ContourStatistics.cxx 
  ContourStatistics.h 
  CubicSplineSolver.cxx 
  CubicSplineSolver.h 
  Polynomial.cxx 
  Polynomial.h 
  ExtraFunctions.cxx 
  ExtraFunctions.h   
  )
  
TARGET_LINK_LIBRARIES(ProstateErrorMapRendering 
  KWWidgets 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKFEM 
  ITKSTATISTICS 
  ITKAlgorithms 
  vtkRendering 
  vtkGraphics 
  vtkHybrid 
  vtkImaging 
  vtkIO 
  vtkFiltering 
  vtkCommon
  )

##########################

ADD_EXECUTABLE(GenerateShapeImage 
  GenerateShapeImage.cxx
  Random.cxx 
  Random.h  
  )
  
TARGET_LINK_LIBRARIES(GenerateShapeImage 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  )

##########################

ADD_EXECUTABLE(ExtractSurface
  ExtractSurface.cxx  
  vtkSMeshWriter.cxx
  vtkSMeshWriter.h
  )
  
TARGET_LINK_LIBRARIES(ExtractSurface
  vtkIO  
  vtkGraphics
  vtksys
  vtkHybrid
  vtkImaging 
  )

##########################

ADD_EXECUTABLE(CreateFemModel 
  CreateFemModel.cxx
  vtkNetgenMeshReader.cxx
  vtkNetgenMeshReader.h
  vtkMeditMeshReader.cxx
  vtkMeditMeshReader.h  
  vtkFEBioWriter.cxx
  vtkFEBioWriter.h
  Random.cxx 
  Random.h  
  )
  
TARGET_LINK_LIBRARIES(CreateFemModel 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  vtkIO  
  vtkGraphics
  )

##########################

ADD_EXECUTABLE(AddImageSegmentationError 
  AddImageSegmentationError.cxx
  )
  
TARGET_LINK_LIBRARIES(AddImageSegmentationError 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  )

##########################

ADD_EXECUTABLE(GetDeformedImage 
  GetDeformedImage.cxx
  vtkLsDynaBinaryPlotReader.cxx
  vtkLsDynaBinaryPlotReader.h
  )
  
TARGET_LINK_LIBRARIES(GetDeformedImage  
  vtksys 
  vtkIO  
  vtkGraphics
  vtkImaging
  vtkHybrid 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  )
  
##########################

ADD_EXECUTABLE(RegisterImages 
  RegisterImages.cxx
  Random.cxx 
  Random.h  
  )
  
TARGET_LINK_LIBRARIES(RegisterImages 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  )
  
##########################

ADD_EXECUTABLE(EvaluateRegistrationError 
  EvaluateRegistrationError.cxx
  )
  
TARGET_LINK_LIBRARIES(EvaluateRegistrationError 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  )
  
##########################

ADD_EXECUTABLE(ComputeMultiImageStatistics 
  ComputeMultiImageStatistics.cxx
  )
  
TARGET_LINK_LIBRARIES(ComputeMultiImageStatistics 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  )
