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
