ADD_LIBRARY(dreIO
	IO/vtkLsDynaBinaryPlotReader.cxx 
	IO/vtkMeditMeshReader.cxx 
	IO/vtkNetgenMeshReader.cxx 
	IO/vtkSMeshWriter.cxx
	IO/vtkFEBioWriter.cxx 
	)

##########################

ADD_EXECUTABLE(CreateShapeModel 
  CreateShapeModel.cxx
  )
  
TARGET_LINK_LIBRARIES(CreateShapeModel 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys
  )

##########################

ADD_EXECUTABLE(GenerateShapeImage 
  GenerateShapeImage.cxx
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
  )
  
TARGET_LINK_LIBRARIES(ExtractSurface
  vtkIO  
  vtkGraphics
  vtksys
  vtkHybrid
  vtkImaging 
  dreIO
  )

##########################

ADD_EXECUTABLE(CreateFemModel 
  CreateFemModel.cxx
  )
  
TARGET_LINK_LIBRARIES(CreateFemModel 
  ITKCommon 
  ITKIO 
  ITKBasicFilters 
  ITKAlgorithms 
  vtksys 
  vtkIO  
  vtkGraphics
  dreIO
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
  dreIO 
  )
  
##########################

ADD_EXECUTABLE(RegisterImages 
  RegisterImages.cxx
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
