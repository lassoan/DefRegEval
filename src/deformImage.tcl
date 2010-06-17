#######################################
# ProstateRegEval uses LPS coordinate system (according to ITK and DICOM standard).
# If metaimages are opened in Slicer they are correctly converted to LPS.
# To properly display STL files in Slicer an LPStoRAS transform shall be applied
# on the model (-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1).
# LPS: ITK, DICOM, meta files, ITKSnap ITK
# RAS: Slicer, NIFTI, ITKSnap NIFTI
# AnatomicalOrientation tag (such as RAI) in meta files are ignored (they are
# assumed to be LPS)
#######################################

proc CreateVolumeMesh { inputFn } {
  log "CreateVolumeMesh"
  set startTime "[clock seconds]"
  set result [ catch {exec "$::tetGenDir/tetgen.exe" -NEFpq1.6a100.0gA $inputFn } ]
  if {$result} {
    log " Warning: CreateVolumeMesh returned with error (usually it is not a problem)"
  } 
  log " [expr [clock seconds]-$startTime]s"
}

proc ExtractSurface { inputImageFn outputSurfaceFn } {
  log "ExtractSurface"
  set startTime "[clock seconds]"
    
  exec "$::binDir/ExtractSurface.exe" \
    --ProstateImageInputFn=$inputImageFn \
    --CombinedSurfMeshOutputFn=$outputSurfaceFn \
    --SupportPosition -5 -85 160 \
    --SupportRadius=120 \
    --ProstateImageThreshold=0.5
    
  log " [expr [clock seconds]-$startTime]s"
}

proc CreateFemModel { objectVolumeMeshFn outputModelFn fx fy fz} {
  log "CreateFemModel" 
  log " f = $fx $fy $fz"
  set startTime "[clock seconds]"
  exec "$::binDir/CreateFemModel" \
    --CombinedVolMeshInputFn=$objectVolumeMeshFn \
    --FemModelOutputFn=$outputModelFn \
    --SupportFixedPosition -40 -10 160 \
    --SupportFixedRadius=65 \
    --ProbePosition1 -12 -120 80 \
    --ProbePosition2 -12 -120 200 \
    --ProbeRadius=10 \
    --ProbeForce $fx $fy $fz \
    --SolverTimeSteps=$::solverTimeSteps
  log " [expr [clock seconds]-$startTime]s"

#    --ProbePosition1 -5 -110 145 
#    --ProbePosition2 -5 -110 190 

}

proc SolveFemModel { inputFn } {
  log "SolveFemModel"
  set startTime "[clock seconds]"
  exec "$::feBioDir/FEBio.exe" -i $inputFn
  log " [expr [clock seconds]-$startTime]s"
} 

proc GetDeformedImage { inputModelFn referenceImageFn outputDeformedImageFn outputDeformedMeshFn outputDeformationFieldImageFn outputFullDeformedImageFn} {
  log "GetDeformedImage"
  set startTime "[clock seconds]"
  exec "$::binDir/GetDeformedImage.exe" --inputModel=$inputModelFn --referenceImage=$referenceImageFn \
    --outputImage=$outputDeformedImageFn \
    --outputFullImage=$outputFullDeformedImageFn \
    --outputDeformationFieldImage=$outputDeformationFieldImageFn \
    --outputVolumeMesh=$outputDeformedMeshFn
  log " [expr [clock seconds]-$startTime]s"
}

proc RegisterImages { fixedImageFn movingImageFn outputRegisteredImageFn outputDeformationFieldImageFn } {
  log "RegisterImages"
  set startTime "[clock seconds]"
  exec "$::binDir/RegisterImages.exe" \
    --fixedImage=$fixedImageFn --movingImage=$movingImageFn \
    --outputRegisteredImage=$outputRegisteredImageFn \
    --outputDeformationField=$outputDeformationFieldImageFn \
    --numberOfIterations=20
  log " [expr [clock seconds]-$startTime]s"
}  

proc EvaluateRegistrationError { femDefFieldFn regDefFieldFn maskImageFn diffDefFieldFn diffMagDefFieldFn } {
  log "EvaluateRegistrationError"
  set startTime "[clock seconds]"
  exec "$::binDir/EvaluateRegistrationError.exe" \
  --deformationReference=$femDefFieldFn --deformationComputed=$regDefFieldFn \
  --deformationDifferenceVector=$diffDefFieldFn --deformationDifferenceMagnitude=$diffMagDefFieldFn \
  --mask=$maskImageFn
  log " [expr [clock seconds]-$startTime]s"
} 
 

#######################################



set solverTimeSteps 5

proc log { txt } {
  puts "$txt"
  set logFile [open $::logFilename a]
  puts $logFile "$txt"
  close $logFile
}

# Source: http://code.activestate.com/recipes/143083/ 
# {{{ Recipe 143083 (r1): Normally Distributed Random Numbers 
# The following procedure generates two independent normally distributed 
# random numbers with mean 0 and vaviance stdDev^2.
proc randNormal {mean stdDev} {
    global dRandNormal1
    global dRandNormal2
    set u1 rand()
    set u2 rand()
    return [expr $mean + $stdDev * sqrt(-2 * log($u1)) * cos(2 * 3.14159 * $u2)]
    # the second number is not used
    # set dRandNormal2 [expr $stdDev * sqrt(-2 * log($u1)) * sin(2 * 3.14159 * $u2)]
}

proc randomBetween { minVal maxVal } {
  return [expr $minVal+($maxVal-$minVal)*rand() ]
}

#######################################

proc start { workDir { perturbForce 0 } { perturbSegmentation 0 } } {

  set ::logFilename "$workDir/results.txt"

  set rawShapeImageFn "$workDir/rawShapeImage.mha"
  set preOpImageFn "$workDir/PreOpImage.mha"
  set preOpGrayscaleImageFn "$workDir/PreOpGrayscaleImage.mha"
  set preOpSegmentedImageFn "$workDir/PreOpSegmentedImage.mha"
  set preOpSurfaceMeshFn "$workDir/PreOpMesh.smesh"
  set preOpVolumeMeshFn "$workDir/PreOpMesh.1.mesh"
  set femModelFn "$workDir/FemDeform.feb"
  set femResultFn "$workDir/FemDeform.plt"
  set intraOpImageFn "$workDir/IntraOpImage.mha"
  set intraOpGrayscaleImageFn "$workDir/IntraOpGrayscaleImage.mha"
  set intraOpMeshFn "$workDir/DeformFieldVolumeMesh.vtu"
  set intraOpSegmentedImageFn "$workDir/IntraOpSegmentedImage.mha"
  set femDefFieldFn "$workDir/DeformFieldFemImage.mha"
  set regDefFieldFn "$workDir/DeformFieldRegImage.mha"
  set intraOpRegisteredImageFn "$workDir/IntraOpSegmentedRegisteredImage.mha"
  set diffDefFieldFn "$workDir/DeformFieldDiffImage.mha"
  set diffMagDefFieldFn "$workDir/DeformFieldDiffMagImage.mha"

  log "Start"

  # Extract and smooth contour
  ExtractSurface $rawShapeImageFn $preOpSurfaceMeshFn

  # Create a tetrahedron mesh in abaqus format
  CreateVolumeMesh $preOpSurfaceMeshFn

  # TODO: Randomize force params
  
  if {$perturbForce} {
    set fx [randNormal 0 5]
    set fy [randNormal -30 20]
    set fz [randNormal 0 2]
  } else {
    set fx -40
    set fy 180
    set fz 2
  }  
  CreateFemModel $preOpVolumeMeshFn $femModelFn $fx $fy $fz
  
  # $::femResultFn is generated from $::femModelFn
  SolveFemModel $femModelFn
  
  GetDeformedImage $femResultFn $preOpGrayscaleImageFn $intraOpImageFn $intraOpMeshFn $femDefFieldFn $intraOpGrayscaleImageFn
  
  log "Done. Exit now."
  
  return 
    
  RegisterImages $preOpGrayscaleImageFn $intraOpGrayscaleImageFn $intraOpRegisteredImageFn $regDefFieldFn

  EvaluateRegistrationError $femDefFieldFn $regDefFieldFn $preOpImageFn $diffDefFieldFn $diffMagDefFieldFn 

  log "End"
}

###########################################################

#set workDir "c:/Users/andras/Documents/ProstateRegEval/test"
set workDir "c:/devel/ProstateTracker/trunk/DeformableReg/debug"
#set binDir "c:/Users/andras/devel/ProstateRegEval-bin/Release"
set binDir "c:/devel/ProstateRegEval-bin/release" 

set shapeModelDir "$workDir/ShapeModel" 
set netGenDir "$binDir/ng431_rel"
set tetGenDir "$binDir/tetgen"
set feBioDir "$binDir/FEBio1p2"

start $workDir 
