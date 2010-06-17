#ifndef __vtkKWProstateErrorMapRenderingWidget_h
#define __vtkKWProstateErrorMapRenderingWidget_h

#include "ContourStatistics.h"

//Includes for VTK
#include "vtkImageData.h"
#include "vtkRenderer.h"
#include "vtkObjectFactory.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkVolumeMapper.h"
#include "vtkMarchingCubes.h"
#include "vtkDataSetMapper.h"

//Includes for KWWidgets
#include "vtkKWRenderWidget.h"
#include "vtkKWCompositeWidget.h"
#include "vtkKWApplication.h"
#include "vtkKWFrame.h"
#include "vtkKWFileBrowserDialog.h"
#include "vtkKWWindowBase.h"
#include "vtkKWMenu.h"
#include "vtkKWScale.h"
#include "vtkKWLabel.h"
#include <vtksys/SystemTools.hxx>

//Includes for ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkImagePCAShapeModelEstimator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkLBFGSOptimizer.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkImageRegistrationMethod.h"
#include "itkLevelSetMotionRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include <itksys/SystemTools.hxx>

//Other includes
#include <sys/stat.h> 
#include <windows.h>
#include <fstream>
#include "Random.h"
#include "ExtraFunctions.h"
#include "Polynomial.h"
#include "CubicSplineSolver.h"

class vtkImageViewer;
class vtkKWRenderWidget;

//declarations for typedefs, consts, and structs
const int DEFAULT_WINDOW_WIDTH = 1268;
const int DEFAULT_WINDOW_HEIGHT = 852;
const double PI = 3.1415926535;
const double TO_RADIANS = (double)PI/(double)180;

static const unsigned int  Dimension = 3;
typedef float  PixelType;
typedef float VectorComponentType;
const unsigned int SplineOrder = 3;
typedef double CoordinateRepType;
typedef itk::Image <PixelType, Dimension> InternalImageType;
typedef itk::Image <unsigned char, 3> OutputImageType;
typedef itk::Vector<VectorComponentType, Dimension > VectorPixelType;
typedef itk::Image<VectorPixelType,  Dimension> DeformationFieldType;
typedef itk::ImageToVTKImageFilter<OutputImageType> ITKImageToVTKFilterType;
typedef itk::ImageToVTKImageFilter<DeformationFieldType> ITKDeformationFieldToVTKFilterType;

template <class InputPixelType>
class ContourOfSignedDistanceMapFunctor
{
public:
 
   InputPixelType operator()( InputPixelType input )
       {
		   if( input <= 0)
			    return 1;
		   else
				return 0;
       }
};

struct Point
{
	int x;
	int y;
};

struct TargetPoint //voxel space
{
	int x;
	int y;
	int z;
};

struct PointPolar
{
	float r;
	float theta;
};

//for the targeting error heat map
struct VoxelTargetingErrorStruct
{
	bool targeted;
	double minError;
	double maxError;
	double meanError;
	double standardDeviation;
	std::vector<double> errors;
};

struct ReportStruct
{
	int totalTargets;
	int numInsignificant;
	int numSignificant;
	double percentInsignificant;
	double percentSignificant;
	double minError;
	double maxError;
	double meanError;
	double meanDeviation;
	double standardDeviation;
};

//Main class of the simulator. Handles all user events, generates the mean shape
//and modes of variation, runs the simulation, records targeting errors, and
//runs the display
class vtkKWProstateErrorMapRenderingWidget : public vtkKWCompositeWidget
{
	public:
		struct TargetPointStruct
		{
			std::vector<TargetPoint> targetPointsITK;
			std::vector<TargetPoint> targetPointsVTK;
		};

		static vtkKWProstateErrorMapRenderingWidget* New();
		vtkTypeRevisionMacro(vtkKWProstateErrorMapRenderingWidget, vtkKWCompositeWidget);

		void SetParentWindow(vtkKWWindowBase*);

		PointPolar ConvertToPolarCoordinates(Point, Point);

	protected:
	  vtkKWProstateErrorMapRenderingWidget();
	  ~vtkKWProstateErrorMapRenderingWidget();

	  struct SimulationInfoStruct
	  {
		  double sigmaFactor[3];
	  };

	  //constants

	  //strings in the menus
	  static const char* LOAD_PROSTATE_CONTOUR_MENU_STRING;
	  static const char* GENERATE_MEAN_SHAPE_MENU_STRING;
	  static const char* RUN_SIMULATION_MENU_STRING;
	  static const char* EXIT_MENU_STRING;

	  //GUI constants
	  static const int RENDER_WIDGET_PADDING = 2;
	  static const int MAIN_RENDER_WIDGET_WIDTH = 1266;
	  static const int MAIN_RENDER_WIDGET_HEIGHT = 774;
	  static const int SMALL_RENDER_WIDGET_WIDTH = 315;
	  static const int SMALL_RENDER_WIDGET_HEIGHT = 206;

	  //initial azimuth angle for viewing the prostate contours (on the screen
	  //with multiple contours)
	  #ifdef _DEBUG
		static const int RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE = -30;
	  #else
		static const int RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE = 150;
	  #endif

	  static const int RENDER_WIDGET_INITIAL_AZIMUTH_SINGLE = -30; //initial azimuth angle for viewing the mean shape
	  static const int RENDER_WIDGET_INITIAL_ELEVATION = 15; //initial azimuth angle for viewing the mean shape
	  static const double MIN_WEIGHT_MODE_1; //weights for the modes of variatioin
	  static const double MAX_WEIGHT_MODE_1;
	  static const double MIN_WEIGHT_MODE_2;
	  static const double MAX_WEIGHT_MODE_2;
	  static const double MIN_WEIGHT_MODE_3;
	  static const double MAX_WEIGHT_MODE_3;
	  static const int SEGMENTATION_ERROR_SEGMENTS = 31; //parameter "N" in the paper
	  static const int MIN_FLUCTUATION_POINTS = 5; //parameter "M" in the paper
	  static const double MAX_SEGMENTATION_ERROR_DIFFERENCE; //parameter "D" in the paper
	  static const int DEFORMATION_SEGMENTS = 3; //not used
	  static const double DEFORMATION_PROBABILITY; //probability of deformation "p" in the paper
	  static const int PERCENT_BASE = 25;
	  static const int PERCENT_APEX = 25;
	  static const int PERCENT_MID = 100 - PERCENT_BASE - PERCENT_APEX;
	  static const double MIN_MID_SEGMENTATION_ERROR; //mm     //parameters dmin and dmax in the paper
	  static const double MAX_MID_SEGMENTATION_ERROR; //mm
	  static const double MIN_BASE_SEGMENTATION_ERROR; //mm
	  static const double MAX_BASE_SEGMENTATION_ERROR; //mm
	  static const double MIN_APEX_SEGMENTATION_ERROR; //mm
	  static const double MAX_APEX_SEGMENTATION_ERROR; //mm
	  static const double MIN_DEFORMATION_MAGNITUDE; //mm     //parameter "md" in the paper
	  static const double MAX_DEFORMATION_MAGNITUDE; //mm     //parameter "Md" in the paper
	  static const int PROSTATE_PADDING = 20; //padding placed around the prostate in the ITK image
	  static const bool TARGET_ALL = true; //target surface as well as the interior of the prostate
	  static const bool TARGET_SURFACE = true; //targeting only the surface of the prostate
	  static const bool RANDOM_TARGETING = false; //use randomized targeting. If false, use systematic targeting
	  static const double ANGLE_NOT_TARGETED; //parameter theta in the section "Targeting parameters and targeting error"
	  static const int MAX_TARGET_POINTS_TO_RENDER = 8000; //maximum number of target points to render (target points are spheres so there is a slow down if too many are drawn)
	  static const double RANDOM_TARGET_PERCENT; //percentage of points to target when using random targeting
	  static const double PERIPHERAL_ZONE_PERCENT; //parameter "P" in the section "Targeting parameters and targeting error"
	  static const int BIOPSY_SAMPLING_Z = 3; //for systematic targeting, how many targets there are along the z-axis (inferior-superior)
	  static const int NUM_TARGETS[BIOPSY_SAMPLING_Z]; //number of target points along the left-right axis for each section in the targeting
	  static const double MAX_INSIGNIFICANT_ERROR; //mm   //for the targeting error heat map...anything below "MAX_INSIGNIFICANT_ERROR" is insignificant
	  static const double MIN_SIGNIFICANT_ERROR; //mm     //for the targeting error heat map...anything above "MIN_SIGNIFICANT_ERROR" is significant
	  static const int SIMULATION_CASES = 1; //how many simulation cases to run
	  static const bool USING_DEMONS_REGISTRATION = false;
	  static const bool USING_BSPLINES_REGISTRATION = false;
	  static const bool USING_LEVEL_SETS_REGISTRATION = true;
	  static const bool DO_REGISTRATION = true;
	  static const bool DO_DEFORMATION = true;

	  //instance variables
	  std::string baseSegmentationErrorString;
	  std::string midSegmentationErrorString;
	  std::string apexSegmentationErrorString;
	  bool multipleImageFramesShowing;
	  vtkKWWindowBase* ParentWindow;
	  ITKImageToVTKFilterType::Pointer ITKImageToVTKFilter;
	  InternalImageType::Pointer MeanShape;
	  OutputImageType::Pointer MeanShapeDisplay;
	  std::vector<InternalImageType::Pointer> EigenModes;
	  std::vector<SimulationInfoStruct> SimulationInfo;
	  ContourStatistics MeanShapeStatistics;
	  std::vector<ContourStatistics> GroundTruthContourStatistics;
	  std::vector<ContourStatistics> DeformationContourStatistics;
	  std::vector<ContourStatistics> RegistrationResultsContourStatistics;
	  std::vector<OutputImageType::Pointer> ITKContours; //ground truth pre-op contours
	  std::vector<OutputImageType::Pointer> ITKContoursSegmentationError; //pre-op contours
	  std::vector<OutputImageType::Pointer> ITKContoursSegmentationErrorDifferences;
	  std::vector<OutputImageType::Pointer> ITKContoursWithDeformation;//ground truth intra-op contours
	  std::vector<OutputImageType::Pointer> ITKContoursWithDeformationAndSegmentationError;//intra-op contours
	  std::vector<OutputImageType::Pointer> ITKContoursWithDeformationAndSegmentationErrorDifferences;
	  std::vector<OutputImageType::Pointer> ITKContoursRegistrationResults;//registered contours
	  std::vector<OutputImageType::Pointer> ITKContoursRegistrationDifferences;
	  std::vector<DeformationFieldType::Pointer> GroundTruthDeformationFields;
	  std::vector<DeformationFieldType::Pointer> RegistrationResultDeformationFields;
	  std::vector<std::vector<TargetPoint>> GroundTruthPreOpTargetLocationsITK;
	  std::vector<std::vector<TargetPoint>> GroundTruthPreOpTargetLocationsVTK;
	  std::vector<std::vector<TargetPoint>> GroundTruthIntraOpTargetLocationsITK;
	  std::vector<std::vector<TargetPoint>> GroundTruthIntraOpTargetLocationsVTK;
	  std::vector<std::vector<TargetPoint>> ResultTargetLocationsITK;
	  std::vector<std::vector<TargetPoint>> ResultTargetLocationsVTK;
	  std::vector<ReportStruct> TargetingErrorReports; //targeting error reports for each simulation
	  VoxelTargetingErrorStruct*** TargetingErrorStructure;

	  //widgets
	  vtkKWRenderWidget* RenderWidget;
	  vtkKWRenderWidget* SegmentationErrorRenderWidget;
	  vtkKWRenderWidget* SegmentationErrorDifferenceRenderWidget;
	  vtkKWRenderWidget* DeformationRenderWidget;
	  vtkKWRenderWidget* DeformationSegmentationErrorRenderWidget;
	  vtkKWRenderWidget* DeformationSegmentationErrorDifferenceRenderWidget;
	  vtkKWRenderWidget* RegistrationResultRenderWidget;
	  vtkKWRenderWidget* RegistrationDifferenceRenderWidget;
	  vtkKWFrame* DoubleImageFrame;
	  vtkKWFrame* DeformationDoubleImageFrame;
	  vtkKWFrame* ImageFrame;
	  vtkKWFrame* SegmentationErrorImageFrame;
	  vtkKWFrame* SegmentationErrorDifferenceImageFrame;
	  vtkKWFrame* DeformationImageFrame;
	  vtkKWFrame* DeformationSegmentationErrorImageFrame;
	  vtkKWFrame* DeformationSegmentationErrorDifferenceImageFrame;
	  vtkKWFrame* RegistrationResultImageFrame;
	  vtkKWFrame* RegistrationDifferenceImageFrame;
	  vtkKWFrame* GlobalFrame;
	  vtkKWFrame* LabelFrame;
	  vtkKWFrame* InfoLabelFrame;
	  vtkKWFrame* ReportLabelFrame;
	  vtkKWMenu* FileMenu;
	  vtkKWMenu* DataMenu;
	  vtkKWScale* ImageScale;
	  vtkKWLabel* FirstModeLabel;
	  vtkKWLabel* SecondModeLabel;
	  vtkKWLabel* ThirdModeLabel;
	  vtkKWLabel* BaseSegmentationErrorLabel;
	  vtkKWLabel* MidSegmentationErrorLabel;
	  vtkKWLabel* ApexSegmentationErrorLabel;
	  vtkKWLabel* TotalTargetsLabel;
	  vtkKWLabel* InsignificantLabel;
	  vtkKWLabel* SignificantLabel;
	  vtkKWLabel* MinErrorLabel;
	  vtkKWLabel* MaxErrorLabel;
	  vtkKWLabel* MeanErrorLabel;
	  vtkKWLabel* MeanDeviationLabel;
	  vtkKWLabel* StandardDeviationLabel;

	  //methods
	  virtual void CreateWidget();
	  std::string GetStringFromDouble(double,int);
	  std::string GetStringFromInteger(int);
	  std::string GetTimeString(int);
	  bool FileExists(std::string);
	  OutputImageType::Pointer GetITKImageFromFile(std::string);
	  vtkImageData* GetVTKImageData(OutputImageType::Pointer, bool);
	  std::vector<Point> GetDiscretizedPath(Point startPoint, Point endPoint);
	  bool IsExtremalContourVoxel(OutputImageType::Pointer, OutputImageType::IndexType, int, int, int, int, int, int);
	  bool IsExtremalVoxel(OutputImageType::Pointer, OutputImageType::IndexType);
	  OutputImageType::Pointer Subtract(OutputImageType::Pointer,OutputImageType::Pointer);
	  OutputImageType::Pointer Add(OutputImageType::Pointer,OutputImageType::Pointer);
	  Point ConvertFromPolarCoordinates(PointPolar, Point);
	  Point GetSliceExtremalPoint(OutputImageType::Pointer, int, int, int, float);
	  std::vector<Point> GetBresenhamPath(Point, Point);
	  void FillSlice(OutputImageType::Pointer, int, int, int);
	  void GrowInitialContour(OutputImageType::Pointer);
	  void GatherContourStatistics(OutputImageType::Pointer, ContourStatistics*);
	  int GetPosteriorYCoordinate(OutputImageType::Pointer, ContourStatistics*, int, int);
	  int GetAnteriorYCoordinate(OutputImageType::Pointer, ContourStatistics*, int, int);
	  TargetPointStruct GetAllTargetPoints(int);
	  void PlaceGroundTruthPreOpTargetLocationsAll(int);
	  void PlaceGroundTruthPreOpTargetLocationsAllRandom(int);
	  void PlaceGroundTruthPreOpTargetLocations(int);
	  void DetermineGroundTruthIntraOpTargetLocations(int);
	  void DetermineResultTargetLocations(int);
	  TargetPoint GetMeanShapeCorrespondencePoint(TargetPoint, ContourStatistics*);
	  TargetPoint ConvertITKPointToVTKPoint(TargetPoint);
	  void RecordTargetingErrors(int);
	  OutputImageType::Pointer CreateContourWithSegmentationError(OutputImageType::Pointer, ContourStatistics&);
	  OutputImageType::Pointer DeformImage(OutputImageType::Pointer, DeformationFieldType::Pointer);
	  OutputImageType::Pointer CreateContourWithDeformationField(OutputImageType::Pointer, ContourStatistics&);
	  DeformationFieldType::Pointer SolveDeformableRegistrationDemons(OutputImageType::Pointer,OutputImageType::Pointer);
	  DeformationFieldType::Pointer SolveDeformableRegistrationLevelSets(OutputImageType::Pointer,OutputImageType::Pointer);
	  DeformationFieldType::Pointer SolveDeformableRegistrationBSplines(OutputImageType::Pointer,OutputImageType::Pointer);
	  bool HasInternalTargetPoints(OutputImageType::Pointer, ContourStatistics*, std::vector<TargetPoint>&);
	  void DisplayVTKImage(vtkImageData*,vtkKWRenderWidget**);
	  void DisplayVTKImage(vtkImageData*,vtkKWRenderWidget**, OutputImageType::Pointer, ContourStatistics*, std::vector<TargetPoint>&, std::vector<TargetPoint>&);
	  void SetupGUIWithMultipleImageFrames(int);
	  void ResetGUIToOneImageFrame();
	  void AddDefaultGUIObjects();
	  void AddObservers();
	  void GenerateMeanShape();
	  void LoadProstateContour();
	  bool IsSurfaceVoxel(OutputImageType::IndexType, OutputImageType::Pointer, ContourStatistics*);
	  void DisplayContour(std::string);
	  void DisplayMeanShapeWithHeatMap();
	  void DisplaySimulationCase(int);
	  void UpdateLabels(int);
	  void SetupTargetingErrorDataStructure();
	  void DestroyTargetingErrorDataStructure();
	  void RecordFinalTargetingErrorStatistics();
	  ReportStruct GetFinalReport();
	  void RunSimulation();

	  // Description:
	  // Processes the events that are passed through CallbackCommand (or others).
	  virtual void ProcessCallbackCommandEvents(vtkObject *caller, unsigned long event, void *calldata);

	private:
	  vtkKWProstateErrorMapRenderingWidget(const vtkKWProstateErrorMapRenderingWidget&);   // Not implemented.
	  void operator=(const vtkKWProstateErrorMapRenderingWidget&);  // Not implemented.
};

#endif