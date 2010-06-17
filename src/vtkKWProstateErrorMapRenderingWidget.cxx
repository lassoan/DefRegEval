#define _CRT_SECURE_NO_WARNINGS

#include "vtkKWProstateErrorMapRenderingWidget.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkKWProstateErrorMapRenderingWidget);
vtkCxxRevisionMacro(vtkKWProstateErrorMapRenderingWidget, "$Revision: 1.1 $");
//----------------------------------------------------------------------------

const char* vtkKWProstateErrorMapRenderingWidget::GENERATE_MEAN_SHAPE_MENU_STRING = "Generate Mean Shape...";
const char* vtkKWProstateErrorMapRenderingWidget::LOAD_PROSTATE_CONTOUR_MENU_STRING = "Load Prostate Contour...";
const char* vtkKWProstateErrorMapRenderingWidget::EXIT_MENU_STRING = "Exit";
const char* vtkKWProstateErrorMapRenderingWidget::RUN_SIMULATION_MENU_STRING = "Run Simulation";

const double vtkKWProstateErrorMapRenderingWidget::MIN_WEIGHT_MODE_1 = -0.2;
const double vtkKWProstateErrorMapRenderingWidget::MAX_WEIGHT_MODE_1 = 0;
const double vtkKWProstateErrorMapRenderingWidget::MIN_WEIGHT_MODE_2 = -1.5;
const double vtkKWProstateErrorMapRenderingWidget::MAX_WEIGHT_MODE_2 = 1.5;
const double vtkKWProstateErrorMapRenderingWidget::MIN_WEIGHT_MODE_3 = 0;
const double vtkKWProstateErrorMapRenderingWidget::MAX_WEIGHT_MODE_3 = 0.2;

const double vtkKWProstateErrorMapRenderingWidget::MIN_MID_SEGMENTATION_ERROR = 0; //mm
const double vtkKWProstateErrorMapRenderingWidget::MAX_MID_SEGMENTATION_ERROR = 1; //mm
const double vtkKWProstateErrorMapRenderingWidget::MIN_BASE_SEGMENTATION_ERROR = 0; //mm
const double vtkKWProstateErrorMapRenderingWidget::MAX_BASE_SEGMENTATION_ERROR = 2; //mm
const double vtkKWProstateErrorMapRenderingWidget::MIN_APEX_SEGMENTATION_ERROR = 0; //mm
const double vtkKWProstateErrorMapRenderingWidget::MAX_APEX_SEGMENTATION_ERROR = 2; //mm
const double vtkKWProstateErrorMapRenderingWidget::MAX_SEGMENTATION_ERROR_DIFFERENCE = 0.3; //mm

const double vtkKWProstateErrorMapRenderingWidget::DEFORMATION_PROBABILITY = 0.25;
const double vtkKWProstateErrorMapRenderingWidget::MIN_DEFORMATION_MAGNITUDE = 0.5; //mm
const double vtkKWProstateErrorMapRenderingWidget::MAX_DEFORMATION_MAGNITUDE = 1; //mm

const double vtkKWProstateErrorMapRenderingWidget::ANGLE_NOT_TARGETED = 120;
const double vtkKWProstateErrorMapRenderingWidget::RANDOM_TARGET_PERCENT = 50;
const double vtkKWProstateErrorMapRenderingWidget::PERIPHERAL_ZONE_PERCENT = 15;
const int vtkKWProstateErrorMapRenderingWidget::NUM_TARGETS[vtkKWProstateErrorMapRenderingWidget::BIOPSY_SAMPLING_Z] = {3, 2, 3};


const double vtkKWProstateErrorMapRenderingWidget::MAX_INSIGNIFICANT_ERROR = 2.5; //mm
const double vtkKWProstateErrorMapRenderingWidget::MIN_SIGNIFICANT_ERROR = 5; //mm

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//Constructor
vtkKWProstateErrorMapRenderingWidget::vtkKWProstateErrorMapRenderingWidget()
{
	this->RenderWidget = NULL;
	this->SegmentationErrorRenderWidget = NULL;
	this->SegmentationErrorDifferenceRenderWidget = NULL;
	this->DeformationRenderWidget = NULL;
	this->DeformationSegmentationErrorRenderWidget = NULL;
	this->DeformationSegmentationErrorDifferenceRenderWidget = NULL;
	this->RegistrationResultRenderWidget = NULL;
	this->RegistrationDifferenceRenderWidget = NULL;
	this->ImageFrame = NULL;
	this->FileMenu = NULL;
	this->DataMenu = NULL;
	this->MeanShape = NULL;
	this->ImageScale = NULL;
	this->GlobalFrame = NULL;
	this->DoubleImageFrame = NULL;
	this->DeformationDoubleImageFrame = NULL;
	this->SegmentationErrorImageFrame = NULL;
	this->SegmentationErrorDifferenceImageFrame = NULL;
	this->DeformationImageFrame = NULL;
	this->DeformationSegmentationErrorImageFrame = NULL;
	this->DeformationSegmentationErrorDifferenceImageFrame = NULL;
	this->RegistrationResultImageFrame = NULL;
	this->RegistrationDifferenceImageFrame = NULL;
	this->LabelFrame = NULL;
	this->InfoLabelFrame = NULL;
	this->ReportLabelFrame = NULL;
	this->FirstModeLabel = NULL;
	this->SecondModeLabel = NULL;
	this->ThirdModeLabel = NULL;
	this->BaseSegmentationErrorLabel = NULL;
	this->MidSegmentationErrorLabel = NULL;
	this->ApexSegmentationErrorLabel = NULL;
	this->TotalTargetsLabel = NULL;
	this->InsignificantLabel = NULL;
	this->SignificantLabel = NULL;
	this->MinErrorLabel = NULL;
	this->MaxErrorLabel = NULL;
	this->MeanErrorLabel = NULL;
	this->MeanDeviationLabel = NULL;
	this->StandardDeviationLabel = NULL;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//Destructor
vtkKWProstateErrorMapRenderingWidget::~vtkKWProstateErrorMapRenderingWidget()
{
	if(this->RenderWidget)
	{
		this->RenderWidget->SetParent(NULL);
		this->RenderWidget->Delete();
		this->RenderWidget = NULL;
	}
	if(this->SegmentationErrorRenderWidget)
	{
		this->SegmentationErrorRenderWidget->SetParent(NULL);
		this->SegmentationErrorRenderWidget->Delete();
		this->SegmentationErrorRenderWidget = NULL;
	}
	if(this->SegmentationErrorDifferenceRenderWidget)
	{
		this->SegmentationErrorDifferenceRenderWidget->SetParent(NULL);
		this->SegmentationErrorDifferenceRenderWidget->Delete();
		this->SegmentationErrorDifferenceRenderWidget = NULL;
	}
	if(this->DeformationRenderWidget)
	{
		this->DeformationRenderWidget->SetParent(NULL);
		this->DeformationRenderWidget->Delete();
		this->DeformationRenderWidget = NULL;
	}
	if(this->DeformationSegmentationErrorRenderWidget)
	{
		this->DeformationSegmentationErrorRenderWidget->SetParent(NULL);
		this->DeformationSegmentationErrorRenderWidget->Delete();
		this->DeformationSegmentationErrorRenderWidget = NULL;
	}
	if(this->DeformationSegmentationErrorDifferenceRenderWidget)
	{
		this->DeformationSegmentationErrorDifferenceRenderWidget->SetParent(NULL);
		this->DeformationSegmentationErrorDifferenceRenderWidget->Delete();
		this->DeformationSegmentationErrorDifferenceRenderWidget = NULL;
	}
	if(this->RegistrationResultRenderWidget)
	{
		this->RegistrationResultRenderWidget->SetParent(NULL);
		this->RegistrationResultRenderWidget->Delete();
		this->RegistrationResultRenderWidget = NULL;
	}
	if(this->RegistrationDifferenceRenderWidget)
	{
		this->RegistrationDifferenceRenderWidget->SetParent(NULL);
		this->RegistrationDifferenceRenderWidget->Delete();
		this->RegistrationDifferenceRenderWidget = NULL;
	}
	if(this->FileMenu)
	{
		this->FileMenu->Delete();
		this->FileMenu = NULL;
	}
	if(this->DataMenu)
	{
		this->DataMenu->Delete();
		this->DataMenu = NULL;
	}
	if(this->ImageScale)
	{
		this->ImageScale->Delete();
		this->ImageScale = NULL;
	}
	if(this->GlobalFrame)
	{
		this->GlobalFrame->Delete();
		this->GlobalFrame = NULL;
	}
	if(this->DoubleImageFrame)
	{
		this->DoubleImageFrame->Delete();
		this->DoubleImageFrame = NULL;
	}
	if(this->DeformationDoubleImageFrame)
	{
		this->DeformationDoubleImageFrame->Delete();
		this->DeformationDoubleImageFrame = NULL;
	}
	if(this->ImageFrame)
	{
		this->ImageFrame->Delete();
		this->ImageFrame = NULL;
	}
	if(this->SegmentationErrorImageFrame)
	{
		this->SegmentationErrorImageFrame->Delete();
		this->SegmentationErrorImageFrame = NULL;
	}
	if(this->SegmentationErrorDifferenceImageFrame)
	{
		this->SegmentationErrorDifferenceImageFrame->Delete();
		this->SegmentationErrorDifferenceImageFrame = NULL;
	}
	if(this->DeformationImageFrame)
	{
		this->DeformationImageFrame->Delete();
		this->DeformationImageFrame = NULL;
	}
	if(this->DeformationSegmentationErrorImageFrame)
	{
		this->DeformationSegmentationErrorImageFrame->Delete();
		this->DeformationSegmentationErrorImageFrame = NULL;
	}
	if(this->DeformationSegmentationErrorDifferenceImageFrame)
	{
		this->DeformationSegmentationErrorDifferenceImageFrame->Delete();
		this->DeformationSegmentationErrorDifferenceImageFrame = NULL;
	}
	if(this->RegistrationResultImageFrame)
	{
		this->RegistrationResultImageFrame->Delete();
		this->RegistrationResultImageFrame = NULL;
	}
	if(this->RegistrationDifferenceImageFrame)
	{
		this->RegistrationDifferenceImageFrame->Delete();
		this->RegistrationDifferenceImageFrame = NULL;
	}
	if(this->LabelFrame)
	{
		this->LabelFrame->Delete();
		this->LabelFrame = NULL;
	}
	if(this->InfoLabelFrame)
	{
		this->InfoLabelFrame->Delete();
		this->InfoLabelFrame = NULL;
	}
	if(this->ReportLabelFrame)
	{
		this->ReportLabelFrame->Delete();
		this->ReportLabelFrame = NULL;
	}
	if(this->FirstModeLabel)
	{
		this->FirstModeLabel->Delete();
		this->FirstModeLabel = NULL;
	}
	if(this->SecondModeLabel)
	{
		this->SecondModeLabel->Delete();
		this->SecondModeLabel = NULL;
	}
	if(this->ThirdModeLabel)
	{
		this->ThirdModeLabel->Delete();
		this->ThirdModeLabel = NULL;
	}
	if(this->BaseSegmentationErrorLabel)
	{
		this->BaseSegmentationErrorLabel->Delete();
		this->BaseSegmentationErrorLabel = NULL;
	}
	if(this->MidSegmentationErrorLabel)
	{
		this->MidSegmentationErrorLabel->Delete();
		this->MidSegmentationErrorLabel = NULL;
	}
	if(this->ApexSegmentationErrorLabel)
	{
		this->ApexSegmentationErrorLabel->Delete();
		this->ApexSegmentationErrorLabel = NULL;
	}
	if(this->TotalTargetsLabel)
	{
		this->TotalTargetsLabel->Delete();
		this->TotalTargetsLabel = NULL;
	}
	if(this->InsignificantLabel)
	{
		this->InsignificantLabel->Delete();
		this->InsignificantLabel = NULL;
	}
	if(this->SignificantLabel)
	{
		this->SignificantLabel->Delete();
		this->SignificantLabel = NULL;
	}
	if(this->MinErrorLabel)
	{
		this->MinErrorLabel->Delete();
		this->MinErrorLabel = NULL;
	}
	if(this->MaxErrorLabel)
	{
		this->MaxErrorLabel->Delete();
		this->MaxErrorLabel = NULL;
	}
	if(this->MeanErrorLabel)
	{
		this->MeanErrorLabel->Delete();
		this->MeanErrorLabel = NULL;
	}
	if(this->MeanDeviationLabel)
	{
		this->MeanDeviationLabel->Delete();
		this->MeanDeviationLabel = NULL;
	}
	if(this->StandardDeviationLabel)
	{
		this->StandardDeviationLabel->Delete();
		this->StandardDeviationLabel = NULL;
	}

	if(DO_REGISTRATION)
	{
		this->DestroyTargetingErrorDataStructure();
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts a double with number of decimal places "precision" to a string
std::string vtkKWProstateErrorMapRenderingWidget::GetStringFromDouble(double num, int precision)
{
	if(num != 0)
	{
		char *str;
		int  dec, sign;
		if(num < 0)
		{
			sign = -1;
		}
		else
		{
			sign = 1;
		}

		//round to precision
		float epsilon = (float)1.2/(float)power(10, precision+2);
		float epsilon2 = (float)5/(float)power(10, precision+1);
		int roundedNum = (num + sign*epsilon + sign*epsilon2) * power(10, precision);
		num = (double)roundedNum/(double)power(10, precision);

		str = _fcvt(num, precision, &dec, &sign );

		std::string beforeDecimal(str);
		std::string afterDecimal(str);

		if(dec <= 0)
		{
			beforeDecimal = "0";

			if(dec < 0)
			{
				std::string num = afterDecimal;

				afterDecimal = "";

				for(int i=1; i<= abs(dec); i++)
				{
					afterDecimal.append("0");
				}
				
				afterDecimal.append(num);
			}
			else
			{
				afterDecimal = afterDecimal.substr(dec, strlen(str) - 1 - dec + 1);
			}
		}
		else
		{
			beforeDecimal = beforeDecimal.substr(0, dec);
			afterDecimal = afterDecimal.substr(dec, strlen(str) - 1 - dec + 1);
		}
		
		
		std::string doubleString(beforeDecimal);
		doubleString.append(".");
		doubleString.append(afterDecimal);

		if(sign == 1)
		{
			std::string withMinus = "-";
			withMinus.append(doubleString);
			return withMinus;
		}
		else
		{
			return doubleString;
		}
	}
	else
	{
		return "0";
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts an integer to a string
std::string vtkKWProstateErrorMapRenderingWidget::GetStringFromInteger(int num)
{
	char buffer [1000];
	return (std::string)_itoa(num,buffer,10);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts a time, in number of seconds, to a sentence, like "5 minutes, 3 seconds."
std::string vtkKWProstateErrorMapRenderingWidget::GetTimeString(int timeElapsed)
{
	int minutes = timeElapsed/60;
	int seconds = timeElapsed - minutes*60;

	std::string timeString(this->GetStringFromInteger(minutes));
	if(minutes == 1)
	{
		timeString.append(" minute, ");
	}
	else
	{
		timeString.append(" minutes, ");
	}
	
	timeString.append(this->GetStringFromInteger(seconds));
	if(seconds == 1)
	{
		timeString.append(" second.");
	}
	else
	{
		timeString.append(" seconds.");
	}
	

	return timeString;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returs true if a file with absolute path "file" exists, false otherwise
bool vtkKWProstateErrorMapRenderingWidget::FileExists(std::string file)
{
	struct stat stFileInfo;

	// Attempt to get the file attributes
	int intStat = stat(file.c_str(),&stFileInfo);

	return (intStat == 0);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function reads-in and returns an ITK image at the absolute file "filePath"
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::GetITKImageFromFile(std::string filePath)
{
	typedef itk::ImageFileReader<OutputImageType > ImageReaderType;
	ImageReaderType::Pointer  ImageReader  = ImageReaderType::New();

	ImageReader->SetFileName(filePath.c_str());

	try
	{
		ImageReader->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		exit(-1);
	}

	return ImageReader->GetOutput();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts an ITK image to a VTK image
vtkImageData* vtkKWProstateErrorMapRenderingWidget::GetVTKImageData(OutputImageType::Pointer ITKImage, bool flipImage)
{
	//get the inverted the volume
	typedef itk::FlipImageFilter<OutputImageType>  FlipFilterType;
	typedef FlipFilterType::FlipAxesArrayType     FlipAxesArrayType;

	FlipFilterType::Pointer flipImageFilter = FlipFilterType::New();

	FlipAxesArrayType flipArray;
	flipArray[0] = 0;
	flipArray[1] = 1;
	flipImageFilter->SetFlipAxes(flipArray);

	flipImageFilter->SetInput(ITKImage);

	OutputImageType::Pointer flippedImage = flipImageFilter->GetOutput();

	
	//convert the volume to VTK
	this->ITKImageToVTKFilter = ITKImageToVTKFilterType::New();
	if(flipImage)
	{
		this->ITKImageToVTKFilter->SetInput(flippedImage);
	}
	else
	{
		this->ITKImageToVTKFilter->SetInput(ITKImage);
	}
	

	try
	{
		this->ITKImageToVTKFilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "Error: "<<err<< std::endl;
		exit(-1);
	}

	return this->ITKImageToVTKFilter->GetOutput();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function gets the discretized path from the point startPoint to the point endPoint. This path is as if someone would
//be walking on a grid.
std::vector<Point> vtkKWProstateErrorMapRenderingWidget::GetDiscretizedPath(Point startPoint, Point endPoint)
{
	std::vector<Point> path;
	path.push_back(startPoint);

	int difX = endPoint.x - startPoint.x;
	int difY = endPoint.y - startPoint.y;

	int signDifX = 0;
	int signDifY = 0;

	if(difX != 0)
	{
		signDifX = difX/abs(difX);
	}
	if(difY != 0)
	{
		signDifY = difY/abs(difY);
	}

	if(difX == 0)
	{
		int y = startPoint.y + signDifY;

		while(y != endPoint.y)
		{
			Point pointInPath;
			pointInPath.x = startPoint.x;
			pointInPath.y = y;

			path.push_back(pointInPath);

			y = y + signDifY;
		}
	}
	else if(difY == 0)
	{
		int x = startPoint.x + signDifX;

		while(x != endPoint.x)
		{
			Point pointInPath;
			pointInPath.x = x;
			pointInPath.y = startPoint.y;

			path.push_back(pointInPath);

			x = x + signDifX;
		}
	}
	else
	{
		//get line going from the center of the start voxel to the
		//center of the end voxel
		double linePointX = startPoint.x + signDifX*0.5;
		double linePointY = startPoint.y + signDifY*0.5;

		double lineVectorX = endPoint.x + signDifX*0.5 - linePointX;
		double lineVectorY = endPoint.y + signDifY*0.5 - linePointY;

		//get vector going orthogonal to the line
		lineVectorY = (double)lineVectorY/(double)lineVectorX;
		lineVectorY = (double)-1/(double)lineVectorY;
		lineVectorX = 1;
		double lineVectorLength = sqrt(lineVectorX*lineVectorX + lineVectorY*lineVectorY);

		for(int j=0; j<=abs(difY); j++)
		{
			for(int i=0; i<=abs(difX); i++)
			{
				if(((i!=0) || (j!=0)) && ((i!=abs(difX)) || (j!=abs(difY))))
				{
					Point cornerPoint1;
					Point cornerPoint2;

					cornerPoint1.x = startPoint.x + signDifX*i + signDifX;
					cornerPoint1.y = startPoint.y + signDifY*j;

					cornerPoint2.x = cornerPoint1.x - signDifX;
					cornerPoint2.y = cornerPoint1.y + signDifY;

					//test if the corner points are on opposite sides of the line

					double vector1X = cornerPoint1.x - linePointX;
					double vector1Y = cornerPoint1.y - linePointY;
					double vector2X = cornerPoint2.x - linePointX;
					double vector2Y = cornerPoint2.y - linePointY;

					double signedDistanceToLine1 = (double)(vector1X*lineVectorX + vector1Y*lineVectorY)/(double)lineVectorLength;
					double signedDistanceToLine2 = (double)(vector2X*lineVectorX + vector2Y*lineVectorY)/(double)lineVectorLength;

					if(((signedDistanceToLine1 < 0) && (signedDistanceToLine2 > 0)) || ((signedDistanceToLine1 > 0) && (signedDistanceToLine2 < 0)))
					{
						//line goes through this voxel, voxel is in the path
						Point pointInPath;
						pointInPath.x = startPoint.x + signDifX*i;
						pointInPath.y = startPoint.y + signDifY*j;
						path.push_back(pointInPath);
					}
				}
			}
		}
	}

	path.push_back(endPoint);

	return path;

}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns true if the point at index "index" in the ITK image "image" is an extremal contour voxel (i.e. it is on the surface of
//the contour)
bool vtkKWProstateErrorMapRenderingWidget::IsExtremalContourVoxel(OutputImageType::Pointer image, OutputImageType::IndexType index, int minX, int maxX, int minY, int maxY, int minZ, int maxZ)
{
	int x = index[0];
	int y = index[1];
	int z = index[2];

	bool thisVoxel = (bool)image->GetPixel(index);

	if(!thisVoxel)
	{
		return false;
	}

	//left voxel
	OutputImageType::IndexType leftIndex;
	leftIndex[0] = x - 1;
	leftIndex[1] = y;
	leftIndex[2] = z;

	//right voxel
	OutputImageType::IndexType rightIndex;
	rightIndex[0] = x + 1;
	rightIndex[1] = y;
	rightIndex[2] = z;

	//up voxel
	OutputImageType::IndexType upIndex;
	upIndex[0] = x;
	upIndex[1] = y + 1;
	upIndex[2] = z;

	//down voxel
	OutputImageType::IndexType downIndex;
	downIndex[0] = x;
	downIndex[1] = y - 1;
	downIndex[2] = z;

	bool voxelLeft = (bool)image->GetPixel(leftIndex);
	bool voxelRight = (bool)image->GetPixel(rightIndex);
	bool voxelUp = (bool)image->GetPixel(upIndex);
	bool voxelDown = (bool)image->GetPixel(downIndex);

	if(voxelLeft && voxelRight && voxelUp && voxelDown)
	{
		return false;
	}

	bool alongX = false;
	int end = 0;
	int start = 0;

	if(!voxelLeft)
	{
		start = index[0];
		end = minX;
		alongX = true;
	}
	else if(!voxelRight)
	{
		start = index[0];
		end = maxX;
		alongX = true;
	}
	else if(!voxelUp)
	{
		start = index[1];
		end = maxY;
	}
	else if(!voxelDown)
	{
		start = index[1];
		end = minY;
	}

	int dif = end - start;
	if(dif < 0)
	{
		dif = -1;
	}
	else
	{
		dif = 1;
	}

	for(int i=1; i<=abs(end - start); i++)
	{
		OutputImageType::IndexType currentIndex;

		if(alongX)
		{
			currentIndex[0] = start + i*dif;
			currentIndex[1] = y;
		}
		else
		{
			currentIndex[0] = x;
			currentIndex[1] = start + i*dif;
		}

		currentIndex[2] = z;

		bool currentVoxel = (bool)image->GetPixel(currentIndex);

		if(currentVoxel)
		{
			return false;
		}
	}

	return true;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns true if the point at index "index" in the ITK image "contour" is an extremal contour voxel (i.e. it is on the surface of
//the contour)

bool vtkKWProstateErrorMapRenderingWidget::IsExtremalVoxel(OutputImageType::Pointer contour, OutputImageType::IndexType index)
{
	if(((int) contour->GetPixel(index)) == 0)
	{
		return false;
	}

	OutputImageType::IndexType checkIndex;

	//right voxel
	checkIndex[0] = index[0] + 1;
	checkIndex[1] = index[1];
	checkIndex[2] = index[2];

	if(((int) contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//left voxel
	checkIndex[0] = index[0] - 1;
	checkIndex[1] = index[1];
	checkIndex[2] = index[2];

	if(((int) contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//up voxel
	checkIndex[0] = index[0];
	checkIndex[1] = index[1] + 1;
	checkIndex[2] = index[2];

	if(((int) contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//down voxel
	checkIndex[0] = index[0];
	checkIndex[1] = index[1] - 1;
	checkIndex[2] = index[2];

	if(((int) contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	return false;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns the image "image1 - image2"
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::Subtract(OutputImageType::Pointer image1, OutputImageType::Pointer image2)
{
	typedef itk::SubtractImageFilter<OutputImageType, OutputImageType, OutputImageType> SubtractFilterType;
    SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();

	subtractFilter->SetInput1(image1);
	subtractFilter->SetInput2(image2);

	subtractFilter->Update();

	return subtractFilter->GetOutput();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns the image "image1 + image2"
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::Add(OutputImageType::Pointer image1, OutputImageType::Pointer image2)
{
	typedef itk::AddImageFilter<OutputImageType,OutputImageType,OutputImageType> AddFilterType;
    AddFilterType::Pointer addFilter = AddFilterType::New();

	addFilter->SetInput1(image1);
	addFilter->SetInput2(image2);

	addFilter->Update();

	return addFilter->GetOutput();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns the surface point of the slice "z" of the ITK image "image" at angle "angle" (in radians) from the gravity point
Point vtkKWProstateErrorMapRenderingWidget::GetSliceExtremalPoint(OutputImageType::Pointer image, int gravityX, int gravityY, int z, float angle)
{
	float stepX = cos(angle);
	float stepY = sin(angle);

	float x = gravityX;
	float y = gravityY;

	int maxSteps = 100;
	int steps = 0;

	OutputImageType::IndexType index;
	Point sliceExtremalPoint;

	while(steps < maxSteps)
	{

		x = x + stepX;
		y = y + stepY;

		index[0] = (int)floor(x);
		index[1] = (int)floor(y);
		index[2] = z;

		if(((int)image->GetPixel(index)) == 0)
		{
			
			sliceExtremalPoint.x = (int) floor(x - stepX);
			sliceExtremalPoint.y = (int) floor(y - stepY);
			break;
		}
		
		steps++;
	}

	return sliceExtremalPoint;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns the bresenham path from point p1 to point p2 (these points aren't tied to any particular contour)
std::vector<Point> vtkKWProstateErrorMapRenderingWidget::GetBresenhamPath(Point p1, Point p2)
{
	std::vector<Point> bresenhamPath;

	int x0 = p1.x;
	int y0 = p1.y;
	int x1 = p2.x;
	int y1 = p2.y;

	int deltaX = x1 - x0;
	int deltaY = y1 - y0;

	Point p;

	if(deltaX == 0)
	{
		if(deltaY != 0)
		{
			int stepY = deltaY/abs(deltaY);
			int y = y0;
			while(y != y1)
			{
			  y = y + stepY;
			  
			  p.x = p1.x;
			  p.y = y;
			  bresenhamPath.push_back(p);
			}
		}
	}
	else if(deltaY == 0)
	{
		if(deltaX != 0)
		{
			int stepX = deltaX/abs(deltaX);
			int x = x0;
			while(x != x1)
			{
			  x = x + stepX;
			  
			  p.x = x;
			  p.y = p1.y;
			  bresenhamPath.push_back(p);
			}
		}
	}
	else if(abs(deltaY) > abs(deltaX))
	{
		float x = x0 + 0.5;
		int y = y0;

		int stepY = deltaY / abs(deltaY);
		float stepX = (float)deltaX / (float)abs(deltaY);

		while(y != y1)
		{
		  x = x + stepX;
		  y = y + stepY;
		  
		  p.x = x;
		  p.y = y;
		  bresenhamPath.push_back(p);
		}
	}
	else
	{
		float y = y0 + 0.5;
		int x = x0;

		int stepX = deltaX / abs(deltaX);
		float stepY = (float)deltaY / (float)abs(deltaX);

		while(x != x1)
		{
		  x = x + stepX;
		  y = y + stepY;
		  
		  p.x = x;
		  p.y = y;
		  bresenhamPath.push_back(p);
		}
	}

	return bresenhamPath;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts a 2D point "pointXY" with origin "origin" to polar coordinates
PointPolar vtkKWProstateErrorMapRenderingWidget::ConvertToPolarCoordinates(Point pointXY, Point origin)
{
	float xDif = pointXY.x - origin.x;
	float yDif = pointXY.y - origin.y;

	PointPolar pointInPolarCoordinates;

	pointInPolarCoordinates.r = sqrt(xDif*xDif + yDif*yDif);
	pointInPolarCoordinates.theta = atan(abs((float)yDif/(float)xDif));

	if((xDif < 0) && (yDif >= 0))
	{
		//second quadrant
		pointInPolarCoordinates.theta = PI - pointInPolarCoordinates.theta;
	}
	else if((xDif <= 0) && (yDif < 0))
	{
		//third quadrant
		pointInPolarCoordinates.theta = PI + pointInPolarCoordinates.theta;
	}
	else if((xDif > 0) && (yDif < 0))
	{
		//fourth quadrant
		pointInPolarCoordinates.theta = 2*PI - pointInPolarCoordinates.theta;
	}

	return pointInPolarCoordinates;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts a 2D point in polar coordinates "pointPolar" with origin "origin" to cartesian coordinates
Point vtkKWProstateErrorMapRenderingWidget::ConvertFromPolarCoordinates(PointPolar pointPolar, Point origin)
{
	float x = pointPolar.r*cos(pointPolar.theta);
	float y = pointPolar.r*sin(pointPolar.theta);

	Point pointXY;
	pointXY.x = floor(x) + origin.x;
	pointXY.y = floor(y) + origin.y;

	return pointXY;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function fills a slice of the ITK contour "image" starting at the point (x, y, z). This fills the contour
//in the x-y plane.
void vtkKWProstateErrorMapRenderingWidget::FillSlice(OutputImageType::Pointer image, int x, int y, int z)
{
	OutputImageType::IndexType index;

	index[0] = x;
	index[1] = y;
	index[2] = z;

	image->SetPixel(index, 1);

	OutputImageType::IndexType indexUp;
	indexUp[0] = x;
	indexUp[1] = y - 1;
	indexUp[2] = z;

	if(((int)image->GetPixel(indexUp)) == 0)
	{
		this->FillSlice(image, x, y - 1, z);
	}

	OutputImageType::IndexType indexRight;
	indexRight[0] = x + 1;
	indexRight[1] = y;
	indexRight[2] = z;

	if(((int)image->GetPixel(indexRight)) == 0)
	{
		this->FillSlice(image, x + 1, y, z);
	}

	OutputImageType::IndexType indexDown;
	indexDown[0] = x;
	indexDown[1] = y + 1;
	indexDown[2] = z;

	if(((int)image->GetPixel(indexDown)) == 0)
	{
		this->FillSlice(image, x, y + 1, z);
	}

	OutputImageType::IndexType indexLeft;
	indexLeft[0] = x - 1;
	indexLeft[1] = y;
	indexLeft[2] = z;

	if(((int)image->GetPixel(indexLeft)) == 0)
	{
		this->FillSlice(image, x - 1, y, z);
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function adds padding around a contour (to make space for adding deformation and segmentation error)
void vtkKWProstateErrorMapRenderingWidget::GrowInitialContour(OutputImageType::Pointer initialContour)
{
	//set up iterator to loop through the image
	OutputImageType::Pointer oldImage = initialContour;
	typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
	IteratorType oldImageIterator(oldImage, oldImage->GetLargestPossibleRegion());

	//create a new empty contour with size slightly bigger than the original
	OutputImageType::Pointer newImage = OutputImageType::New();

	//get the size of the original image
	OutputImageType::SizeType originalSize = oldImage->GetLargestPossibleRegion().GetSize();

	//set the size to be bigger
	OutputImageType::SizeType newSize;
	newSize[0] = originalSize[0] + 2 * PROSTATE_PADDING;
	newSize[1] = originalSize[1] + 2 * PROSTATE_PADDING;
	newSize[2] = originalSize[2] + 2 * PROSTATE_PADDING;

	//get the origin of the original image
	OutputImageType::PointType origin = oldImage->GetOrigin();

	OutputImageType::IndexType start;
	start[0] = origin[0];
	start[1] = origin[1];
	start[2] = origin[2];

	//set up the region
	OutputImageType::RegionType region;
	region.SetSize(newSize);
	region.SetIndex(start);

	newImage->SetRegions(region);
	newImage->Allocate();

	//set spacing
	newImage->SetSpacing(oldImage->GetSpacing());

	//set origin
	newImage->SetOrigin(oldImage->GetOrigin());

	int originOffset[3];
	originOffset[0] = PROSTATE_PADDING;
	originOffset[1] = PROSTATE_PADDING;
	originOffset[2] = PROSTATE_PADDING;

	IteratorType newImageIterator(newImage, newImage->GetLargestPossibleRegion());

	//initialize all pixels to zero
	for(newImageIterator.GoToBegin(); !newImageIterator.IsAtEnd(); ++newImageIterator)
	{
		OutputImageType::IndexType index = newImageIterator.GetIndex();
		OutputImageType::PixelType pixel = newImageIterator.Get();

		newImage->SetPixel(index, 0);
	}

	//fill in all the white pixels
	for(oldImageIterator.GoToBegin(); !oldImageIterator.IsAtEnd(); ++oldImageIterator)
	{
		OutputImageType::IndexType index = oldImageIterator.GetIndex();
		OutputImageType::PixelType pixel = oldImageIterator.Get();

		if(((int)pixel) == 1)
		{
			OutputImageType::IndexType newIndex;
			newIndex[0] = index[0] + originOffset[0];
			newIndex[1] = index[1] + originOffset[1];
			newIndex[2] = index[2] + originOffset[2];

			newImage->SetPixel(newIndex, 1);
		}
	}

	this->ITKContours[this->ITKContours.size() - 1] = newImage;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function gathers contour statistics such as ranges for x,y,z, gravity centers, and surface points for the ITK image "contour"
void vtkKWProstateErrorMapRenderingWidget::GatherContourStatistics(OutputImageType::Pointer contour, ContourStatistics* statistics)
{
	//set up iterator to loop through the contour
	typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
	IteratorType iterator(contour, contour->GetLargestPossibleRegion());

	//gather contour statistics
	int contourPixels = 0;

	int totalX = 0;
	int totalY = 0;
	int totalZ = 0;
	int minX = 10000;
	int minY = 10000;
	int minZ = 10000;
	int maxX = -1;
	int maxY = -1;
	int maxZ = -1;

	for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
	{
		OutputImageType::IndexType index = iterator.GetIndex();
		OutputImageType::PixelType pixel = iterator.Get();

		if(((int)pixel) == 1)
		{
			contourPixels++;

			int x = index[0];
			int y = index[1];
			int z = index[2];

			totalX = totalX + x;
			totalY = totalY + y;
			totalZ = totalZ + z;

			if(x < minX)
			{
				minX = x;
			}
			if(x > maxX)
			{
				maxX = x;
			}

			if(y < minY)
			{
				minY = y;
			}
			if(y > maxY)
			{
				maxY = y;
			}

			if(z < minZ)
			{
				minZ = z;
			}
			if(z > maxZ)
			{
				maxZ = z;
			}
		}
	}

	int xRange = maxX - minX + 1;
	int yRange = maxY - minY + 1;
	int zRange = maxZ - minZ + 1;
	int prostateLength = zRange;
	int halfProstateLength = zRange/2;
	int midX = minX + (maxX - minX)/2;
	int midY = minY + (maxY - minY)/2;
	int midZ = minZ + halfProstateLength;
	int apexMinZ = maxZ - (double)PERCENT_APEX/(double)100 * (zRange);
	int baseMaxZ = minZ + (double)PERCENT_BASE/(double)100 * (zRange);

	//get slice statistics
	int numSlices = maxZ - minZ + 1;
	std::vector<int> sliceMinX;
	std::vector<int> sliceMaxX;
	std::vector<int> sliceMinY;
	std::vector<int> sliceMaxY;

	//initialize slice statistics
	for(int i=0; i<numSlices; i++)
	{
		sliceMinX.push_back(10000);
		sliceMinY.push_back(10000);
		sliceMaxX.push_back(-1);
		sliceMaxY.push_back(-1);
	}

	//initialize min max pairs
	ContourStatistics::MinMaxPair** minMaxZXY;
	ContourStatistics::MinMaxPair** minMaxYXZ;
	ContourStatistics::MinMaxPair** minMaxXYZ;

	minMaxZXY = new ContourStatistics::MinMaxPair*[xRange];
	minMaxYXZ = new ContourStatistics::MinMaxPair*[xRange];
	for(int x=0; x<xRange; x++)
	{
		minMaxZXY[x] = new ContourStatistics::MinMaxPair[yRange];
		minMaxYXZ[x] = new ContourStatistics::MinMaxPair[zRange];

		for(int y=0; y<yRange; y++)
		{
			minMaxZXY[x][y].min = minZ;
			minMaxZXY[x][y].max = maxZ;
		}

		for(int z=0; z<zRange; z++)
		{
			minMaxYXZ[x][z].min = minY;
			minMaxYXZ[x][z].max = maxY;
		}
	}

	minMaxXYZ = new ContourStatistics::MinMaxPair*[yRange];
	for(int y=0; y<yRange; y++)
	{
		minMaxXYZ[y] = new ContourStatistics::MinMaxPair[zRange];

		for(int z=0; z<zRange; z++)
		{
			minMaxXYZ[y][z].min = minX;
			minMaxXYZ[y][z].max = maxX;
		}
	}

	//initialize surface points
	std::vector<ContourStatistics::Point>* surfacePointsPerSlice = new std::vector<ContourStatistics::Point>[numSlices];

	for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
	{
		OutputImageType::IndexType index = iterator.GetIndex();
		OutputImageType::PixelType pixel = iterator.Get();

		if(((int)pixel) == 1)
		{
			int x = index[0];
			int y = index[1];
			int z = index[2];
			
			if(x < sliceMinX[z - minZ])
			{
				sliceMinX[z - minZ] = x;
			}
			if(x > sliceMaxX[z - minZ])
			{
				sliceMaxX[z - minZ] = x;
			}

			if(y < sliceMinY[z - minZ])
			{
				sliceMinY[z - minZ] = y;
			}
			if(y > sliceMaxY[z - minZ])
			{
				sliceMaxY[z - minZ] = y;
			}

			if((x >= minX) && (x <= maxX) && (y >= minY) && (y <= maxY) && (z >= minZ) && (z <= maxZ))
			{
				if(x < minMaxXYZ[y - minY][z - minZ].min)
				{
					minMaxXYZ[y - minY][z - minZ].min = x;
				}
				if(x > minMaxXYZ[y - minY][z - minZ].max)
				{
					minMaxXYZ[y - minY][z - minZ].max = x;
				}
			}

			if((x >= minX) && (x <= maxX) && (y >= minY) && (y <= maxY) && (z >= minZ) && (z <= maxZ))
			{
				if(y < minMaxYXZ[x - minX][z - minZ].min)
				{
					minMaxYXZ[x - minX][z - minZ].min = y;
				}
				if(y > minMaxYXZ[x - minX][z - minZ].max)
				{
					minMaxYXZ[x - minX][z - minZ].max = y;
				}
			}

			if((x >= minX) && (x <= maxX) && (y >= minY) && (y <= maxY) && (z >= minZ) && (z <= maxZ))
			{
				if(z < minMaxZXY[x - minX][y - minY].min)
				{
					minMaxZXY[x - minX][y - minY].min = z;
				}
				if(z > minMaxZXY[x - minX][y - minY].max)
				{
					minMaxZXY[x - minX][y - minY].max = z;
				}
			}

			if(this->IsExtremalVoxel(contour, index))
			{
				ContourStatistics::Point p;
				p.x = x;
				p.y = y;
				p.z = z;
				surfacePointsPerSlice[z - minZ].push_back(p);
			}
		}
	}

	//store surface points in vector
	std::vector<std::vector<ContourStatistics::Point>> surfacePoints;
	for(int z=0; z<numSlices; z++)
	{
		surfacePoints.push_back(surfacePointsPerSlice[z]);
	}

	statistics->setNumSlices(numSlices);
	statistics->setMinX(minX);
	statistics->setMaxX(maxX);
	statistics->setMinY(minY);
	statistics->setMaxY(maxY);
	statistics->setMinZ(minZ);
	statistics->setMaxZ(maxZ);
	statistics->setMinXForSlices(sliceMinX);
	statistics->setMaxXForSlices(sliceMaxX);
	statistics->setMinYForSlices(sliceMinY);
	statistics->setMaxYForSlices(sliceMaxY);
	statistics->determineGravityCenterForSlices();
	statistics->setMinMaxXYZ(minMaxXYZ);
	statistics->setMinMaxZXY(minMaxZXY);
	statistics->setMinMaxYXZ(minMaxYXZ);
	statistics->setSurfacePoints(surfacePoints);

	for(int x=0; x<xRange; x++)
	{
		delete [] minMaxZXY[x];
		delete [] minMaxYXZ[x];
	}
	delete [] minMaxZXY;
	delete [] minMaxYXZ;

	for(int y=0; y<yRange; y++)
	{
		delete [] minMaxXYZ[y];
	}
	delete [] minMaxXYZ;

	delete [] surfacePointsPerSlice;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function converts a point in ITK to a point in VTK (for display)
TargetPoint vtkKWProstateErrorMapRenderingWidget::ConvertITKPointToVTKPoint(TargetPoint p_ITK)
{
	TargetPoint p_VTK;
	p_VTK.x = p_ITK.x/1.8;
	p_VTK.y = -p_ITK.y/2 - 8;
	
	#ifdef _DEBUG
		p_VTK.z = -3*(p_ITK.z);
	#else
		p_VTK.z = 3*(p_ITK.z);
	#endif

	return p_VTK;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function finds the y-coordinate on the posterior side of the ITK image "contour" given the x and z coordinates
int vtkKWProstateErrorMapRenderingWidget::GetPosteriorYCoordinate(OutputImageType::Pointer contour, ContourStatistics* statistics, int x, int z)
{
	const int prostateVoxelThreshold = 3;
	const int loopThreshold = 200;

	OutputImageType::IndexType index;

	int y = statistics->getMinYForSlice(z - statistics->getMinZ());

	index[0] = x;
	index[1] = y;
	index[2] = z;

	if(((int)contour->GetPixel(index)) == 0)
	{
		//go down until the surface is reached
		bool contourReached = false;

		while(!contourReached)
		{
			y++; //y-axis faces down in ITK
			index[1] = y;
			contourReached = (((int)contour->GetPixel(index)) == 1);
		}
	}

	//go down until the surface is passed
	bool contourPassed = false;

	int numProstateVoxels = 1;
	int lastProstateY = y;
	int loops = 0;

	while(!contourPassed)
	{
		y++; //y-axis faces down in ITK
		index[1] = y;
		int pixelValue = (int)contour->GetPixel(index);
		contourPassed = (pixelValue == 0)  && (numProstateVoxels >= prostateVoxelThreshold);

		if(pixelValue == 1)
		{
			numProstateVoxels++;
			lastProstateY = y;
		}

		loops++;

		if(loops >= loopThreshold)
		{
			break;
		}
	}

	y--;

	if(loops >= loopThreshold)
	{
		y = lastProstateY;
	}

	return y;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function finds the y-coordinate on the anterior side of the ITK image "contour" given the x and z coordinates
int vtkKWProstateErrorMapRenderingWidget::GetAnteriorYCoordinate(OutputImageType::Pointer contour, ContourStatistics* statistics, int x, int z)
{
	const int prostateVoxelThreshold = 3;
	const int loopThreshold = 200;

	OutputImageType::IndexType index;

	int y = statistics->getMaxYForSlice(z - statistics->getMinZ());

	index[0] = x;
	index[1] = y;
	index[2] = z;

	if(((int)contour->GetPixel(index)) == 0)
	{
		//go up until the surface is reached
		bool contourReached = false;

		while(!contourReached)
		{
			y--; //y-axis faces down in ITK
			index[1] = y;
			contourReached = (((int)contour->GetPixel(index)) == 1);
		}
	}

	//go up until the surface is passed
	bool contourPassed = false;

	int numProstateVoxels = 1;
	int lastProstateY = y;
	int loops = 0;

	while(!contourPassed)
	{
		y--; //y-axis faces down in ITK
		index[1] = y;
		int pixelValue = (int)contour->GetPixel(index);
		contourPassed = (pixelValue == 0)  && (numProstateVoxels >= prostateVoxelThreshold);

		if(pixelValue == 1)
		{
			numProstateVoxels++;
			lastProstateY = y;
		}

		loops++;

		if(loops >= loopThreshold)
		{
			break;
		}
	}

	y++;

	if(loops >= loopThreshold)
	{
		y = lastProstateY;
	}

	return y;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns all possible target points for the contour at index "contourIndex". If TARGET_SURFACE is chosen, all surface
//points are returned. Otherwise, all points are returned. This selection excludes any point in the sweep of angle ANGLE_NOT_TARGETED
//above the gravity center of each slice.
vtkKWProstateErrorMapRenderingWidget::TargetPointStruct vtkKWProstateErrorMapRenderingWidget::GetAllTargetPoints(int contourIndex)
{
	TargetPointStruct allTargetPoints;

	OutputImageType::Pointer contour = this->ITKContours.at(contourIndex);
	OutputImageType::IndexType index;

	//get contour statistics
	ContourStatistics* statistics = &this->GroundTruthContourStatistics.at(contourIndex);

	int numSlices = statistics->getNumSlices();
	int minZ = statistics->getMinZ();

	double angleStart = ((double)3*PI/(double)2) - ((double)(TO_RADIANS * ANGLE_NOT_TARGETED)/(double)2);
	double angleEnd = ((double)3*PI/(double)2) + ((double)(TO_RADIANS * ANGLE_NOT_TARGETED)/(double)2);

	if(TARGET_SURFACE)
	{
		for(int z=0; z<numSlices; z++)
		{
			std::vector<ContourStatistics::Point> surfacePoints = statistics->getSurfacePointsForSlice(z);

			for(int p=0; p<surfacePoints.size(); p++)
			{
				ContourStatistics::Point surfacePoint = surfacePoints.at(p);
				double angle = statistics->angleOfSurfacePoint(surfacePoint, z);

				if((angle <= angleStart) || (angle >= angleEnd))
				{
					TargetPoint p_ITK;
					p_ITK.x = surfacePoint.x;
					p_ITK.y = surfacePoint.y;
					p_ITK.z = z + minZ;

					TargetPoint p_VTK = this->ConvertITKPointToVTKPoint(p_ITK);

					allTargetPoints.targetPointsITK.push_back(p_ITK);
					allTargetPoints.targetPointsVTK.push_back(p_VTK);
				}
			}
		}
	}
	else
	{
		for(int z=0; z<numSlices; z++)
		{
			int minXForSlice = statistics->getMinXForSlice(z);
			int minYForSlice = statistics->getMinYForSlice(z);
			int xRange = statistics->getMaxXForSlice(z) - minXForSlice + 1;
			int yRange = statistics->getMaxYForSlice(z) - minYForSlice + 1;

			bool** targeted = new bool*[xRange];
			for(int x=0; x<xRange; x++)
			{
				targeted[x] = new bool[yRange];

				for(int y=0; y<yRange; y++)
				{
					targeted[x][y] = false;
				}
			}

			Point gravityPoint;
			gravityPoint.x = statistics->getGravityXForSlice(z);
			gravityPoint.y = statistics->getGravityYForSlice(z);

			std::vector<ContourStatistics::Point> surfacePoints = statistics->getSurfacePointsForSlice(z);

			for(int p=0; p<surfacePoints.size(); p++)
			{
				ContourStatistics::Point surfacePoint = surfacePoints.at(p);
				double angle = statistics->angleOfSurfacePoint(surfacePoint, z);

				if((angle <= angleStart) || (angle >= angleEnd))
				{
					Point surfacePoint2D;
					surfacePoint2D.x = surfacePoint.x;
					surfacePoint2D.y = surfacePoint.y;

					std::vector<Point> bresenhamPath = this->GetBresenhamPath(gravityPoint, surfacePoint2D);

					for(int i=0; i<bresenhamPath.size(); i++)
					{
						Point checkTargetPoint = bresenhamPath.at(i);
						
						if(!targeted[checkTargetPoint.x - minXForSlice][checkTargetPoint.y - minYForSlice])
						{
							TargetPoint p_ITK;
							p_ITK.x = checkTargetPoint.x;
							p_ITK.y = checkTargetPoint.y;
							p_ITK.z = z + minZ;

							TargetPoint p_VTK = this->ConvertITKPointToVTKPoint(p_ITK);

							allTargetPoints.targetPointsITK.push_back(p_ITK);
							allTargetPoints.targetPointsVTK.push_back(p_VTK);

							targeted[checkTargetPoint.x - minXForSlice][checkTargetPoint.y - minYForSlice] = true;

							//
							//fill in voxels not hit by the bresenham path
							//
							//try a few voxels below
							//
							for(int i=1; i<=2; i++)
							{
								p_ITK.y = p_ITK.y + i;

								index[0] = p_ITK.x;
								index[1] = p_ITK.y;
								index[2] = p_ITK.z;

								if((p_ITK.y < statistics->getMaxY(p_ITK.x, p_ITK.z)) && (((int)contour->GetPixel(index)) == 1) && !(targeted[p_ITK.x - minXForSlice][p_ITK.y - minYForSlice]))
								{
									p_VTK = this->ConvertITKPointToVTKPoint(p_ITK);

									allTargetPoints.targetPointsITK.push_back(p_ITK);
									allTargetPoints.targetPointsVTK.push_back(p_VTK);

									targeted[p_ITK.x - minXForSlice][p_ITK.y - minYForSlice] = true;
								}
							}
						}
					}
				}
			}

			for(int x=0; x<xRange; x++)
			{
				delete [] targeted[x];
			}
			delete [] targeted;
		}
	}
	

	return allTargetPoints;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function places the target points for a particular contour as "all random". i.e. a certain percentage of all points in the contour
//are added to the target location vectors for the contour at index "contourIndex"
void vtkKWProstateErrorMapRenderingWidget::PlaceGroundTruthPreOpTargetLocationsAllRandom(int contourIndex)
{
	TargetPointStruct allTargetPoints = this->GetAllTargetPoints(contourIndex);

	int numTargetPoints = ((double)RANDOM_TARGET_PERCENT/(double)100)*allTargetPoints.targetPointsITK.size();

	std::vector<TargetPoint> currentGroundTruthPreOpTargetLocationsITK;
	std::vector<TargetPoint> currentGroundTruthPreOpTargetLocationsVTK;

	Random r;
	for(int p=0; p<numTargetPoints; p++)
	{
		int randomIndex = r.randomInt(0, allTargetPoints.targetPointsITK.size() - 1);
		currentGroundTruthPreOpTargetLocationsITK.push_back(allTargetPoints.targetPointsITK.at(randomIndex));
		currentGroundTruthPreOpTargetLocationsVTK.push_back(allTargetPoints.targetPointsVTK.at(randomIndex));

		//remove the points that were added
		allTargetPoints.targetPointsITK.erase(allTargetPoints.targetPointsITK.begin() + randomIndex);
		allTargetPoints.targetPointsVTK.erase(allTargetPoints.targetPointsVTK.begin() + randomIndex);
	}

	this->GroundTruthPreOpTargetLocationsITK.push_back(currentGroundTruthPreOpTargetLocationsITK);
	this->GroundTruthPreOpTargetLocationsVTK.push_back(currentGroundTruthPreOpTargetLocationsVTK);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function places all possible target points (excludes any point in the sweep of angle ANGLE_NOT_TARGETED
//above the gravity center of each slice) into the target location vectors for the contour at index "contourIndex"
void vtkKWProstateErrorMapRenderingWidget::PlaceGroundTruthPreOpTargetLocationsAll(int contourIndex)
{
	TargetPointStruct allTargetPoints = this->GetAllTargetPoints(contourIndex);

	this->GroundTruthPreOpTargetLocationsITK.push_back(allTargetPoints.targetPointsITK);
	this->GroundTruthPreOpTargetLocationsVTK.push_back(allTargetPoints.targetPointsVTK);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function places target points SYSTEMATICALLY on the contour at index "contourIndex"
void vtkKWProstateErrorMapRenderingWidget::PlaceGroundTruthPreOpTargetLocations(int contourIndex)
{
	OutputImageType::Pointer contour = this->ITKContours.at(contourIndex);

	//get contour statistics
	ContourStatistics* statistics = &this->GroundTruthContourStatistics.at(contourIndex);

	std::vector<TargetPoint> currentGroundTruthPreOpTargetLocationsITK;
	std::vector<TargetPoint> currentGroundTruthPreOpTargetLocationsVTK;

	int minX = statistics->getMinX();
	int maxX = statistics->getMaxX();
	int minY = statistics->getMinY();
	int minZ = statistics->getMinZ();
	int maxZ = statistics->getMaxZ();
	int numSlices = statistics->getNumSlices();

	int** regionsMinX;
	int** regionsMaxX;
	int** regionsMinZ;
	int** regionsMaxZ;
	Random r;

	//determine boundaries for regions
	regionsMinX = new int*[BIOPSY_SAMPLING_Z];
	regionsMaxX = new int*[BIOPSY_SAMPLING_Z];
	regionsMinZ = new int*[BIOPSY_SAMPLING_Z];
	regionsMaxZ = new int*[BIOPSY_SAMPLING_Z];

	for(int i=0;i<BIOPSY_SAMPLING_Z; i++)
	{
		regionsMinX[i] = new int[NUM_TARGETS[i]];
		regionsMaxX[i] = new int[NUM_TARGETS[i]];
		regionsMinZ[i] = new int[NUM_TARGETS[i]];
		regionsMaxZ[i] = new int[NUM_TARGETS[i]];

		for(int j=0; j<NUM_TARGETS[i]; j++)
		{
			regionsMinX[i][j] = minX + ((double)j/(double)NUM_TARGETS[i])*(maxX - minX);
			regionsMaxX[i][j] = minX + ((double)(j + 1)/(double)NUM_TARGETS[i])*(maxX - minX);
			regionsMinZ[i][j] = minZ + ((double)i/(double)BIOPSY_SAMPLING_Z)*(maxZ - minZ);
			regionsMaxZ[i][j] = maxZ + ((double)(i + 1)/(double)BIOPSY_SAMPLING_Z)*(maxZ - minZ);
		}
	}

	for(int z=0; z<BIOPSY_SAMPLING_Z; z++)
	{
		int zValue = minZ + ((double)(z + 1)/(double)(BIOPSY_SAMPLING_Z + 1))*(maxZ - minZ);
		int xValue = 0;

		for(int x=0; x<NUM_TARGETS[z]; x++)
		{
			if(RANDOM_TARGETING)
			{
				zValue = r.randomInt(regionsMinZ[z][x], regionsMaxZ[z][x]);
				xValue = r.randomInt(regionsMinX[z][x], regionsMaxX[z][x]);
			}
			else
			{
				xValue = minX + ((double)(x + 1)/(double)(NUM_TARGETS[z] + 1))*(maxX - minX);
			}

			if(xValue < statistics->getMinXForSlice(zValue - minZ))
			{
				xValue = r.randomInt(statistics->getMinXForSlice(zValue - minZ), regionsMaxX[z][x]);
			}
			else if(xValue > statistics->getMaxXForSlice(zValue - minZ))
			{
				xValue = r.randomInt(regionsMinX[z][x], statistics->getMaxXForSlice(zValue - minZ));
			}

			TargetPoint p_ITK;
			p_ITK.x = xValue;
			p_ITK.z = zValue;

			int posteriorYCoordinate = this->GetPosteriorYCoordinate(contour, statistics, xValue, zValue);

			if(TARGET_SURFACE)
			{
				p_ITK.y = posteriorYCoordinate;
			}
			else
			{
				int anteriorYCoordinate = this->GetAnteriorYCoordinate(contour, statistics, xValue, zValue);
				
				int minPeripheralZoneY = posteriorYCoordinate - (int)(((double)PERIPHERAL_ZONE_PERCENT/(double)100)*(posteriorYCoordinate - anteriorYCoordinate));
				int maxPeripheralZoneY = posteriorYCoordinate;

				if(RANDOM_TARGETING)
				{
					p_ITK.y = r.randomInt(minPeripheralZoneY, maxPeripheralZoneY);
				}
				else
				{
					//put y-coordinate in the middle of the peripheral zone
					p_ITK.y = minPeripheralZoneY + (int)(0.5*(maxPeripheralZoneY - minPeripheralZoneY));
				}
			}

			TargetPoint p_VTK = this->ConvertITKPointToVTKPoint(p_ITK);

			currentGroundTruthPreOpTargetLocationsITK.push_back(p_ITK);
			currentGroundTruthPreOpTargetLocationsVTK.push_back(p_VTK);
		}
	}

	this->GroundTruthPreOpTargetLocationsITK.push_back(currentGroundTruthPreOpTargetLocationsITK);
	this->GroundTruthPreOpTargetLocationsVTK.push_back(currentGroundTruthPreOpTargetLocationsVTK);

	//delete regions
	for(int i=0;i<BIOPSY_SAMPLING_Z; i++)
	{
		delete [] regionsMinX[i];
		delete [] regionsMaxX[i];
		delete [] regionsMinZ[i];
		delete [] regionsMaxZ[i];
	}
	delete [] regionsMinX;
	delete [] regionsMaxX;
	delete [] regionsMinZ;
	delete [] regionsMaxZ;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function determines the ground truth intra-op target locations for the contour at index "contourIndex". The ground-truth 
//intra-op target locations are determined by applying the ground-truth deformation field to the already placed target points.
void vtkKWProstateErrorMapRenderingWidget::DetermineGroundTruthIntraOpTargetLocations(int contourIndex)
{
	std::vector<TargetPoint> groundTruthPreOpTargetLocationsITK = this->GroundTruthPreOpTargetLocationsITK.at(contourIndex);
	DeformationFieldType::Pointer groundTruthDeformationField = this->GroundTruthDeformationFields.at(contourIndex);

	std::vector<TargetPoint> currentGroundTruthIntraOpTargetLocationsITK;
	std::vector<TargetPoint> currentGroundTruthIntraOpTargetLocationsVTK;
	OutputImageType::IndexType index;

	for(int p=0; p<groundTruthPreOpTargetLocationsITK.size(); p++)
	{
		TargetPoint preOpPoint_ITK = groundTruthPreOpTargetLocationsITK.at(p);

		//get translation vector for the pre-op target point
		index[0] = preOpPoint_ITK.x;
		index[1] = preOpPoint_ITK.y;
		index[2] = preOpPoint_ITK.z;

		VectorPixelType translationVector = groundTruthDeformationField->GetPixel(index);

		TargetPoint intraOpPoint_ITK;
		intraOpPoint_ITK.x = preOpPoint_ITK.x + translationVector[0];
		intraOpPoint_ITK.y = preOpPoint_ITK.y + translationVector[1];
		intraOpPoint_ITK.z = preOpPoint_ITK.z + translationVector[2];

		TargetPoint intraOpPoint_VTK = this->ConvertITKPointToVTKPoint(intraOpPoint_ITK);

		currentGroundTruthIntraOpTargetLocationsITK.push_back(intraOpPoint_ITK);
		currentGroundTruthIntraOpTargetLocationsVTK.push_back(intraOpPoint_VTK);
	}

	this->GroundTruthIntraOpTargetLocationsITK.push_back(currentGroundTruthIntraOpTargetLocationsITK);
	this->GroundTruthIntraOpTargetLocationsVTK.push_back(currentGroundTruthIntraOpTargetLocationsVTK);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function determines the resulting intra-op target locations for the contour at index "contourIndex" after the registration. The 
//intra-op target locations are determined by applying the resulting deformation field to the already placed target points.

void vtkKWProstateErrorMapRenderingWidget::DetermineResultTargetLocations(int contourIndex)
{
	std::vector<TargetPoint> groundTruthPreOpTargetLocationsITK = this->GroundTruthPreOpTargetLocationsITK.at(contourIndex);
	DeformationFieldType::Pointer registrationResultDeformationField = this->RegistrationResultDeformationFields.at(contourIndex);

	std::vector<TargetPoint> currentResultTargetLocationsITK;
	std::vector<TargetPoint> currentResultTargetLocationsVTK;
	OutputImageType::IndexType index;

	for(int p=0; p<groundTruthPreOpTargetLocationsITK.size(); p++)
	{
		TargetPoint preOpPoint_ITK = groundTruthPreOpTargetLocationsITK.at(p);

		//get translation vector for the pre-op target point
		index[0] = preOpPoint_ITK.x;
		index[1] = preOpPoint_ITK.y;
		index[2] = preOpPoint_ITK.z;

		VectorPixelType translationVector = registrationResultDeformationField->GetPixel(index);

		TargetPoint resultPoint_ITK;
		resultPoint_ITK.x = preOpPoint_ITK.x + translationVector[0];
		resultPoint_ITK.y = preOpPoint_ITK.y + translationVector[1];
		resultPoint_ITK.z = preOpPoint_ITK.z + translationVector[2];

		TargetPoint resultPoint_VTK =  this->ConvertITKPointToVTKPoint(resultPoint_ITK);
 
		currentResultTargetLocationsITK.push_back(resultPoint_ITK);
		currentResultTargetLocationsVTK.push_back(resultPoint_VTK);
	}

	this->ResultTargetLocationsITK.push_back(currentResultTargetLocationsITK);
	this->ResultTargetLocationsVTK.push_back(currentResultTargetLocationsVTK);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns the point on the mean shape corresponding to the 3D point "point".
//***This function has a bug. For the most part, the targeting error heat map has holes in it. This is due to the correspondence points
//***not being found correctly.
TargetPoint vtkKWProstateErrorMapRenderingWidget::GetMeanShapeCorrespondencePoint(TargetPoint point, ContourStatistics* statistics)
{
	TargetPoint meanShapeCorrespondencePoint;

	int minX = statistics->getMinX();
	int maxX = statistics->getMaxX();
	int minY = statistics->getMinX();
	int maxY = statistics->getMaxX();
	int minZ = statistics->getMinZ();
	int maxZ = statistics->getMaxZ();

	
	int minXSlice = statistics->getMinXForSlice(point.z - minZ);
	int maxXSlice = statistics->getMaxXForSlice(point.z - minZ);
	int minYSlice = statistics->getMinYForSlice(point.z - minZ);
	int maxYSlice = statistics->getMaxYForSlice(point.z - minZ);
	int minXLine = statistics->getMinX(point.y, point.z);
	int maxXLine = statistics->getMaxX(point.y, point.z);
	int minYLine = statistics->getMinY(point.x, point.z);
	int maxYLine = statistics->getMaxY(point.x, point.z);

	float percentX = (float)(point.x - minXLine)/(float)(maxXLine - minXLine);
	float percentY = (float)(point.y - minYLine)/(float)(maxYLine - minYLine);

	int meanShapeMinX = this->MeanShapeStatistics.getMinX();
	int meanShapeMaxX = this->MeanShapeStatistics.getMaxX();
	int meanShapeMinY = this->MeanShapeStatistics.getMinY();
	int meanShapeMaxY = this->MeanShapeStatistics.getMaxY();
	int meanShapeMinZ = this->MeanShapeStatistics.getMinZ();
	int meanShapeMaxZ = this->MeanShapeStatistics.getMaxZ();
	
	int meanShapeZ = meanShapeMinZ + point.z - minZ;

	int meanShapeMinXSlice = this->MeanShapeStatistics.getMinXForSlice(meanShapeZ - meanShapeMinZ);
	int meanShapeMaxXSlice = this->MeanShapeStatistics.getMaxXForSlice(meanShapeZ - meanShapeMinZ);
	int meanShapeMinYSlice = this->MeanShapeStatistics.getMinYForSlice(meanShapeZ - meanShapeMinZ);
	int meanShapeMaxYSlice = this->MeanShapeStatistics.getMaxYForSlice(meanShapeZ - meanShapeMinZ);

	int meanShapeX = meanShapeMinXSlice + ((float)(point.x - minXSlice)/(float)(maxXSlice - minXSlice))*(meanShapeMaxXSlice - meanShapeMinXSlice);
	int meanShapeY = meanShapeMinYSlice + ((float)(point.y - minYSlice)/(float)(maxYSlice - minYSlice))*(meanShapeMaxYSlice - meanShapeMinYSlice);

	if(meanShapeX < this->MeanShapeStatistics.getMinX(meanShapeY, meanShapeZ))
	{
		meanShapeX = this->MeanShapeStatistics.getMinX(meanShapeY, meanShapeZ);
	}
	else if(meanShapeX > this->MeanShapeStatistics.getMaxX(meanShapeY, meanShapeZ))
	{
		meanShapeX = this->MeanShapeStatistics.getMaxX(meanShapeY, meanShapeZ);
	}

	if(meanShapeY < this->MeanShapeStatistics.getMinY(meanShapeX, meanShapeZ))
	{
		meanShapeY = this->MeanShapeStatistics.getMinY(meanShapeX, meanShapeZ);
	}
	else if(meanShapeY > this->MeanShapeStatistics.getMaxY(meanShapeX, meanShapeZ))
	{
		meanShapeY = this->MeanShapeStatistics.getMaxY(meanShapeX, meanShapeZ);
	}

	meanShapeCorrespondencePoint.x = meanShapeX;
	meanShapeCorrespondencePoint.y = meanShapeY;
	meanShapeCorrespondencePoint.z = meanShapeZ;

	return meanShapeCorrespondencePoint;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function loops through the target points for the contour at index "contourIndex". It records the targeting error statistics associated
//with this contour.
void vtkKWProstateErrorMapRenderingWidget::RecordTargetingErrors(int contourIndex)
{
	double spacingX = this->MeanShapeDisplay->GetSpacing()[0];
	double spacingY = this->MeanShapeDisplay->GetSpacing()[1];
	double spacingZ = this->MeanShapeDisplay->GetSpacing()[2];
	int meanShapeMinX = this->MeanShapeStatistics.getMinX();
	int meanShapeMinY = this->MeanShapeStatistics.getMinY();
	int meanShapeMinZ = this->MeanShapeStatistics.getMinZ();

	ContourStatistics* currentContourStatistics = &this->GroundTruthContourStatistics.at(contourIndex);
	std::vector<TargetPoint> groundTruthPreOpTargetLocationsITK = this->GroundTruthPreOpTargetLocationsITK.at(contourIndex);
	std::vector<TargetPoint> groundTruthIntraOpTargetLocationsITK = this->GroundTruthIntraOpTargetLocationsITK.at(contourIndex);
	std::vector<TargetPoint> resultTargetLocationsITK = this->ResultTargetLocationsITK.at(contourIndex);

	std::vector<double> errors;

	for(int p=0; p<groundTruthPreOpTargetLocationsITK.size(); p++)
	{
		TargetPoint groundTruthPreOpPoint = groundTruthPreOpTargetLocationsITK.at(p);
		TargetPoint groundTruthIntraOpPoint = groundTruthIntraOpTargetLocationsITK.at(p);
		TargetPoint resultTargetPoint = resultTargetLocationsITK.at(p);

		TargetPoint meanShapePoint = this->GetMeanShapeCorrespondencePoint(groundTruthPreOpPoint, currentContourStatistics);

		TargetPoint differenceVoxelSpace;
		differenceVoxelSpace.x = resultTargetPoint.x - groundTruthIntraOpPoint.x;
		differenceVoxelSpace.y = resultTargetPoint.y - groundTruthIntraOpPoint.y;
		differenceVoxelSpace.z = resultTargetPoint.z - groundTruthIntraOpPoint.z;

		double differenceMetricSpaceX; //in mm
		double differenceMetricSpaceY; //in mm
		double differenceMetricSpaceZ; //in mm
		differenceMetricSpaceX = differenceVoxelSpace.x * spacingX;
		differenceMetricSpaceY = differenceVoxelSpace.y * spacingY;
		differenceMetricSpaceZ = differenceVoxelSpace.z * spacingZ;

		double euclideanDistanceToGroundTruth = sqrt(differenceMetricSpaceX*differenceMetricSpaceX + differenceMetricSpaceY*differenceMetricSpaceY + differenceMetricSpaceZ*differenceMetricSpaceZ);

		int xInDataStructure = meanShapePoint.x - meanShapeMinX;
		int yInDataStructure = meanShapePoint.y - meanShapeMinY;
		int zInDataStructure = meanShapePoint.z - meanShapeMinZ;

		this->TargetingErrorStructure[xInDataStructure][yInDataStructure][zInDataStructure].targeted = true;
		this->TargetingErrorStructure[xInDataStructure][yInDataStructure][zInDataStructure].errors.push_back(euclideanDistanceToGroundTruth);
		errors.push_back(euclideanDistanceToGroundTruth);
	}

	//record targeting error statistics for this test case
	double totalError = 0;
	double minError = 99999;
	double maxError = 0;
	int n = errors.size();
	int totalInsignificant = 0;
	int totalSignificant = 0;

	//
	//find min, max, mean
	//
	for(int e=0; e<n; e++)
	{
		double error = errors.at(e);

		if(error < minError)
		{
			minError = error;
		}

		if(error > maxError)
		{
			maxError = error;
		}

		if(error <= MAX_INSIGNIFICANT_ERROR)
		{
			totalInsignificant++;
		}

		if(error >= MIN_SIGNIFICANT_ERROR)
		{
			totalSignificant++;
		}

		totalError = totalError + error;
	}

	double meanError = (double)totalError/(double)n;

	//
	//find variance, standard deviation, and mean deviation
	//

	double varianceSum = 0;
	double meanDeviationSum = 0;

	for(int e=0; e<n; e++)
	{
		double error = errors.at(e);
		varianceSum = varianceSum + (error - meanError)*(error - meanError);
		meanDeviationSum = meanDeviationSum + abs(error - meanError);
	}

	double variance = (double)varianceSum/(double)(n - 1);
	double standardDeviation = sqrt(variance);
	double meanDeviation = (double)meanDeviationSum/(double)n;
	double percentInsignificant = 100*((double)totalInsignificant/(double)n);
	double percentSignificant = 100*((double)totalSignificant/(double)n);

	ReportStruct report;
	report.totalTargets = n;
	report.numInsignificant = totalInsignificant;
	report.numSignificant = totalSignificant;
	report.percentInsignificant = percentInsignificant;
	report.percentSignificant = percentSignificant;
	report.minError = minError;
	report.maxError = maxError;
	report.meanError = meanError;
	report.meanDeviation = meanDeviation;
	report.standardDeviation = standardDeviation;

	this->TargetingErrorReports.push_back(report);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function records the targeting error statistics with respect to the entire simulation. i.e. it averages out all the errors at each
//voxel on the mean shape. It stores the errors in the "TargetingErrorStructure"
void vtkKWProstateErrorMapRenderingWidget::RecordFinalTargetingErrorStatistics()
{
	int totalX = this->MeanShapeStatistics.getMaxX() - this->MeanShapeStatistics.getMinX() + 1;
	int totalY = this->MeanShapeStatistics.getMaxY() - this->MeanShapeStatistics.getMinY() + 1;
	int totalZ = this->MeanShapeStatistics.getMaxZ() - this->MeanShapeStatistics.getMinZ() + 1;

	for(int x=0; x<totalX; x++)
	{
		for(int y=0; y<totalY; y++)
		{
			for(int z=0; z<totalZ; z++)
			{
				if(this->TargetingErrorStructure[x][y][z].targeted)
				{
					double totalError = 0;
					double minError = 99999;
					double maxError = 0;
					int n = this->TargetingErrorStructure[x][y][z].errors.size();

					//
					//find min, max, mean
					//
					for(int e=0; e<n; e++)
					{
						double error = this->TargetingErrorStructure[x][y][z].errors.at(e);

						if(error < minError)
						{
							minError = error;
						}

						if(error > maxError)
						{
							maxError = error;
						}

						totalError = totalError + error;
					}

					double meanError = (double)totalError/(double)n;

					//
					//find variance, population variance, standard deviation, population standard deviation
					//

					double varianceSum = 0;

					for(int e=0; e<n; e++)
					{
						double error = this->TargetingErrorStructure[x][y][z].errors.at(e);
						varianceSum = varianceSum + (error - meanError)*(error - meanError);
					}

					double variance = (double)varianceSum/(double)(n - 1);
					
					this->TargetingErrorStructure[x][y][z].minError = minError;
					this->TargetingErrorStructure[x][y][z].maxError = maxError;
					this->TargetingErrorStructure[x][y][z].meanError = meanError;
					this->TargetingErrorStructure[x][y][z].standardDeviation = sqrt(variance);
				}
			}
		}
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function obtains the report for the entire simulation. The report contains targeting errror statistics for the entire simulation.
//The report is obtained from the results in the "TargetingErrorStructure"
ReportStruct vtkKWProstateErrorMapRenderingWidget::GetFinalReport()
{
	ReportStruct finalReport;

	int totalX = this->MeanShapeStatistics.getMaxX() - this->MeanShapeStatistics.getMinX() + 1;
	int totalY = this->MeanShapeStatistics.getMaxY() - this->MeanShapeStatistics.getMinY() + 1;
	int totalZ = this->MeanShapeStatistics.getMaxZ() - this->MeanShapeStatistics.getMinZ() + 1;

	int totalTargets = 0;
	double globalMinError = 99999;
	double globalMaxError = 0;
	double totalMean = 0;
	double totalStandardDeviation = 0;
	int totalInsignificant = 0;
	int totalSignificant = 0;

	for(int x=0; x<totalX; x++)
	{
		for(int y=0; y<totalY; y++)
		{
			for(int z=0; z<totalZ; z++)
			{
				if(this->TargetingErrorStructure[x][y][z].targeted)
				{
					totalTargets++;

					//update global min/max
					if(this->TargetingErrorStructure[x][y][z].minError < globalMinError)
					{
						globalMinError = this->TargetingErrorStructure[x][y][z].minError;
					}
					if(this->TargetingErrorStructure[x][y][z].maxError > globalMaxError)
					{
						globalMaxError = this->TargetingErrorStructure[x][y][z].maxError;
					}

					//update total insignificant/significant
					if(this->TargetingErrorStructure[x][y][z].meanError <= MAX_INSIGNIFICANT_ERROR)
					{
						totalInsignificant++;
					}
					else if(this->TargetingErrorStructure[x][y][z].meanError >= MIN_SIGNIFICANT_ERROR)
					{
						totalSignificant++;
					}

					totalMean = totalMean + this->TargetingErrorStructure[x][y][z].meanError;
				}
			}
		}
	}

	double meanOfMeanErrors = (double)totalMean/(double)totalTargets;
	double percentInsignificant = 100*((double)totalInsignificant/(double)totalTargets);
	double percentSignificant = 100*((double)totalSignificant/(double)totalTargets);


	//determine mean deviation and standard deviation
	double varianceSum = 0;
	double meanDeviationSum = 0;

	for(int x=0; x<totalX; x++)
	{
		for(int y=0; y<totalY; y++)
		{
			for(int z=0; z<totalZ; z++)
			{
				if(this->TargetingErrorStructure[x][y][z].targeted)
				{
					varianceSum = varianceSum + (this->TargetingErrorStructure[x][y][z].meanError - meanOfMeanErrors)*(this->TargetingErrorStructure[x][y][z].meanError - meanOfMeanErrors);
					meanDeviationSum = meanDeviationSum + abs(this->TargetingErrorStructure[x][y][z].meanError - meanOfMeanErrors);
				}
			}
		}
	}

	double meanStandardDeviation = sqrt((double)varianceSum/(double)(totalTargets - 1));
	double globalMeanDeviation = (double)meanDeviationSum/(double)totalTargets;

	finalReport.totalTargets = totalTargets;
	finalReport.numInsignificant = totalInsignificant;
	finalReport.numSignificant = totalSignificant;
	finalReport.percentInsignificant = percentInsignificant;
	finalReport.percentSignificant = percentSignificant;
	finalReport.minError = globalMinError;
	finalReport.maxError = globalMaxError;
	finalReport.meanError = meanOfMeanErrors;
	finalReport.meanDeviation = globalMeanDeviation;
	finalReport.standardDeviation = meanStandardDeviation;

	return finalReport;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function applies segmentation error to the ITK image "contour" according to the method in the paper
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::CreateContourWithSegmentationError(OutputImageType::Pointer contour, ContourStatistics& statistics)
{
	int zRange = statistics.getMaxZ() - statistics.getMinZ();
	int halfProstateLength = zRange/2;
	int midZ = statistics.getMinZ() + halfProstateLength;
	int apexMinZ = statistics.getMaxZ() - (double)PERCENT_APEX/(double)100 * (zRange);
	int baseMaxZ = statistics.getMinZ() + (double)PERCENT_BASE/(double)100 * (zRange);

	double spacingX = contour->GetSpacing()[0];
	double spacingY = contour->GetSpacing()[1];

	//create a new empty contour, with the same size as the original
	OutputImageType::Pointer segmentedContour = OutputImageType::New();

	//get the size of the original image
	OutputImageType::SizeType originalSize = contour->GetLargestPossibleRegion().GetSize();

	//set the size to be the same
	OutputImageType::SizeType newSize;
	newSize[0] = originalSize[0];
	newSize[1] = originalSize[1];
	newSize[2] = originalSize[2];

	//get the origin of the original image
	OutputImageType::PointType origin = contour->GetOrigin();

	OutputImageType::IndexType start;
	start[0] = origin[0];
	start[1] = origin[1];
	start[2] = origin[2];

	//set up the region
	OutputImageType::RegionType region;
	region.SetSize(newSize);
	region.SetIndex(start);

	segmentedContour->SetRegions(region);
	segmentedContour->Allocate();

	//set spacing
	segmentedContour->SetSpacing(contour->GetSpacing());

	//set origin
	segmentedContour->SetOrigin(contour->GetOrigin());

	typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
	IteratorType iterator2(segmentedContour, segmentedContour->GetLargestPossibleRegion());

	//initialize all pixels to zero
	for(iterator2.GoToBegin(); !iterator2.IsAtEnd(); ++iterator2)
	{
		OutputImageType::IndexType index = iterator2.GetIndex();
		OutputImageType::PixelType pixel = iterator2.Get();

		segmentedContour->SetPixel(index, 0);
	}

	//introduce segmentation error
	Random r;

	//NEW WAY
	for(int i=0; i<statistics.getNumSlices(); i++)
	{
		OutputImageType::IndexType index;

		int z = i + statistics.getMinZ();
		int distanceToMid = abs(midZ - z);

		std::vector<float> contourErrorPointsR;
		std::vector<float> contourErrorPointsTheta;

		Point gravityCenterPoint;
		gravityCenterPoint.x = statistics.getGravityXForSlice(i);
		gravityCenterPoint.y = statistics.getGravityYForSlice(i);

		//gather points with segmentation error at each segment
		float angleStep = (float)2*PI/(float)SEGMENTATION_ERROR_SEGMENTS;
		float angle = -angleStep; //start at zero

		bool currentSign; //true if positive, false if negative
		int currentSignPoints = 0;
		float lastSegError;

		for(int s=0; s<SEGMENTATION_ERROR_SEGMENTS; s++)
		{
			angle = angle + angleStep;

			Point extremalPoint = this->GetSliceExtremalPoint(contour, gravityCenterPoint.x, gravityCenterPoint.y, z, angle);

			float segError = 0;

			float minError = 0;
			float maxError = 0;

			if((s == 0) || (currentSignPoints >= MIN_FLUCTUATION_POINTS))
			{
				if(z < midZ) //slice is closer to the base
				{
					minError = MIN_MID_SEGMENTATION_ERROR + (((float)distanceToMid/(float)halfProstateLength)*abs(MIN_BASE_SEGMENTATION_ERROR - MIN_MID_SEGMENTATION_ERROR));
					maxError = MAX_MID_SEGMENTATION_ERROR + (((float)distanceToMid/(float)halfProstateLength)*abs(MAX_BASE_SEGMENTATION_ERROR - MAX_MID_SEGMENTATION_ERROR));
				}
				else //slice is closer to the apex
				{
					minError = MIN_MID_SEGMENTATION_ERROR + (((float)distanceToMid/(float)halfProstateLength)*abs(MIN_APEX_SEGMENTATION_ERROR - MIN_MID_SEGMENTATION_ERROR));
					maxError = MAX_MID_SEGMENTATION_ERROR + (((float)distanceToMid/(float)halfProstateLength)*abs(MAX_APEX_SEGMENTATION_ERROR - MAX_MID_SEGMENTATION_ERROR));
				}
			}
			else
			{
				if(currentSign)
				{
					if((lastSegError - MAX_SEGMENTATION_ERROR_DIFFERENCE) < 0)
					{
						minError = 0;
					}
					else
					{
						minError = lastSegError - MAX_SEGMENTATION_ERROR_DIFFERENCE;
					}
					maxError = lastSegError;
				}
				else
				{
					minError = lastSegError;

					if(lastSegError + MAX_SEGMENTATION_ERROR_DIFFERENCE > 0)
					{
						maxError = 0;
					}
					else
					{
						maxError = lastSegError + MAX_SEGMENTATION_ERROR_DIFFERENCE;
					}	
				}
			}


			segError = r.randomFloat(minError, maxError);
			
			if((s == 0) || (currentSignPoints >= MIN_FLUCTUATION_POINTS))
			{
				if(r.flipCoin(0.5))
				{
					segError = -segError;
				}
			}

			lastSegError = segError;

			currentSign = (segError >= 0);

			if(currentSignPoints >= MIN_FLUCTUATION_POINTS)
			{
				currentSignPoints = 0;
			}

			currentSignPoints++;

			//stretch or compress the contour
			float normalizedX = 0;
			float normalizedY = 0;

			float difX = extremalPoint.x - gravityCenterPoint.x;
			float difY = extremalPoint.y - gravityCenterPoint.y;

			if((difX != 0) || (difY != 0))
			{
				normalizedX = (double)difX/(double)sqrt((double)(difX*difX + difY*difY));
				normalizedY = (double)difY/(double)sqrt((double)(difX*difX + difY*difY));
			}

			int xSegError = (int)((double)segError*normalizedX/(double)spacingX);
			int ySegError = (int)((double)segError*normalizedY/(double)spacingY);

			Point extremalPointSegError;
			extremalPointSegError.x = extremalPoint.x + xSegError;
			extremalPointSegError.y = extremalPoint.y + ySegError;

			PointPolar extremalPointSegErrorInPolarCoordinates = this->ConvertToPolarCoordinates(extremalPointSegError, gravityCenterPoint);

			contourErrorPointsR.push_back(extremalPointSegErrorInPolarCoordinates.r);
			contourErrorPointsTheta.push_back(extremalPointSegErrorInPolarCoordinates.theta);
		}

		contourErrorPointsR.push_back(contourErrorPointsR.at(0));
		contourErrorPointsTheta.push_back(2*PI);

		CubicSplineSolver* cubicSplineSolver = new CubicSplineSolver(contourErrorPointsTheta, contourErrorPointsR);

		//fill in new contour
		float step = 0.001;
		int steps = (float)360/(float)step;

		for(int s=0; s<steps; s++)
		{
			PointPolar contourPointInPolarCoordinates;

			float angleRadians = s*step*TO_RADIANS;
			contourPointInPolarCoordinates.theta = angleRadians;
			contourPointInPolarCoordinates.r = cubicSplineSolver->evaluate(angleRadians);

			Point contourPoint = this->ConvertFromPolarCoordinates(contourPointInPolarCoordinates, gravityCenterPoint);

			index[0] = contourPoint.x;
			index[1] = contourPoint.y;
			index[2] = z;

			segmentedContour->SetPixel(index, 1);
		}

		//fill in rest of slice through region growing
		this->FillSlice(segmentedContour, gravityCenterPoint.x, gravityCenterPoint.y, z);

		delete cubicSplineSolver;
	}

	return segmentedContour;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function deforms the ITK contour "image" with the deformation field "deformationField"
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::DeformImage(OutputImageType::Pointer image, DeformationFieldType::Pointer deformationField)
{
	//create the filter that applies the deformation
	typedef itk::WarpImageFilter<OutputImageType, OutputImageType, DeformationFieldType>  WarpImageFilterType;
	WarpImageFilterType::Pointer filter = WarpImageFilterType::New();

	//use linear interpolation method for float-valued outputs of the deformation
	typedef itk::LinearInterpolateImageFunction<OutputImageType, double>  InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	filter->SetInterpolator(interpolator);
	filter->SetOutputSpacing(deformationField->GetSpacing());
	filter->SetOutputOrigin(deformationField->GetOrigin());
	filter->SetDeformationField(deformationField);

	filter->SetInput(image);

	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		exit(-1);
	}

	return filter->GetOutput();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function deforms the ITK contour "contour" according to the method in the paper
OutputImageType::Pointer vtkKWProstateErrorMapRenderingWidget::CreateContourWithDeformationField(OutputImageType::Pointer contour, ContourStatistics& statistics)
{
	int numSlices = statistics.getNumSlices();
	int minX = statistics.getMinX();
	int maxX = statistics.getMaxX();
	int minY = statistics.getMinY();
	int maxY = statistics.getMaxY();
	int minZ = statistics.getMinZ();
	int maxZ = statistics.getMaxZ();

	double spacingX = contour->GetSpacing()[0];
	double spacingY = contour->GetSpacing()[1];

	bool ***setDeformationVector;
	setDeformationVector = new bool**[maxX - minX + 1];
	for(int i=0; i<(maxX - minX + 1); i++)
	{
		setDeformationVector[i] = new bool*[maxY - minY + 1];

		for(int j=0; j<(maxY - minY + 1); j++)
		{
			setDeformationVector[i][j] = new bool[maxZ - minZ + 1];

			for(int k=0; k<(maxZ - minZ + 1); k++)
			{
				setDeformationVector[i][j][k] = false;
			}
		}
	}

	//create the deformation field
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();

	//get the size of the input contour
	OutputImageType::SizeType inputImageSize = contour->GetLargestPossibleRegion().GetSize();

	//set the size of the deformation field to be the same as the original contour
	DeformationFieldType::SizeType deformationFieldSize;
	deformationFieldSize[0] = inputImageSize[0];
	deformationFieldSize[1] = inputImageSize[1];
	deformationFieldSize[2] = inputImageSize[2];

	//get the origin of the input contour
	OutputImageType::PointType origin = contour->GetOrigin();

	DeformationFieldType::IndexType start;
	start[0] = origin[0];
	start[1] = origin[1];
	start[2] = origin[2];

	//set up the region
	DeformationFieldType::RegionType region;
	region.SetSize(deformationFieldSize);
	region.SetIndex(start);

	deformationField->SetRegions(region);
	deformationField->Allocate();

	//set spacing
	deformationField->SetSpacing(contour->GetSpacing());

	//set origin
	deformationField->SetOrigin(contour->GetOrigin());

	//random number generator
	Random r;

	//set up iterator to loop through the image
	typedef itk::ImageRegionIteratorWithIndex<OutputImageType> IteratorType;
	IteratorType iterator(contour, contour->GetLargestPossibleRegion());

	//iterate through the image, if there's a zero (i.e. no prostate), set the deformation
	//vector to zero, otherwise assign a non-zero deformation vector
	for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
	{
		OutputImageType::IndexType index = iterator.GetIndex();
		OutputImageType::PixelType pixel = iterator.Get();

		VectorPixelType deformationVector;

		DeformationFieldType::IndexType deformationFieldIndex;
		deformationFieldIndex[0] = index[0];
		deformationFieldIndex[1] = index[1];
		deformationFieldIndex[2] = index[2];

		if(((int)pixel) == 1)
		{
			//have deformation applied in-plane at the edges of the contour 
			//(in order for segmentation error simulation to work)

			if(this->IsExtremalContourVoxel(contour, index, minX, maxX, minY, maxY, minZ, maxZ))
			{
				//deform this contour voxel with probability: DEFORMATION_PROBABILITY
				if(r.flipCoin(DEFORMATION_PROBABILITY))
				{
					int x = index[0];
					int y = index[1];
					int z = index[2];

					int difX = x - statistics.getGravityXForSlice(z - minZ);
					int difY = y - statistics.getGravityYForSlice(z - minZ);

					double deformationMagnitude = r.randomInt(MIN_DEFORMATION_MAGNITUDE, MAX_DEFORMATION_MAGNITUDE);

					if(r.flipCoin(0.5))
					{
						deformationMagnitude = deformationMagnitude*(-1);
					}

					double normalizedX = 0;
					double normalizedY = 0;

					if((difX != 0) || (difY != 0))
					{
						normalizedX = (double)difX/(double)sqrt((double)(difX*difX + difY*difY));
						normalizedY = (double)difY/(double)sqrt((double)(difX*difX + difY*difY));
					}

					//double xDeformationReal = ((double)deformationMagnitude*normalizedX/(double)spacingX);
					//double yDeformationReal = ((double)deformationMagnitude*normalizedY/(double)spacingY);
					int xDeformationVoxelized = (int)((double)deformationMagnitude*normalizedX/(double)spacingX);
					int yDeformationVoxelized = (int)((double)deformationMagnitude*normalizedY/(double)spacingY);

					//stretch (or compress) the contour
					Point startPoint;
					startPoint.x = x;
					startPoint.y = y;

					Point endPoint;
					endPoint.x = x + xDeformationVoxelized;
					endPoint.y = y + yDeformationVoxelized;

					//int xLength = abs(endPoint.x - startPoint.x);
					//int yLength = abs(endPoint.y - startPoint.y);

					std::vector<Point> pointsToDeform = this->GetDiscretizedPath(startPoint, endPoint);

					for(int i=0; i<pointsToDeform.size(); i++)
					{
						Point pointToDeform = pointsToDeform.at(i);

						deformationVector[0] = pointToDeform.x - endPoint.x;
						deformationVector[1] = pointToDeform.y - endPoint.y;
						deformationVector[2] = 0;

						deformationFieldIndex[0] = pointToDeform.x;
						deformationFieldIndex[1] = pointToDeform.y;
						deformationFieldIndex[2] = index[2];

						if(deformationMagnitude < 0)
						{
							setDeformationVector[pointToDeform.x - minX][pointToDeform.y - minY][index[2] - minZ] = true;
						}

						deformationField->SetPixel(deformationFieldIndex, deformationVector);
					}
				}
				else
				{
					deformationVector[0] = 0;
					deformationVector[1] = 0;
					deformationVector[2] = 0;

					deformationField->SetPixel(deformationFieldIndex, deformationVector);
				}
			}
			else if(!setDeformationVector[index[0] - minX][index[1] - minY][index[2] - minZ])
			{
				deformationVector[0] = 0;
				deformationVector[1] = 0;
				deformationVector[2] = 0;

				deformationField->SetPixel(deformationFieldIndex, deformationVector);
			}

		}
		else //no prostate here, no deformation field
		{
			deformationVector[0] = 0;
			deformationVector[1] = 0;
			deformationVector[2] = 0;

			deformationField->SetPixel(deformationFieldIndex, deformationVector);
		}
	}

	this->GroundTruthDeformationFields.push_back(deformationField);

	for(int i=0; i<(maxX - minX + 1); i++)
	{
		for(int j=0; j<(maxY - minY + 1); j++)
		{
			delete [] setDeformationVector[i][j];
		}

		delete [] setDeformationVector[i];
	}

	delete [] setDeformationVector;

	return this->DeformImage(contour, deformationField);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function solves the registration using DEMONS between the moving image "movingImage" and fixed image "fixedImage" and stores the
//result in a deformation field.
//NOTE: This function may not work...It doesn't give good results for the input contours
DeformationFieldType::Pointer vtkKWProstateErrorMapRenderingWidget::SolveDeformableRegistrationDemons(OutputImageType::Pointer movingImage, OutputImageType::Pointer fixedImage)
{
	typedef itk::DemonsRegistrationFilter<OutputImageType,OutputImageType,DeformationFieldType> DeformableRegistrationFilterType;
	DeformableRegistrationFilterType::Pointer deformableRegistrationFilter = DeformableRegistrationFilterType::New();

	deformableRegistrationFilter->SetFixedImage(fixedImage);
	deformableRegistrationFilter->SetMovingImage(movingImage);
	deformableRegistrationFilter->SetNumberOfIterations(150);

	try
	{
		deformableRegistrationFilter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		exit(-1);
	}

	return deformableRegistrationFilter->GetDeformationField();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function solves the registration using LEVEL SETS between the moving image "movingImage" and fixed image "fixedImage" and stores the
//result in a deformation field.

DeformationFieldType::Pointer vtkKWProstateErrorMapRenderingWidget::SolveDeformableRegistrationLevelSets(OutputImageType::Pointer movingImage, OutputImageType::Pointer fixedImage)
{
	typedef itk::CastImageFilter<OutputImageType, InternalImageType > FixedImageCasterType;
	typedef itk::CastImageFilter<OutputImageType, InternalImageType > MovingImageCasterType;

	FixedImageCasterType::Pointer fixedImageCaster   = FixedImageCasterType::New();
	MovingImageCasterType::Pointer movingImageCaster = MovingImageCasterType::New();

	fixedImageCaster->SetInput(fixedImage);
	movingImageCaster->SetInput(movingImage); 

	typedef itk::Vector<float, Dimension>    VectorPixelType;
	typedef itk::LevelSetMotionRegistrationFilter<InternalImageType,InternalImageType,DeformationFieldType>   RegistrationFilterType;
	RegistrationFilterType::Pointer filter = RegistrationFilterType::New();


	filter->SetFixedImage(fixedImageCaster->GetOutput());
	filter->SetMovingImage(movingImageCaster->GetOutput());


	filter->SetNumberOfIterations(50);
	filter->SetGradientSmoothingStandardDeviations(4);


	filter->Update();

	return filter->GetDeformationField();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function solves the registration using B SPLINES between the moving image "movingImage" and fixed image "fixedImage" and stores the
//result in a deformation field.
//NOTE: This function may not work...It doesn't give good results for the input contours
DeformationFieldType::Pointer vtkKWProstateErrorMapRenderingWidget::SolveDeformableRegistrationBSplines(OutputImageType::Pointer movingImage, OutputImageType::Pointer fixedImage)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef itk::CastImageFilter<OutputImageType, InternalImageType > FixedImageCasterType;
	typedef itk::CastImageFilter<OutputImageType, InternalImageType > MovingImageCasterType;

	FixedImageCasterType::Pointer fixedImageCaster   = FixedImageCasterType::New();
	MovingImageCasterType::Pointer movingImageCaster = MovingImageCasterType::New();

	fixedImageCaster->SetInput(fixedImage);
	movingImageCaster->SetInput(movingImage); 

	fixedImageCaster->Update();
	movingImageCaster->Update();

	InternalImageType::Pointer internalFixedImage = fixedImageCaster->GetOutput();
	InternalImageType::Pointer internalMovingImage = movingImageCaster->GetOutput();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef itk::BSplineDeformableTransform<CoordinateRepType, Dimension, SplineOrder> TransformType;
	typedef itk::LBFGSOptimizer  OptimizerType;
	typedef itk::MeanSquaresImageToImageMetric<InternalImageType, InternalImageType >  MetricType;
	typedef itk::LinearInterpolateImageFunction<InternalImageType, double>   InterpolatorType;
	typedef itk::ImageRegistrationMethod<InternalImageType, InternalImageType >   RegistrationType;

	MetricType::Pointer         metric        = MetricType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetInterpolator(interpolator);

	////////////////////////////////////////////////////////////////////////

	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform(transform);

	registration->SetFixedImage(internalFixedImage);
	registration->SetMovingImage(internalMovingImage);

	InternalImageType::RegionType fixedRegion = internalFixedImage->GetBufferedRegion();

	registration->SetFixedImageRegion(fixedRegion);

	////////////////////////////////////////////////////////////////////////

	typedef TransformType::RegionType RegionType;
	RegionType bsplineRegion;
	RegionType::SizeType   gridSizeOnImage;
	RegionType::SizeType   gridBorderSize;
	RegionType::SizeType   totalGridSize;

	gridSizeOnImage.Fill(5);
	gridBorderSize.Fill(3);    // Border for spline order = 3 ( 1 lower, 2 upper )
	totalGridSize = gridSizeOnImage + gridBorderSize;

	bsplineRegion.SetSize(totalGridSize);

	typedef TransformType::SpacingType SpacingType;
	SpacingType spacing = internalFixedImage->GetSpacing();

	typedef TransformType::OriginType OriginType;
	OriginType origin = internalFixedImage->GetOrigin();;

	OutputImageType::SizeType fixedImageSize = fixedRegion.GetSize();

	for(unsigned int r=0; r<Dimension; r++)
	{
		spacing[r] *= static_cast<double>(fixedImageSize[r] - 1)  / 
					  static_cast<double>(gridSizeOnImage[r] - 1);
		origin[r]  -=  spacing[r]; 
	}

	transform->SetGridSpacing(spacing);
	transform->SetGridOrigin(origin);
	transform->SetGridRegion(bsplineRegion);
	transform->SetGridDirection(internalFixedImage->GetDirection());

	////////////////////////////////////////////////////////////////////////

	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters = transform->GetNumberOfParameters();

	ParametersType parameters( numberOfParameters );

	parameters.Fill( 0.0 );

	transform->SetParameters( parameters );

	registration->SetInitialTransformParameters( transform->GetParameters() );

	////////////////////////////////////////////////////////////////////////

	optimizer->SetGradientConvergenceTolerance( 0.05 );
	optimizer->SetLineSearchAccuracy( 0.9 );
	optimizer->SetDefaultStepLength( 1.5 );
	optimizer->TraceOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations( 1000 );

	try 
    { 
		registration->StartRegistration();
    } 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cout << "ExceptionObject caught !" << std::endl; 
		std::cout << err << std::endl; 
		exit(-1);
	}

	OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

	transform->SetParameters( finalParameters );

	////////////////////////////////////////////////////////////////////////

    DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
    deformationField->SetRegions(fixedRegion);
    deformationField->SetOrigin(internalFixedImage->GetOrigin() );
    deformationField->SetSpacing(internalFixedImage->GetSpacing());
    deformationField->Allocate();

    typedef itk::ImageRegionIterator<DeformationFieldType> FieldIterator;
    FieldIterator fi(deformationField, fixedRegion);

    fi.GoToBegin();

    TransformType::InputPointType  fixedPoint;
    TransformType::OutputPointType movingPoint;
    DeformationFieldType::IndexType index;

    VectorPixelType displacement;

    while(!fi.IsAtEnd())
	{
		index = fi.GetIndex();
		deformationField->TransformIndexToPhysicalPoint(index, fixedPoint);
		movingPoint = transform->TransformPoint(fixedPoint);
		displacement[0] = movingPoint[0] - fixedPoint[0];
		displacement[1] = movingPoint[1] - fixedPoint[1];
		displacement[2] = movingPoint[2] - fixedPoint[2];
		fi.Set(displacement);
		++fi;
	}

	return deformationField;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns true if any target point in the vector "targetPoints" is within the ITK contour "contour" (i.e. if any target point is
//NOT on the surface of the contour)
bool vtkKWProstateErrorMapRenderingWidget::HasInternalTargetPoints(OutputImageType::Pointer contour, ContourStatistics* statistics, std::vector<TargetPoint>& targetPoints)
{
	OutputImageType::IndexType index;

	for(int p=0; p<targetPoints.size(); p++)
	{
		TargetPoint targetPoint = targetPoints.at(p);

		index[0] = targetPoint.x;
		index[1] = targetPoint.y;
		index[2] = targetPoint.z;

		if(!this->IsSurfaceVoxel(index, contour, statistics))
		{
			return true;
		}
	}

	return false;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function displays the VTK image stored in "VTKImage" in the render widget "aRenderWidget"
void vtkKWProstateErrorMapRenderingWidget::DisplayVTKImage(vtkImageData* VTKImage, vtkKWRenderWidget** aRenderWidget)
{
	//create a grayscale lookup table
	vtkLookupTable *lookupTable = vtkLookupTable::New();
	lookupTable->SetRange(0, 255); // image intensity range
	lookupTable->SetValueRange(0.0, 1.0); // from black to white
	lookupTable->SetSaturationRange(0.0, 0.0); // no color saturation
	lookupTable->SetRampToLinear();
	lookupTable->Build();

	//Render the Prostate Volume
	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();

	vtkMarchingCubes* iso = vtkMarchingCubes::New();
	iso->SetInput(VTKImage);
	iso->SetValue(0, 1);

	mapper->SetInput(iso->GetOutput());
	mapper->SetLookupTable(lookupTable);

	(*aRenderWidget)->GetRenderer()->RemoveAllViewProps();

	vtkActor* prostateContourActor = vtkActor::New();
	prostateContourActor->SetMapper(mapper);
	(*aRenderWidget)->GetRenderer()->AddActor(prostateContourActor);

	lookupTable->Delete();
	mapper->Delete();
	iso->Delete();
	prostateContourActor->Delete();

	(*aRenderWidget)->Reset();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function displays the VTK image stored in "VTKImage" in the render widget "aRenderWidget" along with any target points
void vtkKWProstateErrorMapRenderingWidget::DisplayVTKImage(vtkImageData* VTKImage, vtkKWRenderWidget** aRenderWidget, OutputImageType::Pointer contour, ContourStatistics* statistics, std::vector<TargetPoint>& targetPointsITK, std::vector<TargetPoint>& targetPointsVTK)
{
	//create a grayscale lookup table
	vtkLookupTable *lookupTable = vtkLookupTable::New();
	lookupTable->SetRange(0, 255); // image intensity range
	lookupTable->SetValueRange(0.0, 1.0); // from black to white
	lookupTable->SetSaturationRange(0.0, 0.0); // no color saturation
	lookupTable->SetRampToLinear();
	lookupTable->Build();

	//Render the Prostate Volume
	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();

	vtkMarchingCubes* iso = vtkMarchingCubes::New();
	iso->SetInput(VTKImage);
	iso->SetValue(0, 1);

	mapper->SetInput(iso->GetOutput());
	mapper->SetLookupTable(lookupTable);

	(*aRenderWidget)->GetRenderer()->RemoveAllViewProps();

	vtkActor* prostateContourActor = vtkActor::New();
	prostateContourActor->SetMapper(mapper);

	if(this->HasInternalTargetPoints(contour, statistics, targetPointsITK))
	{
		prostateContourActor->GetProperty()->SetOpacity(0.3);
	}

	(*aRenderWidget)->GetRenderer()->AddActor(prostateContourActor);

	//create sphere actors for each target point and add them to the render widget
	if(targetPointsVTK.size() > MAX_TARGET_POINTS_TO_RENDER)
	{
		// get MAX_TARGET_POINTS_TO_RENDER random points from the points to show

		std::vector<TargetPoint> targetPointsCopyVTK = targetPointsVTK;
		Random r;
		for(int p=0; p<MAX_TARGET_POINTS_TO_RENDER; p++)
		{
			int randomIndex = r.randomInt(0, targetPointsCopyVTK.size() - 1);
			TargetPoint targetPoint = targetPointsCopyVTK.at(randomIndex);

			vtkSphereSource* sphereSource = vtkSphereSource::New();
			vtkPolyDataMapper* sphereMapper = vtkPolyDataMapper::New();
			vtkActor* sphereActor = vtkActor::New();

			sphereSource->SetCenter(targetPoint.x, targetPoint.y, targetPoint.z);
			sphereSource->SetRadius(2.0);

			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

			sphereActor->SetMapper(sphereMapper);
			sphereActor->GetProperty()->SetColor(0, 0, 1);
			
			(*aRenderWidget)->GetRenderer()->AddActor(sphereActor);

			sphereSource->Delete();
			sphereMapper->Delete();
			sphereActor->Delete();

			//delete the target point
			targetPointsCopyVTK.erase(targetPointsCopyVTK.begin()+randomIndex);
		}
	}
	else
	{
		for(int p=0; p<targetPointsVTK.size(); p++)
		{
			TargetPoint targetPoint = targetPointsVTK.at(p);

			vtkSphereSource* sphereSource = vtkSphereSource::New();
			vtkPolyDataMapper* sphereMapper = vtkPolyDataMapper::New();
			vtkActor* sphereActor = vtkActor::New();

			sphereSource->SetCenter(targetPoint.x, targetPoint.y, targetPoint.z);
			sphereSource->SetRadius(2.0);

			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

			sphereActor->SetMapper(sphereMapper);
			sphereActor->GetProperty()->SetColor(0, 0, 1);
			
			(*aRenderWidget)->GetRenderer()->AddActor(sphereActor);

			sphereSource->Delete();
			sphereMapper->Delete();
			sphereActor->Delete();
		}
	}

	lookupTable->Delete();
	mapper->Delete();
	iso->Delete();
	prostateContourActor->Delete();

	(*aRenderWidget)->Reset();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function generates the mean shape, modes of variation, and eigenvalues, and stores them in the directory 
//"C:\ProstateErrorMapRendering\Output"
//The user is prompted to select a directory with at least one aligned shape stored as a ".mha" file. To select
//the directory, the user must select one of the ".mha" files stored in the directory

void vtkKWProstateErrorMapRenderingWidget::GenerateMeanShape()
{
	vtkKWFileBrowserDialog* dialog = vtkKWFileBrowserDialog::New();
	dialog->SetParent(this);
	dialog->Create();
	dialog->Invoke();

	if(dialog->GetStatus() == vtkKWDialog::StatusOK)
	{
		//get chosen directory
		std::string fileNamePath(dialog->GetFileName());
		int indexOfLastSlash = fileNamePath.find_last_of("/");
		std::string directory = fileNamePath.substr(0, indexOfLastSlash + 1);

		//replace forward slashes with backslashes
		bool done = false;
		while(!done)
		{
			int location = directory.find("/");
			done = (location == std::string::npos);

			if(!done)
			{
				directory.replace(location,1,"\\");
			}
		}
		
		//create directory for output
		std::string outputDirectory = "C:\\ProstateErrorMapRendering\\Output";
		CreateDirectory(outputDirectory.c_str(), NULL);

		//find all ".mha" files in the directory
		bool working = true;
		int counter = 0;
		std::string buffer;
		std::vector<std::string> fileNames;

		std::string searchDirectory = directory;
		searchDirectory.append("*");

		WIN32_FIND_DATA findData;
		HANDLE handle = FindFirstFile(searchDirectory.c_str(),&findData);

		if(handle != INVALID_HANDLE_VALUE)
		{
			buffer = findData.cFileName;

			if((buffer.size() > 4) && (buffer.substr(buffer.size() - 4, buffer.size()) == ".mha"))
			{
				fileNames.push_back(directory + buffer);
			}

			while(working)
			{
				  FindNextFile(handle,&findData);
				  if(findData.cFileName != buffer)
				  {
						 buffer = findData.cFileName;

						 if((buffer.size() > 4) && (buffer.substr(buffer.size() - 4, buffer.size()) == ".mha"))
						 {
							fileNames.push_back(directory + buffer);
						 }
				  }
				  else
				  {
						  //end of files reached
						  working = false;
				  }
			}
		}

		std::cout<<"reading aligned shapes"<<std::endl<<std::endl;

		//store aligned shapes and create mean shape

		//vector/array of the distance images that will be generated
		std::vector<InternalImageType::Pointer> distanceImages(fileNames.size());

		//image reader
		typedef itk::ImageFileReader< InternalImageType  > ImageReaderType;
		ImageReaderType::Pointer  ImageReader  = ImageReaderType::New();

		//writer
		typedef itk::ImageFileWriter< InternalImageType >  WriterType;

		typedef itk::ShiftScaleImageFilter<InternalImageType, InternalImageType> DivideFilterType;

		//define the distance map generator filter
		typedef itk::SignedDanielssonDistanceMapImageFilter<InternalImageType, InternalImageType> DistanceMapFilterType;

		//define the PCA Shape Model estimator
		typedef itk::ImagePCAShapeModelEstimator<InternalImageType, InternalImageType> ShapeEstimatorType;
		ShapeEstimatorType::Pointer shapeEstimator = ShapeEstimatorType::New();
		shapeEstimator->SetNumberOfTrainingImages(fileNames.size());
		shapeEstimator->SetNumberOfPrincipalComponentsRequired(3);

		// note that the filter has to process each training image individually
		// whereas the PCA shape estimator will process all the training maps generated
		// simultaneously

		//in a for loop read all the training images one by one
		//and generate distance map for each
		for (int i=0; i<fileNames.size(); i++)
		{	
			std::cout<<"reading aligned shape: "<<fileNames.at(i)<<std::endl;

			ImageReader->SetFileName(fileNames.at(i).c_str());
			DistanceMapFilterType::Pointer filterDaniel = DistanceMapFilterType::New();
			filterDaniel->SetInput(ImageReader->GetOutput());
			filterDaniel->Update();
			distanceImages[i] = filterDaniel->GetOutput();	
			shapeEstimator->SetInput(i, distanceImages[i]);

			std::cout<<"done reading aligned shape"<<std::endl<<std::endl;
		}

		try
		{
		  shapeEstimator->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
			exit(-1);
		}

		std::cout<<std::endl;

		std::cout<<"outputting signed distance map of mean shape"<<std::endl;
		std::cout<<std::endl;

		//output distance map of mean shape
		this->MeanShape = InternalImageType::New();
		this->MeanShape = shapeEstimator->GetOutput(0);

		WriterType::Pointer      writer =  WriterType::New();
		DivideFilterType::Pointer divider = DivideFilterType::New();

		std::string meanShapeDistanceMapFile = outputDirectory;
		meanShapeDistanceMapFile.append("\\MeanShape_DMap.mha");

		writer->SetFileName(meanShapeDistanceMapFile.c_str());
		writer->SetInput(this->MeanShape);
		writer->Update();

		std::cout<<"mean shape signed distance map generation complete"<<std::endl<<std::endl;


		std::cout<<"generating contour of mean shape"<<std::endl;

		//this is still signed distance map, need a contour

		vnl_vector<double> eigenValues(3); 
		eigenValues = shapeEstimator->GetEigenValues();

		typedef itk::UnaryFunctorImageFilter<
					InternalImageType,OutputImageType,ContourOfSignedDistanceMapFunctor<typename InternalImageType::PixelType> > ContourExtractorType;

		ContourExtractorType::Pointer contourExtractor = ContourExtractorType::New();

		contourExtractor->SetInput(this->MeanShape);

		contourExtractor->Update();

		OutputImageType::Pointer MeanImageContour = contourExtractor->GetOutput();

		typedef itk::ImageFileWriter< OutputImageType >  BinaryWriterType;
		BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();

		std::string meanShapeFile = outputDirectory;
		meanShapeFile.append("\\MeanShape.mha");
		binaryWriter->SetFileName(meanShapeFile.c_str());
		binaryWriter->SetInput(MeanImageContour);
		binaryWriter->Update();

		std::cout<<"contour generation complete"<<std::endl<<std::endl;

		std::cout<<"generating modes of variation"<<std::endl<<std::endl;

		char *Suffix = new char[20]; 

		//store modes of variation
		for (int mode=1; mode<=3; mode++)
		{	
			//form the file name
			sprintf(Suffix, "%02d.mha", mode);
			std::string variationModeFile = outputDirectory;
			variationModeFile.append("\\VariationMode_");
			variationModeFile.append(Suffix);

			this->EigenModes.push_back(shapeEstimator->GetOutput(mode));

			writer->SetFileName(variationModeFile.c_str());
			writer->SetInput(this->EigenModes[mode-1]);
			writer->Update();
		}
		delete [] Suffix;

		std::cout<<"generation of modes of variation complete"<<std::endl<<std::endl;

		//store eigen values
		std::cout<<"storing eigenvalues"<<std::endl<<std::endl;
		
		std::string eigenValueFilePath = outputDirectory;
		eigenValueFilePath.append("\\eigenValues.txt");

		std::ofstream eigenValueFile(eigenValueFilePath.c_str());

		if (eigenValueFile.is_open())
		{
			for (int mode=1; mode<=3; mode++)
			{
				eigenValueFile<<eigenValues[mode-1]<<"\n";
			}

			eigenValueFile.close();
		}
		else
		{
			std::cerr << "Fatal Error: Could not store eigenvalues" << std::endl;
			exit(-1);
		}
		std::cout<<"finished storing eigenvalues"<<std::endl<<std::endl;

		//convert mean shape to VTK and display
		this->DisplayContour(meanShapeFile);
	}

	dialog->Delete();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function loads a prostate contour (selected by the user) and displays it
void vtkKWProstateErrorMapRenderingWidget::LoadProstateContour()
{
	vtkKWFileBrowserDialog* dialog = vtkKWFileBrowserDialog::New();
	dialog->SetParent(this);
	dialog->Create();
	dialog->Invoke();

	if(dialog->GetStatus() == vtkKWDialog::StatusOK)
	{
		this->DisplayContour(dialog->GetFileName());
	}

	dialog->Delete();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function returns true if the point at index "index" is along the surface of the ITK contour "contour"
bool vtkKWProstateErrorMapRenderingWidget::IsSurfaceVoxel(OutputImageType::IndexType index, OutputImageType::Pointer contour, ContourStatistics* statistics)
{
	if(((int)contour->GetPixel(index)) == 0)
	{
		return false;
	}

	if((index[2] == statistics->getMinZ()) || (index[2] == statistics->getMaxZ()))
	{
		return true;
	}

	OutputImageType::IndexType checkIndex;

	//right
	checkIndex[0] = index[0] + 1;
	checkIndex[1] = index[1];
	checkIndex[2] = index[2];

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//left
	checkIndex[0] = index[0] - 1;

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//up
	checkIndex[0] = index[0];
	checkIndex[1] = index[1] + 1;

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//down
	checkIndex[1] = index[1] - 1;

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//front
	checkIndex[1] = index[1];
	checkIndex[2] = index[2] + 1;

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	//back
	checkIndex[2] = index[2] - 1;

	if(((int)contour->GetPixel(checkIndex)) == 0)
	{
		return true;
	}

	return false;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function displays the contour stored at the absolute file "filePath" in the main render widget
void vtkKWProstateErrorMapRenderingWidget::DisplayContour(std::string filePath)
{
	//Get the ITK volume
	OutputImageType::Pointer ITKVolume = this->GetITKImageFromFile(filePath);

	//convert to VTK
	vtkImageData* toDisplay = this->GetVTKImageData(ITKVolume, true);

	//display the volume
	this->DisplayVTKImage(toDisplay, &this->RenderWidget);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function displays the mean shape (along with its heat map) in the main render widget
void vtkKWProstateErrorMapRenderingWidget::DisplayMeanShapeWithHeatMap()
{
	OutputImageType::IndexType index;

	int minX = this->MeanShapeStatistics.getMinX();
	int maxX = this->MeanShapeStatistics.getMaxX();
	int minY = this->MeanShapeStatistics.getMinY();
	int maxY = this->MeanShapeStatistics.getMaxY();
	int minZ = this->MeanShapeStatistics.getMinZ();
	int maxZ = this->MeanShapeStatistics.getMaxZ();

	double spacing[3];
	for(int i=0; i<3; i++)
	{
		spacing[i] = this->MeanShapeDisplay->GetSpacing()[i];
	}

	double minSpacing = spacing[0];
	for(int i=1; i<3; i++)
	{
		if(spacing[i] < minSpacing)
		{
			minSpacing = spacing[i];
		}
	}

	//determine spacing ratio for determining length of voxels
	double spacingRatio[3];
	for(int i=0; i<3; i++)
	{
		spacingRatio[i] = (double)spacing[i]/(double)minSpacing;
	}

	double xLength = spacingRatio[0];
	double yLength = spacingRatio[1];
	double zLength = spacingRatio[2];

	//draw the contour, with a cube at each voxel
	for(int x=minX; x<=maxX; x++)
	{
		for(int y=minY; y<=maxY; y++)
		{
			for(int z=minZ; z<=maxZ; z++)
			{
				index[0] = x;
				index[1] = y;
				index[2] = z;

				if(this->IsSurfaceVoxel(index, this->MeanShapeDisplay, &this->MeanShapeStatistics))
				{
					//draw cube at this voxel
					vtkCubeSource* cubeSource = vtkCubeSource::New();
					vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
					vtkActor* cubeActor = vtkActor::New();

					cubeSource->SetCenter((x - minX)*xLength + (double)xLength/(double)2, (maxY - (y - minY))*yLength + (double)yLength/(double)2, (maxZ - (z - minZ))*zLength + (double)zLength/(double)2);
					cubeSource->SetXLength(xLength);
					cubeSource->SetYLength(yLength);
					cubeSource->SetZLength(zLength);

					cubeMapper->SetInputConnection(cubeSource->GetOutputPort());

					cubeActor->SetMapper(cubeMapper);

					double red = 1.0;
					double green = 1.0;
					double blue = 1.0; //default color is white

					bool targeted = false;

					if(DO_REGISTRATION)
					{
						//make heat map

						if(this->TargetingErrorStructure[x - minX][y - minY][z - minZ].targeted)
						{
							targeted = true;

							//get the mean error for this targeted voxel
							//and determine its color on the heat map
							double meanError = this->TargetingErrorStructure[x - minX][y - minY][z - minZ].meanError;

							if(meanError <= MAX_INSIGNIFICANT_ERROR)
							{
								//mean error is insignificant, colour is green
								red = 0;
								green = 1;
								blue = 0;
							}
							else if(meanError >= MIN_SIGNIFICANT_ERROR)
							{
								//mean error is significant, colour is red
								red = 1;
								green = 0;
								blue = 0;
							}
							else
							{
								//mean error is somewhere in between
								ColorHSV colorHSV;
								colorHSV.hue = 120 - ((double)(meanError - MAX_INSIGNIFICANT_ERROR)/(double)(MIN_SIGNIFICANT_ERROR - MAX_INSIGNIFICANT_ERROR))*120;
								colorHSV.saturation = 1.0;
								colorHSV.value = 1.0;

								ColorRGB colorRGB = ConvertFromHSVToRGB(colorHSV);
								red = colorRGB.red;
								green = colorRGB.green;
								blue = colorRGB.blue;
							}
						}
					}

					cubeActor->GetProperty()->SetColor(red, green, blue);
					if(DO_REGISTRATION && !targeted)
					{
						cubeActor->GetProperty()->SetOpacity(0.1);
					}

					this->RenderWidget->GetRenderer()->AddActor(cubeActor);
					
					cubeSource->Delete();
					cubeMapper->Delete();
					cubeActor->Delete();
				}
			}
		}
	}

	this->RenderWidget->Reset();
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function displays one case of the simulation in the GUI (depending on the value of the scale). If the scale value is zero, the mean
//shape is shown, otherwise the simulation case associated with that scale value is shown
void vtkKWProstateErrorMapRenderingWidget::DisplaySimulationCase(int scaleValue)
{
	if(scaleValue == 0)
	{
		if(this->multipleImageFramesShowing)
		{
			this->ResetGUIToOneImageFrame();
		}

		if(DO_REGISTRATION)
		{
			this->DisplayMeanShapeWithHeatMap();
		}
		else
		{
			this->DisplayVTKImage(this->GetVTKImageData(this->MeanShapeDisplay, true), &this->RenderWidget);
		}
	}
	else
	{
		if(!this->multipleImageFramesShowing)
		{
			this->SetupGUIWithMultipleImageFrames(scaleValue);
		}

		this->UpdateLabels(scaleValue);
		
		this->DisplayVTKImage(this->GetVTKImageData(this->ITKContours.at(scaleValue - 1), true), &this->RenderWidget, this->ITKContours.at(scaleValue - 1), &this->GroundTruthContourStatistics.at(scaleValue - 1), this->GroundTruthPreOpTargetLocationsITK.at(scaleValue - 1), this->GroundTruthPreOpTargetLocationsVTK.at(scaleValue - 1));
		this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursSegmentationError.at(scaleValue - 1), true), &this->SegmentationErrorRenderWidget);
		this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursSegmentationErrorDifferences.at(scaleValue - 1), true), &this->SegmentationErrorDifferenceRenderWidget);
		
		if(DO_DEFORMATION)
		{
			this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursWithDeformation.at(scaleValue - 1), true), &this->DeformationRenderWidget, this->ITKContoursWithDeformation.at(scaleValue - 1), &this->DeformationContourStatistics.at(scaleValue - 1), this->GroundTruthIntraOpTargetLocationsITK.at(scaleValue - 1), this->GroundTruthIntraOpTargetLocationsVTK.at(scaleValue - 1));
			this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursWithDeformationAndSegmentationError.at(scaleValue - 1), true), &this->DeformationSegmentationErrorRenderWidget);
			this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursWithDeformationAndSegmentationErrorDifferences.at(scaleValue - 1), true), &this->DeformationSegmentationErrorDifferenceRenderWidget);
		}

		if(DO_REGISTRATION)
		{
			this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursRegistrationResults.at(scaleValue - 1), true), &this->RegistrationResultRenderWidget, this->ITKContoursRegistrationResults.at(scaleValue - 1), &this->RegistrationResultsContourStatistics.at(scaleValue - 1), this->ResultTargetLocationsITK.at(scaleValue - 1), this->ResultTargetLocationsVTK.at(scaleValue - 1));
			this->DisplayVTKImage(this->GetVTKImageData(this->ITKContoursRegistrationDifferences.at(scaleValue - 1), true), &this->RegistrationDifferenceRenderWidget);
		}
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function updates the text in the display according to the simulation case
void vtkKWProstateErrorMapRenderingWidget::UpdateLabels(int scaleValue)
{
	//update labels for simulation info
	SimulationInfoStruct currentSimulationInfo = this->SimulationInfo.at(scaleValue - 1);

	std::string modeFactorStrings[3];

	for(int mode=1; mode<=3; mode++)
	{
		modeFactorStrings[mode-1] = "Weight mode ";
		modeFactorStrings[mode-1].append(this->GetStringFromInteger(mode));
		modeFactorStrings[mode-1].append(": ");

		if(currentSimulationInfo.sigmaFactor[mode-1] != 0)
		{
			if(currentSimulationInfo.sigmaFactor[mode-1] > 0)//is positive
			{
				modeFactorStrings[mode-1].append("+");
			}
			else//is negative
			{
				modeFactorStrings[mode-1].append("-");
			}
		}

		modeFactorStrings[mode-1].append(this->GetStringFromDouble(abs(currentSimulationInfo.sigmaFactor[mode-1]), 2));

		if(currentSimulationInfo.sigmaFactor[mode-1] != 0)
		{
			modeFactorStrings[mode-1].append(" Sigma");
		}
	}

	this->FirstModeLabel->SetText(modeFactorStrings[0].c_str());
	this->SecondModeLabel->SetText(modeFactorStrings[1].c_str());
	this->ThirdModeLabel->SetText(modeFactorStrings[2].c_str());

	this->BaseSegmentationErrorLabel->SetText(baseSegmentationErrorString.c_str());
	this->MidSegmentationErrorLabel->SetText(midSegmentationErrorString.c_str());
	this->ApexSegmentationErrorLabel->SetText(apexSegmentationErrorString.c_str());

	if(DO_REGISTRATION)
	{
		//update labels for targeting error report
		std::string totalTargetsString = "Total Targets: ";
		totalTargetsString.append(this->GetStringFromInteger(this->TargetingErrorReports.at(scaleValue - 1).totalTargets));
		this->TotalTargetsLabel->SetText(totalTargetsString.c_str());

		std::string insignificantString = "Total Insignificant ";
		insignificantString.append("(<= ");
		insignificantString.append(this->GetStringFromDouble(MAX_INSIGNIFICANT_ERROR, 1));
		insignificantString.append(" mm): ");
		insignificantString.append(this->GetStringFromInteger(this->TargetingErrorReports.at(scaleValue - 1).numInsignificant));
		insignificantString.append(" (");
		insignificantString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).percentInsignificant, 1));
		insignificantString.append("%)");
		this->InsignificantLabel->SetText(insignificantString.c_str());

		std::string significantString = "Total Significant ";
		significantString.append("(>= ");
		significantString.append(this->GetStringFromDouble(MIN_SIGNIFICANT_ERROR, 1));
		significantString.append(" mm): ");
		significantString.append(this->GetStringFromInteger(this->TargetingErrorReports.at(scaleValue - 1).numSignificant));
		significantString.append(" (");
		significantString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).percentSignificant, 1));
		significantString.append("%)");
		this->SignificantLabel->SetText(significantString.c_str());

		std::string minErrorString = "Min Error: ";
		minErrorString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).minError, 1));
		minErrorString.append(" mm");
		this->MinErrorLabel->SetText(minErrorString.c_str());

		std::string maxErrorString = "Max Error: ";
		maxErrorString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).maxError, 1));
		maxErrorString.append(" mm");
		this->MaxErrorLabel->SetText(maxErrorString.c_str());

		std::string meanErrorString = "Mean Error: ";
		meanErrorString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).meanError, 1));
		meanErrorString.append(" mm");
		this->MeanErrorLabel->SetText(meanErrorString.c_str());

		std::string meanDeviationString = "Mean Deviation: ";
		meanDeviationString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).meanDeviation, 1));
		meanDeviationString.append(" mm");
		this->MeanDeviationLabel->SetText(meanDeviationString.c_str());

		std::string standardDeviationString = "Standard Deviation: ";
		standardDeviationString.append(this->GetStringFromDouble(this->TargetingErrorReports.at(scaleValue - 1).standardDeviation, 2));
		this->StandardDeviationLabel->SetText(standardDeviationString.c_str());
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function allocates the targeting error data structure
void vtkKWProstateErrorMapRenderingWidget::SetupTargetingErrorDataStructure()
{
	int totalX = this->MeanShapeStatistics.getMaxX() - this->MeanShapeStatistics.getMinX() + 1;
	int totalY = this->MeanShapeStatistics.getMaxY() - this->MeanShapeStatistics.getMinY() + 1;
	int totalZ = this->MeanShapeStatistics.getMaxZ() - this->MeanShapeStatistics.getMinZ() + 1;

	this->TargetingErrorStructure = new VoxelTargetingErrorStruct**[totalX];

	for(int x=0; x<totalX; x++)
	{
		this->TargetingErrorStructure[x] = new VoxelTargetingErrorStruct*[totalY];

		for(int y=0; y<totalY; y++)
		{
			this->TargetingErrorStructure[x][y] = new VoxelTargetingErrorStruct[totalZ];

			for(int z=0; z<totalZ; z++)
			{
				this->TargetingErrorStructure[x][y][z].targeted = false;
			}
		}
	}
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function destroys the targeting error data structure
void vtkKWProstateErrorMapRenderingWidget::DestroyTargetingErrorDataStructure()
{
	int totalX = this->MeanShapeStatistics.getMaxX() - this->MeanShapeStatistics.getMinX() + 1;
	int totalY = this->MeanShapeStatistics.getMaxY() - this->MeanShapeStatistics.getMinY() + 1;
	int totalZ = this->MeanShapeStatistics.getMaxZ() - this->MeanShapeStatistics.getMinZ() + 1;

	for(int x=0; x<totalX; x++)
	{
		for(int y=0; y<totalY; y++)
		{
			delete [] this->TargetingErrorStructure[x][y];
		}

		delete [] this->TargetingErrorStructure[x];
	}

	delete [] this->TargetingErrorStructure;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function runs the simulation with the workflow as per the paper
void vtkKWProstateErrorMapRenderingWidget::RunSimulation()
{
	int timeBegin = time(NULL);

	typedef itk::ImageFileReader<InternalImageType> ImageReaderType;

	typedef itk::ImageFileWriter<InternalImageType> WriterType;
    WriterType::Pointer  writer =  WriterType::New();

	typedef itk::ImageFileWriter<OutputImageType>  BinaryWriterType;
    BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();

	typedef itk::ShiftScaleImageFilter<InternalImageType, InternalImageType> DivideFilterType;
    DivideFilterType::Pointer divider = DivideFilterType::New();

	typedef itk::AddImageFilter<InternalImageType,InternalImageType,InternalImageType> AddFilterType;
    AddFilterType::Pointer addFilter = AddFilterType::New();

    typedef itk::SubtractImageFilter<InternalImageType, InternalImageType, InternalImageType> SubtractFilterType;
    SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();

	typedef itk::UnaryFunctorImageFilter<
				InternalImageType,OutputImageType,ContourOfSignedDistanceMapFunctor<typename InternalImageType::PixelType> > ContourExtractorType;
	

	std::string outputDirectory = "C:\\ProstateErrorMapRendering\\Output";

	//create directory for simulation cases
	std::string randomContourDirectory = "C:\\ProstateErrorMapRendering\\SimulationCases";
	CreateDirectory(randomContourDirectory.c_str(), NULL);

	//delete all previously created random contours
	bool working = true;
	std::string buffer;
	WIN32_FIND_DATA findData;
	HANDLE handle = FindFirstFile("C:\\ProstateErrorMapRendering\\SimulationCases\\*",&findData);

	if(handle != INVALID_HANDLE_VALUE)
	{
		buffer = findData.cFileName;

		if(buffer.size() >= 4)
		{
			std::string file = randomContourDirectory;
			file.append("\\");
			file.append(buffer);
			DeleteFile(file.c_str());
		}

		while(working)
		{
			  FindNextFile(handle,&findData);
			  if(findData.cFileName != buffer)
			  {
					 buffer = findData.cFileName;

					 if(buffer.size() >= 4)
					 {
						 std::string file = randomContourDirectory;
						 file.append("\\");
						 file.append(buffer);
						 DeleteFile(file.c_str());
					 }
			  }
			  else
			  {
					  //end of files reached
					  working = false;
			  }
		}
	}

	//get eigenvalues
	vnl_vector<double> eigenValues(3);

	std::string eigenValueFile = outputDirectory;
	eigenValueFile.append("\\eigenValues.txt");

	char* line = new char[2000];
	std::ifstream eigenValueInputStream(eigenValueFile.c_str());

	eigenValueInputStream.getline(line, 2000);
	eigenValues[0] = strtod(line, NULL);

	eigenValueInputStream.getline(line, 2000);
	eigenValues[1] = strtod(line, NULL);

	eigenValueInputStream.getline(line, 2000);
	eigenValues[2] = strtod(line, NULL);

	eigenValueInputStream.close();

	//get EigenModes

	char* Suffix = new char[20];

	if(this->EigenModes.size() == 0)
	{
		for(int mode=1; mode<=3; mode++)
		{
			//form the file name
			sprintf(Suffix, "%02d.mha", mode);
			std::string variationModeFile = outputDirectory;
			variationModeFile.append("\\VariationMode_");
			variationModeFile.append(Suffix);

			ImageReaderType::Pointer  ImageReader  = ImageReaderType::New();
			ImageReader->SetFileName(variationModeFile.c_str());

			try
			{
				ImageReader->Update();
			}
			catch( itk::ExceptionObject & err ) 
			{ 
				std::cerr << "ExceptionObject caught !" << std::endl; 
				std::cerr << err << std::endl; 
				exit(-1);
			}

			this->EigenModes.push_back(ImageReader->GetOutput());
		}
	}

	//get mean shape
	if(!this->MeanShape)
	{
		ImageReaderType::Pointer  ImageReader  = ImageReaderType::New();

		std::string meanShapeFile = outputDirectory;
		meanShapeFile.append("\\MeanShape_DMap.mha");
		ImageReader->SetFileName(meanShapeFile.c_str());

		try
		{
			ImageReader->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
			exit(-1);
		}

		this->MeanShape = InternalImageType::New();
		this->MeanShape = ImageReader->GetOutput();
	}

	this->MeanShapeDisplay = this->GetITKImageFromFile("C:\\ProstateErrorMapRendering\\Output\\MeanShape.mha");
	this->GatherContourStatistics(this->MeanShapeDisplay, &this->MeanShapeStatistics);

	if(DO_REGISTRATION)
	{
		this->SetupTargetingErrorDataStructure();
	}

	//generate n random prostate contours
	Random random;

	std::cout<<"Generating "<<SIMULATION_CASES<<" simulation case";
	if(SIMULATION_CASES > 1)
	{
		std::cout<<"s";
	}
	std::cout<<std::endl<<std::endl;

	this->SimulationInfo.clear();
	this->ITKContours.clear();
	this->ITKContoursSegmentationError.clear();
	this->ITKContoursSegmentationErrorDifferences.clear();
	this->ITKContoursWithDeformation.clear();
	this->ITKContoursWithDeformationAndSegmentationError.clear();
	this->ITKContoursWithDeformationAndSegmentationErrorDifferences.clear();
	this->RegistrationResultDeformationFields.clear();
	this->ITKContoursRegistrationResults.clear();

	for(int i=1; i<=SIMULATION_CASES; i++)
	{
		int timeBeginCurrentCase = time(NULL);

		if(SIMULATION_CASES > 1)
		{
		std::cout<<"*************************************"<<std::endl;
		std::cout<<"Generating simulation case "<<i<<std::endl<<std::endl;
		}

		ContourExtractorType::Pointer contourExtractor = ContourExtractorType::New();

		std::cout<<"Generating pre-op prostate contour"<<std::endl;

		//randomize PCA weights
		SimulationInfoStruct currentSimulationInfo;

		currentSimulationInfo.sigmaFactor[0] = random.randomDouble(MIN_WEIGHT_MODE_1, MAX_WEIGHT_MODE_1);
		currentSimulationInfo.sigmaFactor[1] = random.randomDouble(MIN_WEIGHT_MODE_2, MAX_WEIGHT_MODE_2);
		currentSimulationInfo.sigmaFactor[2] = random.randomDouble(MIN_WEIGHT_MODE_3, MAX_WEIGHT_MODE_3);

		std::cout<<"Weight Mode 1: "<<currentSimulationInfo.sigmaFactor[0]<<std::endl;
		std::cout<<"Weight Mode 2: "<<currentSimulationInfo.sigmaFactor[1]<<std::endl;
		std::cout<<"Weight Mode 3: "<<currentSimulationInfo.sigmaFactor[2]<<std::endl<<std::endl;

		//force them to have two decimal places
		currentSimulationInfo.sigmaFactor[0] = (double)((int)(currentSimulationInfo.sigmaFactor[0] * 100))/(double)100;
		currentSimulationInfo.sigmaFactor[1] = (double)((int)(currentSimulationInfo.sigmaFactor[1] * 100))/(double)100;
		currentSimulationInfo.sigmaFactor[2] = (double)((int)(currentSimulationInfo.sigmaFactor[2] * 100))/(double)100;

		this->SimulationInfo.push_back(currentSimulationInfo);

		std::string contourFile = "C:\\ProstateErrorMapRendering\\SimulationCases\\Case_";
		contourFile.append(this->GetStringFromInteger(i));
		contourFile.append(".mha");

		std::string distanceMapFile = "C:\\ProstateErrorMapRendering\\SimulationCases\\Case_";
		distanceMapFile.append(this->GetStringFromInteger(i));
		distanceMapFile.append("_DMap.mha");

		InternalImageType::Pointer image = this->MeanShape;

		for(int mode = 1; mode <= 3; mode++)
		{
			if(currentSimulationInfo.sigmaFactor[mode - 1] > 0) //is positive
			{
				addFilter->SetInput1(image);
				divider->SetInput(this->EigenModes[mode - 1]);
				divider->SetScale(currentSimulationInfo.sigmaFactor[mode - 1]*sqrt(eigenValues[mode - 1]));
				divider->Update();
				addFilter->SetInput2(divider->GetOutput());
				addFilter->Update();
				image = addFilter->GetOutput();
			}
			else//is negative
			{
				subtractFilter->SetInput1(image);
				divider->SetInput(this->EigenModes[mode - 1]);
				divider->SetScale(currentSimulationInfo.sigmaFactor[mode - 1]*sqrt(eigenValues[mode - 1]));
				divider->Update();
				subtractFilter->SetInput2(divider->GetOutput());
				subtractFilter->Update();
				image = subtractFilter->GetOutput();
			}
		}

		writer->SetFileName(distanceMapFile.c_str());
		writer->SetInput(image);
		writer->Update();

		contourExtractor->SetInput(image);
		contourExtractor->Update();

		binaryWriter->SetFileName(contourFile.c_str());
		binaryWriter->SetInput(contourExtractor->GetOutput());
		binaryWriter->Update();

		int baseSegmentationError = 0;
		int midSegmentationError = 0;
		int apexSegmentationError = 0;
		int deformationBaseSegmentationError = 0;
		int deformationMidSegmentationError = 0;
		int deformationApexSegmentationError = 0;

		this->ITKContours.push_back(contourExtractor->GetOutput());

		//grow contour to be in the same coordinate system as the segmented/deformed contours
		this->GrowInitialContour(this->ITKContours.at(i-1));

		ContourStatistics currentGroundTruthContourStatistics;
		this->GatherContourStatistics(this->ITKContours.at(i-1), &currentGroundTruthContourStatistics);
		this->GroundTruthContourStatistics.push_back(currentGroundTruthContourStatistics);

		if(TARGET_ALL && !RANDOM_TARGETING)
		{
			this->PlaceGroundTruthPreOpTargetLocationsAll(i-1);
		}
		else if(TARGET_ALL && RANDOM_TARGETING)
		{
			this->PlaceGroundTruthPreOpTargetLocationsAllRandom(i-1);
		}
		else
		{
			this->PlaceGroundTruthPreOpTargetLocations(i-1);
		}

		std::cout<<"Applying segmentation error to pre-op prostate"<<std::endl;
		this->ITKContoursSegmentationError.push_back(this->CreateContourWithSegmentationError(this->ITKContours.at(i-1), this->GroundTruthContourStatistics.at(i-1)));
		this->ITKContoursSegmentationErrorDifferences.push_back(this->Subtract(this->ITKContoursSegmentationError.at(i-1), this->ITKContours.at(i-1)));

		if(DO_DEFORMATION)
		{
			std::cout<<"Generating intra-op prostate through deformation"<<std::endl;
			this->ITKContoursWithDeformation.push_back(this->CreateContourWithDeformationField(this->ITKContours.at(i-1), this->GroundTruthContourStatistics.at(i-1)));

			ContourStatistics currentDeformationContourStatistics;
			this->GatherContourStatistics(this->ITKContoursWithDeformation.at(i-1), &currentDeformationContourStatistics);
			this->DeformationContourStatistics.push_back(currentDeformationContourStatistics);

			this->DetermineGroundTruthIntraOpTargetLocations(i-1);

			std::cout<<"Applying segmentation error to intra-op prostate"<<std::endl;
			this->ITKContoursWithDeformationAndSegmentationError.push_back(this->CreateContourWithSegmentationError(this->ITKContoursWithDeformation.at(i-1), this->DeformationContourStatistics.at(i-1)));
			this->ITKContoursWithDeformationAndSegmentationErrorDifferences.push_back(this->Subtract(this->ITKContoursWithDeformationAndSegmentationError.at(i-1), this->ITKContoursWithDeformation.at(i-1)));
		}

		if(DO_REGISTRATION)
		{
			if(USING_DEMONS_REGISTRATION)
			{
				std::cout<<"Performing deformable registration: DEMONS"<<std::endl<<std::endl;
				this->RegistrationResultDeformationFields.push_back(this->SolveDeformableRegistrationDemons(this->ITKContoursSegmentationError.at(i-1), this->ITKContoursWithDeformationAndSegmentationError.at(i-1)));
			}
			else if(USING_BSPLINES_REGISTRATION)
			{
				std::cout<<"Performing deformable registration: B-Splines"<<std::endl<<std::endl;
				this->RegistrationResultDeformationFields.push_back(this->SolveDeformableRegistrationBSplines(this->ITKContoursSegmentationError.at(i-1), this->ITKContoursWithDeformationAndSegmentationError.at(i-1)));
			}
			else //using levels sets
			{
				std::cout<<"Performing deformable registration: Level Sets"<<std::endl<<std::endl;
				this->RegistrationResultDeformationFields.push_back(this->SolveDeformableRegistrationLevelSets(this->ITKContoursSegmentationError.at(i-1), this->ITKContoursWithDeformationAndSegmentationError.at(i-1)));
			}

			this->DetermineResultTargetLocations(i-1);
			this->ITKContoursRegistrationResults.push_back(this->DeformImage(this->ITKContours.at(i-1), this->RegistrationResultDeformationFields.at(i-1)));

			ContourStatistics currentRegistrationResultsContourStatistics;
			this->GatherContourStatistics(this->ITKContoursRegistrationResults.at(i-1), &currentRegistrationResultsContourStatistics);
			this->RegistrationResultsContourStatistics.push_back(currentRegistrationResultsContourStatistics);

			this->ITKContoursRegistrationDifferences.push_back(this->Subtract(this->ITKContoursRegistrationResults.at(i-1), this->ITKContours.at(i-1)));
			this->RecordTargetingErrors(i-1);
		}

		int timeEndCurrentCase = time(NULL);

		if(SIMULATION_CASES > 1)
		{
			std::cout<<"Done generating simulation case "<<i<<std::endl<<std::endl;
			
			int timeDiff = timeEndCurrentCase - timeBeginCurrentCase;

			if(timeDiff < 60)
			{
				std::cout<<"Time: "<<timeDiff<<" seconds."<<std::endl;
			}
			else
			{
				std::cout<<"Time: "<<this->GetTimeString(timeDiff).c_str()<<std::endl;
			}
		}

		std::cout<<"*************************************"<<std::endl<<std::endl<<std::endl;
	}

	ReportStruct finalReport;

	if(DO_REGISTRATION)
	{
		this->RecordFinalTargetingErrorStatistics();
		finalReport = this->GetFinalReport();
	}

	int timeEnd = time(NULL);
	int timeDiff = timeEnd - timeBegin;

	std::cout<<"Done simulation"<<std::endl;

	if(timeDiff < 60)
	{
		std::cout<<"Total Time: "<<timeDiff<<" seconds."<<std::endl<<std::endl;
	}
	else
	{
		std::cout<<"Total Time: "<<this->GetTimeString(timeDiff).c_str()<<std::endl<<std::endl;
	}

	std::cout<<"*************************************"<<std::endl<<std::endl;

	if(DO_REGISTRATION)
	{
		// Display Final Report
		std::cout<<"Total Targets: "<<finalReport.totalTargets<<std::endl;
		std::cout<<"Percent Insignificant: "<<this->GetStringFromDouble(finalReport.percentInsignificant, 1)<<std::endl;
		std::cout<<"Percent Significant: "<<this->GetStringFromDouble(finalReport.percentSignificant, 1)<<std::endl;
		std::cout<<"Global Min Targeting Error: "<<this->GetStringFromDouble(finalReport.minError, 1)<<" mm"<<std::endl;
		std::cout<<"Global Max Targeting Error: "<<this->GetStringFromDouble(finalReport.maxError, 1)<<" mm"<<std::endl;
		std::cout<<"Global Mean Targeting Error: "<<this->GetStringFromDouble(finalReport.meanError, 1)<<" mm"<<std::endl;
		std::cout<<"Global Mean Deviation: "<<this->GetStringFromDouble(finalReport.meanDeviation, 1)<<" mm"<<std::endl;
		std::cout<<"Global Standard Deviation: "<<finalReport.standardDeviation<<std::endl;
		std::cout<<std::endl<<"*************************************"<<std::endl<<std::endl;
	}

	//set the scale range
	this->RemoveCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);
	this->ImageScale->SetStateToNormal();
	this->ImageScale->SetRange(0, SIMULATION_CASES);
	this->ImageScale->SetValue(0);
	this->AddCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);

	delete [] Suffix;
	delete [] line;

	//display mean shape
	this->DisplaySimulationCase(0);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//Main event handler
void vtkKWProstateErrorMapRenderingWidget::ProcessCallbackCommandEvents(
  vtkObject *caller, unsigned long event, void *calldata)
{
	vtkKWApplication *app = vtkKWApplication::SafeDownCast(
	this->GetApplication());

	if((caller == this->FileMenu) && (event == vtkKWMenu::MenuItemInvokedEvent))
	{
		if(this->FileMenu->GetItemSelectedState(this->FileMenu->GetIndexOfItem(GENERATE_MEAN_SHAPE_MENU_STRING))) //chose to generate mean shape
		{
			this->GenerateMeanShape();
			this->FileMenu->SetItemSelectedState(this->FileMenu->GetIndexOfItem(GENERATE_MEAN_SHAPE_MENU_STRING), 0);
		}
		else if(this->FileMenu->GetItemSelectedState(this->FileMenu->GetIndexOfItem(LOAD_PROSTATE_CONTOUR_MENU_STRING))) //chose to load prostate volume
		{
			this->LoadProstateContour();
			this->FileMenu->SetItemSelectedState(this->FileMenu->GetIndexOfItem(LOAD_PROSTATE_CONTOUR_MENU_STRING), 0);
		}
		else if (this->FileMenu->GetItemSelectedState(this->FileMenu->GetIndexOfItem(EXIT_MENU_STRING)))//chose to exit
		{
			this->ParentWindow->Close();
		}
		
	}
	else if((caller == this->DataMenu) && (event == vtkKWMenu::MenuItemInvokedEvent))
	{
		if(this->DataMenu->GetItemSelectedState(this->DataMenu->GetIndexOfItem(RUN_SIMULATION_MENU_STRING))) //chose to run simulation
		{
			this->RunSimulation();
			this->DataMenu->SetItemSelectedState(this->DataMenu->GetIndexOfItem(RUN_SIMULATION_MENU_STRING), 0);
		}
	}
	else if(caller == this->ImageScale && event == vtkKWScale::ScaleValueChangingEvent)
	{
		this->DisplaySimulationCase(this->ImageScale->GetValue());
	}
	

	this->Superclass::ProcessCallbackCommandEvents(caller, event, calldata);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

void vtkKWProstateErrorMapRenderingWidget::AddObservers()
{
	this->AddCallbackCommandObserver(this->FileMenu, vtkKWMenu::MenuItemInvokedEvent);
	this->AddCallbackCommandObserver(this->DataMenu, vtkKWMenu::MenuItemInvokedEvent);
	this->AddCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

void vtkKWProstateErrorMapRenderingWidget::SetParentWindow(vtkKWWindowBase* w)
{
	this->ParentWindow = w;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function sets up the GUI to display one simulation case (i.e. with render widgets for each of pre-op and intra-op contours as well
//as contours showing the difference with the ground truth).
void vtkKWProstateErrorMapRenderingWidget::SetupGUIWithMultipleImageFrames(int scaleValue)
{
	this->multipleImageFramesShowing = true;

	//delete widgets in the single-image GUI
	char* renderWidgetName = new char[strlen(this->RenderWidget->GetWidgetName()) + 1];
	strcpy(renderWidgetName, this->RenderWidget->GetWidgetName());
	this->RenderWidget->SetParent(NULL);
	this->RenderWidget->RemoveAllChildren();
	this->RenderWidget->Delete();
	this->RenderWidget = NULL;

	const char* imageFrameName = this->ImageFrame->GetWidgetName();
	this->ImageFrame->Delete();
	this->ImageFrame = NULL;

	this->RemoveCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);
	const char* imageScaleName = this->ImageScale->GetWidgetName();
	this->ImageScale->Delete();
	this->ImageScale = NULL;

	//remove widgets from the GUI
	this->Script("destroy %s", (const char*)renderWidgetName);
	this->Script("destroy %s", imageFrameName);
	this->Script("destroy %s", imageScaleName);

	//Create the frames
	this->GlobalFrame = vtkKWFrame::New();
	this->GlobalFrame->SetParent(this);
	this->GlobalFrame->Create();

	this->DoubleImageFrame = vtkKWFrame::New();
	this->DoubleImageFrame->SetParent(this->GlobalFrame);
	this->DoubleImageFrame->Create();

	this->DeformationDoubleImageFrame = vtkKWFrame::New();
	this->DeformationDoubleImageFrame->SetParent(this->GlobalFrame);
	this->DeformationDoubleImageFrame->Create();

	this->ImageFrame = vtkKWFrame::New();
	this->ImageFrame->SetParent(this->DoubleImageFrame);
	this->ImageFrame->Create();

	this->SegmentationErrorImageFrame = vtkKWFrame::New();
	this->SegmentationErrorImageFrame->SetParent(this->DoubleImageFrame);
	this->SegmentationErrorImageFrame->Create();

	this->SegmentationErrorDifferenceImageFrame = vtkKWFrame::New();
	this->SegmentationErrorDifferenceImageFrame->SetParent(this->DoubleImageFrame);
	this->SegmentationErrorDifferenceImageFrame->Create();

	this->DeformationImageFrame = vtkKWFrame::New();
	this->DeformationImageFrame->SetParent(this->DeformationDoubleImageFrame);
	this->DeformationImageFrame->Create();

	this->DeformationSegmentationErrorImageFrame = vtkKWFrame::New();
	this->DeformationSegmentationErrorImageFrame->SetParent(this->DeformationDoubleImageFrame);
	this->DeformationSegmentationErrorImageFrame->Create();

	this->DeformationSegmentationErrorDifferenceImageFrame = vtkKWFrame::New();
	this->DeformationSegmentationErrorDifferenceImageFrame->SetParent(this->DeformationDoubleImageFrame);
	this->DeformationSegmentationErrorDifferenceImageFrame->Create();

	this->RegistrationResultImageFrame = vtkKWFrame::New();
	this->RegistrationResultImageFrame->SetParent(this->GlobalFrame);
	this->RegistrationResultImageFrame->Create();
	this->RegistrationResultImageFrame->SetBackgroundColor(0,0,0);

	this->RegistrationDifferenceImageFrame = vtkKWFrame::New();
	this->RegistrationDifferenceImageFrame->SetParent(this->GlobalFrame);
	this->RegistrationDifferenceImageFrame->Create();
	this->RegistrationDifferenceImageFrame->SetBackgroundColor(0,0,0);

	this->LabelFrame = vtkKWFrame::New();
	this->LabelFrame->SetParent(this);
	this->LabelFrame->Create();

	this->InfoLabelFrame = vtkKWFrame::New();
	this->InfoLabelFrame->SetParent(this->LabelFrame);
	this->InfoLabelFrame->Create();

	this->ReportLabelFrame = vtkKWFrame::New();
	this->ReportLabelFrame->SetParent(this->LabelFrame);
	this->ReportLabelFrame->Create();

	//create the text labels
	this->FirstModeLabel = vtkKWLabel::New();
	this->FirstModeLabel->SetParent(this->InfoLabelFrame);
	this->FirstModeLabel->Create();
	this->FirstModeLabel->SetText("First Mode");

	this->SecondModeLabel = vtkKWLabel::New();
	this->SecondModeLabel->SetParent(this->InfoLabelFrame);
	this->SecondModeLabel->Create();
	this->SecondModeLabel->SetText("Second Mode");

	this->ThirdModeLabel = vtkKWLabel::New();
	this->ThirdModeLabel->SetParent(this->InfoLabelFrame);
	this->ThirdModeLabel->Create();
	this->ThirdModeLabel->SetText("Third Mode");

	this->BaseSegmentationErrorLabel = vtkKWLabel::New();
	this->BaseSegmentationErrorLabel->SetParent(this->InfoLabelFrame);
	this->BaseSegmentationErrorLabel->Create();
	this->BaseSegmentationErrorLabel->SetText("Base Segmentation Error");

	this->MidSegmentationErrorLabel = vtkKWLabel::New();
	this->MidSegmentationErrorLabel->SetParent(this->InfoLabelFrame);
	this->MidSegmentationErrorLabel->Create();
	this->MidSegmentationErrorLabel->SetText("Mid Segmentation Error");

	this->ApexSegmentationErrorLabel = vtkKWLabel::New();
	this->ApexSegmentationErrorLabel->SetParent(this->InfoLabelFrame);
	this->ApexSegmentationErrorLabel->Create();
	this->ApexSegmentationErrorLabel->SetText("Apex Segmentation Error");

	this->TotalTargetsLabel = vtkKWLabel::New();
	this->TotalTargetsLabel->SetParent(this->ReportLabelFrame);
	this->TotalTargetsLabel->Create();
	this->TotalTargetsLabel->SetText("Total Targets");

	this->InsignificantLabel = vtkKWLabel::New();
	this->InsignificantLabel->SetParent(this->ReportLabelFrame);
	this->InsignificantLabel->Create();
	this->InsignificantLabel->SetText("Total Insignificant");

	this->SignificantLabel = vtkKWLabel::New();
	this->SignificantLabel->SetParent(this->ReportLabelFrame);
	this->SignificantLabel->Create();
	this->SignificantLabel->SetText("Total Significant");

	this->MinErrorLabel = vtkKWLabel::New();
	this->MinErrorLabel->SetParent(this->ReportLabelFrame);
	this->MinErrorLabel->Create();
	this->MinErrorLabel->SetText("Min Targeting Error");

	this->MaxErrorLabel = vtkKWLabel::New();
	this->MaxErrorLabel->SetParent(this->ReportLabelFrame);
	this->MaxErrorLabel->Create();
	this->MaxErrorLabel->SetText("Max Targeting Error");

	this->MeanErrorLabel = vtkKWLabel::New();
	this->MeanErrorLabel->SetParent(this->ReportLabelFrame);
	this->MeanErrorLabel->Create();
	this->MeanErrorLabel->SetText("Mean Error");

	this->MeanDeviationLabel = vtkKWLabel::New();
	this->MeanDeviationLabel->SetParent(this->ReportLabelFrame);
	this->MeanDeviationLabel->Create();
	this->MeanDeviationLabel->SetText("Mean Deviation");

	this->StandardDeviationLabel = vtkKWLabel::New();
	this->StandardDeviationLabel->SetParent(this->ReportLabelFrame);
	this->StandardDeviationLabel->Create();
	this->StandardDeviationLabel->SetText("Standard Deviation");

	//Add a render widget and attach it to the image frame
	this->RenderWidget = vtkKWRenderWidget::New();
    this->RenderWidget->SetParent(this->ImageFrame);
    this->RenderWidget->Create();
	this->RenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->RenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->RenderWidget->CornerAnnotationVisibilityOn();
	this->RenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->RenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add a render widget and attach it to the segmentation error image frame
	this->SegmentationErrorRenderWidget = vtkKWRenderWidget::New();
    this->SegmentationErrorRenderWidget->SetParent(this->SegmentationErrorImageFrame);
    this->SegmentationErrorRenderWidget->Create();
	this->SegmentationErrorRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->SegmentationErrorRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->SegmentationErrorRenderWidget->CornerAnnotationVisibilityOn();
	this->SegmentationErrorRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->SegmentationErrorRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add a render widget and attach it to the segmentation error difference image frame
	this->SegmentationErrorDifferenceRenderWidget = vtkKWRenderWidget::New();
    this->SegmentationErrorDifferenceRenderWidget->SetParent(this->SegmentationErrorDifferenceImageFrame);
    this->SegmentationErrorDifferenceRenderWidget->Create();
	this->SegmentationErrorDifferenceRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->SegmentationErrorDifferenceRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->SegmentationErrorDifferenceRenderWidget->CornerAnnotationVisibilityOn();
	this->SegmentationErrorDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->SegmentationErrorDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add a render widget and attach it to the deformation image frame
	this->DeformationRenderWidget = vtkKWRenderWidget::New();
    this->DeformationRenderWidget->SetParent(this->DeformationImageFrame);
    this->DeformationRenderWidget->Create();
	this->DeformationRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->DeformationRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->DeformationRenderWidget->CornerAnnotationVisibilityOn();
	this->DeformationRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->DeformationRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add a render widget and attach it to the deformation segmentation error image frame
	this->DeformationSegmentationErrorRenderWidget = vtkKWRenderWidget::New();
    this->DeformationSegmentationErrorRenderWidget->SetParent(this->DeformationSegmentationErrorImageFrame);
    this->DeformationSegmentationErrorRenderWidget->Create();
	this->DeformationSegmentationErrorRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->DeformationSegmentationErrorRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->DeformationSegmentationErrorRenderWidget->CornerAnnotationVisibilityOn();
	this->DeformationSegmentationErrorRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->DeformationSegmentationErrorRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add a render widget and attach it to the deformation segmentation error difference image frame
	this->DeformationSegmentationErrorDifferenceRenderWidget = vtkKWRenderWidget::New();
    this->DeformationSegmentationErrorDifferenceRenderWidget->SetParent(this->DeformationSegmentationErrorDifferenceImageFrame);
    this->DeformationSegmentationErrorDifferenceRenderWidget->Create();
	this->DeformationSegmentationErrorDifferenceRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->DeformationSegmentationErrorDifferenceRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->DeformationSegmentationErrorDifferenceRenderWidget->CornerAnnotationVisibilityOn();
	this->DeformationSegmentationErrorDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->DeformationSegmentationErrorDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add the registration result render widget and attach it to the registration result image frame
	this->RegistrationResultRenderWidget = vtkKWRenderWidget::New();
    this->RegistrationResultRenderWidget->SetParent(this->RegistrationResultImageFrame);
    this->RegistrationResultRenderWidget->Create();
	this->RegistrationResultRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->RegistrationResultRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->RegistrationResultRenderWidget->CornerAnnotationVisibilityOn();
	this->RegistrationResultRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->RegistrationResultRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//Add the registration difference render widget and attach it to the Registration difference image frame
	this->RegistrationDifferenceRenderWidget = vtkKWRenderWidget::New();
    this->RegistrationDifferenceRenderWidget->SetParent(this->RegistrationDifferenceImageFrame);
    this->RegistrationDifferenceRenderWidget->Create();
	this->RegistrationDifferenceRenderWidget->SetWidth(SMALL_RENDER_WIDGET_WIDTH);
	this->RegistrationDifferenceRenderWidget->SetHeight(SMALL_RENDER_WIDGET_HEIGHT);
    this->RegistrationDifferenceRenderWidget->CornerAnnotationVisibilityOn();
	this->RegistrationDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_MULTIPLE);
	this->RegistrationDifferenceRenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//pack render widgets
    this->Script("pack %s -expand n -fill both", 
              this->RenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->SegmentationErrorRenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->SegmentationErrorDifferenceRenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->DeformationRenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->DeformationSegmentationErrorRenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->DeformationSegmentationErrorDifferenceRenderWidget->GetWidgetName());

	this->Script("pack %s -expand n -pady %i", 
              this->RegistrationResultRenderWidget->GetWidgetName(), SMALL_RENDER_WIDGET_HEIGHT);

	this->Script("pack %s -expand n -pady %i", 
              this->RegistrationDifferenceRenderWidget->GetWidgetName(), SMALL_RENDER_WIDGET_HEIGHT);

	//re-create the image scale
	this->ImageScale = vtkKWScale::New();
	this->ImageScale->SetParent(this);
	this->ImageScale->Create();
	this->ImageScale->SetResolution(1.0);
	this->ImageScale->SetRange(0, SIMULATION_CASES);
	this->ImageScale->SetValue(scaleValue);

	//pack frames
	this->Script("pack %s -expand n -fill both", 
              this->GlobalFrame->GetWidgetName());

	this->Script("pack %s -expand n -side left -ipady 0", 
              this->DoubleImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side left -ipady 0", 
              this->DeformationDoubleImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -pady 1", 
              this->ImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -padx 1", 
              this->SegmentationErrorImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -padx 1 -pady 1", 
              this->SegmentationErrorDifferenceImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -pady 1", 
              this->DeformationImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -padx 1", 
              this->DeformationSegmentationErrorImageFrame->GetWidgetName());

	this->Script("pack %s -expand n -side top -ipady 0 -padx 1 -pady 1", 
              this->DeformationSegmentationErrorDifferenceImageFrame->GetWidgetName());

	this->Script("pack %s -expand y -side left -ipady 0 -padx 1", 
              this->RegistrationResultImageFrame->GetWidgetName());

	this->Script("pack %s -expand y -side left -ipady 0 -padx 1", 
              this->RegistrationDifferenceImageFrame->GetWidgetName());

	//pack the labels
    this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->FirstModeLabel->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->SecondModeLabel->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->ThirdModeLabel->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->BaseSegmentationErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->MidSegmentationErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -ipady 0",
			  this->ApexSegmentationErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->TotalTargetsLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->InsignificantLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->SignificantLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->MinErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->MaxErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->MeanErrorLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->MeanDeviationLabel->GetWidgetName());

	this->Script("pack %s -expand n -side top -anchor nw -ipady 0",
			  this->StandardDeviationLabel->GetWidgetName());

	this->Script("pack %s -expand n -fill both", 
              this->LabelFrame->GetWidgetName());

	this->Script("pack %s -expand n -anchor nw -side left -ipady 0", 
              this->InfoLabelFrame->GetWidgetName());

	this->Script("pack %s -expand n -side left -ipady 0 -padx 17", 
              this->ReportLabelFrame->GetWidgetName());

	//pack the scale
	this->Script("pack %s -expand n -fill both", 
              this->ImageScale->GetWidgetName());

	this->AddCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);

	delete [] renderWidgetName;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function resets the GUI to one image frame (i.e. for showing only the mean shape)
void vtkKWProstateErrorMapRenderingWidget::ResetGUIToOneImageFrame()
{
	//delete widgets in the multi-image GUI
	char* renderWidgetName = new char[strlen(this->RenderWidget->GetWidgetName()) + 1];
	strcpy(renderWidgetName, this->RenderWidget->GetWidgetName());
	this->RenderWidget->SetParent(NULL);
	this->RenderWidget->RemoveAllChildren();
	this->RenderWidget->Delete();
	this->RenderWidget = NULL;

	char* segmentationErrorRenderWidgetName = new char[strlen(this->SegmentationErrorRenderWidget->GetWidgetName()) + 1];
	strcpy(segmentationErrorRenderWidgetName, this->SegmentationErrorRenderWidget->GetWidgetName());
	this->SegmentationErrorRenderWidget->SetParent(NULL);
	this->SegmentationErrorRenderWidget->RemoveAllChildren();
	this->SegmentationErrorRenderWidget->Delete();
	this->SegmentationErrorRenderWidget = NULL;

	char* segmentationErrorDifferenceRenderWidgetName = new char[strlen(this->SegmentationErrorDifferenceRenderWidget->GetWidgetName()) + 1];
	strcpy(segmentationErrorDifferenceRenderWidgetName, this->SegmentationErrorDifferenceRenderWidget->GetWidgetName());
	this->SegmentationErrorDifferenceRenderWidget->SetParent(NULL);
	this->SegmentationErrorDifferenceRenderWidget->RemoveAllChildren();
	this->SegmentationErrorDifferenceRenderWidget->Delete();
	this->SegmentationErrorDifferenceRenderWidget = NULL;

	char* deformationRenderWidgetName = new char[strlen(this->DeformationRenderWidget->GetWidgetName()) + 1];
	strcpy(deformationRenderWidgetName, this->DeformationRenderWidget->GetWidgetName());
	this->DeformationRenderWidget->SetParent(NULL);
	this->DeformationRenderWidget->RemoveAllChildren();
	this->DeformationRenderWidget->Delete();
	this->DeformationRenderWidget = NULL;

	char* deformationSegmentationErrorRenderWidgetName = new char[strlen(this->DeformationSegmentationErrorRenderWidget->GetWidgetName()) + 1];
	strcpy(deformationSegmentationErrorRenderWidgetName, this->DeformationSegmentationErrorRenderWidget->GetWidgetName());
	this->DeformationSegmentationErrorRenderWidget->SetParent(NULL);
	this->DeformationSegmentationErrorRenderWidget->RemoveAllChildren();
	this->DeformationSegmentationErrorRenderWidget->Delete();
	this->DeformationSegmentationErrorRenderWidget = NULL;

	char* deformationSegmentationErrorDifferenceRenderWidgetName = new char[strlen(this->DeformationSegmentationErrorDifferenceRenderWidget->GetWidgetName()) + 1];
	strcpy(deformationSegmentationErrorDifferenceRenderWidgetName, this->DeformationSegmentationErrorDifferenceRenderWidget->GetWidgetName());
	this->DeformationSegmentationErrorDifferenceRenderWidget->SetParent(NULL);
	this->DeformationSegmentationErrorDifferenceRenderWidget->RemoveAllChildren();
	this->DeformationSegmentationErrorDifferenceRenderWidget->Delete();
	this->DeformationSegmentationErrorDifferenceRenderWidget = NULL;

	char* registrationResultRenderWidgetName = new char[strlen(this->RegistrationResultRenderWidget->GetWidgetName()) + 1];
	strcpy(registrationResultRenderWidgetName, this->RegistrationResultRenderWidget->GetWidgetName());
	this->RegistrationResultRenderWidget->SetParent(NULL);
	this->RegistrationResultRenderWidget->RemoveAllChildren();
	this->RegistrationResultRenderWidget->Delete();
	this->RegistrationResultRenderWidget = NULL;

	char* registrationDifferenceRenderWidgetName = new char[strlen(this->RegistrationDifferenceRenderWidget->GetWidgetName()) + 1];
	strcpy(registrationDifferenceRenderWidgetName, this->RegistrationDifferenceRenderWidget->GetWidgetName());
	this->RegistrationDifferenceRenderWidget->SetParent(NULL);
	this->RegistrationDifferenceRenderWidget->RemoveAllChildren();
	this->RegistrationDifferenceRenderWidget->Delete();
	this->RegistrationDifferenceRenderWidget = NULL;

	const char* firstModeLabelName = this->FirstModeLabel->GetWidgetName();
	this->FirstModeLabel->Delete();
	this->FirstModeLabel = NULL;

	const char* secondModeLabelName = this->SecondModeLabel->GetWidgetName();
	this->SecondModeLabel->Delete();
	this->SecondModeLabel = NULL;

	const char* thirdModeLabelName = this->ThirdModeLabel->GetWidgetName();
	this->ThirdModeLabel->Delete();
	this->ThirdModeLabel = NULL;

	const char* baseSegmentationErrorLabelName = this->BaseSegmentationErrorLabel->GetWidgetName();
	this->BaseSegmentationErrorLabel->Delete();
	this->BaseSegmentationErrorLabel = NULL;

	const char* midSegmentationErrorLabelName = this->MidSegmentationErrorLabel->GetWidgetName();
	this->MidSegmentationErrorLabel->Delete();
	this->MidSegmentationErrorLabel = NULL;

	const char* apexSegmentationErrorLabelName = this->ApexSegmentationErrorLabel->GetWidgetName();
	this->ApexSegmentationErrorLabel->Delete();
	this->ApexSegmentationErrorLabel = NULL;

	const char* totalTargetsLabelName = this->TotalTargetsLabel->GetWidgetName();
	this->TotalTargetsLabel->Delete();
	this->TotalTargetsLabel = NULL;

	const char* insignificantLabelName = this->InsignificantLabel->GetWidgetName();
	this->InsignificantLabel->Delete();
	this->InsignificantLabel = NULL;

	const char* significantLabelName = this->SignificantLabel->GetWidgetName();
	this->SignificantLabel->Delete();
	this->SignificantLabel = NULL;

	const char* minErrorLabelName = this->MinErrorLabel->GetWidgetName();
	this->MinErrorLabel->Delete();
	this->MinErrorLabel = NULL;

	const char* maxErrorLabelName = this->MaxErrorLabel->GetWidgetName();
	this->MaxErrorLabel->Delete();
	this->MaxErrorLabel = NULL;

	const char* meanErrorLabelName = this->MeanErrorLabel->GetWidgetName();
	this->MeanErrorLabel->Delete();
	this->MeanErrorLabel = NULL;

	const char* meanDeviationLabelName = this->MeanDeviationLabel->GetWidgetName();
	this->MeanDeviationLabel->Delete();
	this->MeanDeviationLabel = NULL;

	const char* standardDeviationLabelName = this->StandardDeviationLabel->GetWidgetName();
	this->StandardDeviationLabel->Delete();
	this->StandardDeviationLabel = NULL;

	const char* doubleImageFrameName = this->DoubleImageFrame->GetWidgetName();
	this->DoubleImageFrame->Delete();
	this->DoubleImageFrame = NULL;

	const char* deformationDoubleImageFrameName = this->DeformationDoubleImageFrame->GetWidgetName();
	this->DeformationDoubleImageFrame->Delete();
	this->DeformationDoubleImageFrame = NULL;

	const char* imageFrameName = this->ImageFrame->GetWidgetName();
	this->ImageFrame->Delete();
	this->ImageFrame = NULL;

	const char* segmentationErrorImageFrameName = this->SegmentationErrorImageFrame->GetWidgetName();
	this->SegmentationErrorImageFrame->Delete();
	this->SegmentationErrorImageFrame = NULL;

	const char* segmentationErrorDifferenceImageFrameName = this->SegmentationErrorDifferenceImageFrame->GetWidgetName();
	this->SegmentationErrorDifferenceImageFrame->Delete();
	this->SegmentationErrorDifferenceImageFrame = NULL;

	const char* deformationImageFrameName = this->DeformationImageFrame->GetWidgetName();
	this->DeformationImageFrame->Delete();
	this->DeformationImageFrame = NULL;

	const char* deformationSegmentationErrorImageFrameName = this->DeformationSegmentationErrorImageFrame->GetWidgetName();
	this->DeformationSegmentationErrorImageFrame->Delete();
	this->DeformationSegmentationErrorImageFrame = NULL;

	const char* deformationSegmentationErrorDifferenceImageFrameName = this->DeformationSegmentationErrorDifferenceImageFrame->GetWidgetName();
	this->DeformationSegmentationErrorDifferenceImageFrame->Delete();
	this->DeformationSegmentationErrorDifferenceImageFrame = NULL;

	const char* registrationResultImageFrameName = this->RegistrationResultImageFrame->GetWidgetName();
	this->RegistrationResultImageFrame->Delete();
	this->RegistrationResultImageFrame = NULL;

	const char* registrationDifferenceImageFrameName = this->RegistrationDifferenceImageFrame->GetWidgetName();
	this->RegistrationDifferenceImageFrame->Delete();
	this->RegistrationDifferenceImageFrame = NULL;

	const char* globalFrameName = this->GlobalFrame->GetWidgetName();
	this->GlobalFrame->Delete();
	this->GlobalFrame = NULL;

	const char* infoLabelFrameName = this->InfoLabelFrame->GetWidgetName();
	this->InfoLabelFrame->Delete();
	this->InfoLabelFrame = NULL;

	const char* reportLabelFrameName = this->ReportLabelFrame->GetWidgetName();
	this->ReportLabelFrame->Delete();
	this->ReportLabelFrame = NULL;

	const char* labelFrameName = this->LabelFrame->GetWidgetName();
	this->LabelFrame->Delete();
	this->LabelFrame = NULL;

	this->RemoveCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);
	const char* imageScaleName = this->ImageScale->GetWidgetName();
	this->ImageScale->Delete();
	this->ImageScale = NULL;

	//remove widgets from the GUI
	this->Script("destroy %s", (const char*)renderWidgetName);
	this->Script("destroy %s", (const char*)segmentationErrorRenderWidgetName);
	this->Script("destroy %s", (const char*)segmentationErrorDifferenceRenderWidgetName);
	this->Script("destroy %s", (const char*)deformationRenderWidgetName);
	this->Script("destroy %s", (const char*)deformationSegmentationErrorRenderWidgetName);
	this->Script("destroy %s", (const char*)deformationSegmentationErrorDifferenceRenderWidgetName);
	this->Script("destroy %s", (const char*)registrationResultRenderWidgetName);
	this->Script("destroy %s", (const char*)registrationDifferenceRenderWidgetName);
	this->Script("destroy %s", firstModeLabelName);
	this->Script("destroy %s", secondModeLabelName);
	this->Script("destroy %s", thirdModeLabelName);
	this->Script("destroy %s", baseSegmentationErrorLabelName);
	this->Script("destroy %s", midSegmentationErrorLabelName);
	this->Script("destroy %s", apexSegmentationErrorLabelName);
	this->Script("destroy %s", totalTargetsLabelName);
	this->Script("destroy %s", insignificantLabelName);
	this->Script("destroy %s", significantLabelName);
	this->Script("destroy %s", minErrorLabelName);
	this->Script("destroy %s", maxErrorLabelName);
	this->Script("destroy %s", meanErrorLabelName);
	this->Script("destroy %s", meanDeviationLabelName);
	this->Script("destroy %s", standardDeviationLabelName);
	this->Script("destroy %s", doubleImageFrameName);
	this->Script("destroy %s", deformationDoubleImageFrameName);
	this->Script("destroy %s", segmentationErrorImageFrameName);
	this->Script("destroy %s", segmentationErrorDifferenceImageFrameName);
	this->Script("destroy %s", imageFrameName);
	this->Script("destroy %s", deformationSegmentationErrorImageFrameName);
	this->Script("destroy %s", deformationSegmentationErrorDifferenceImageFrameName);
	this->Script("destroy %s", deformationImageFrameName);
	this->Script("destroy %s", registrationResultImageFrameName);
	this->Script("destroy %s", registrationDifferenceImageFrameName);
	this->Script("destroy %s", globalFrameName);
	this->Script("destroy %s", infoLabelFrameName);
	this->Script("destroy %s", reportLabelFrameName);
	this->Script("destroy %s", labelFrameName);
	this->Script("destroy %s", imageScaleName);

	this->AddDefaultGUIObjects();

	this->ImageScale->SetStateToNormal();

	this->AddCallbackCommandObserver(this->ImageScale, vtkKWScale::ScaleValueChangingEvent);

	delete [] renderWidgetName;
	delete [] segmentationErrorRenderWidgetName;
	delete [] segmentationErrorDifferenceRenderWidgetName;
	delete [] deformationRenderWidgetName;
	delete [] deformationSegmentationErrorRenderWidgetName;
	delete [] deformationSegmentationErrorDifferenceRenderWidgetName;
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

//This function adds the GUI objects which aren't destroyed when switching between views
void vtkKWProstateErrorMapRenderingWidget::AddDefaultGUIObjects()
{
	this->multipleImageFramesShowing = false;

	//Create the frames
	this->ImageFrame = vtkKWFrame::New();
	this->ImageFrame->SetParent(this);
	this->ImageFrame->Create();

	// Add a render widget, attach it to the image frame, and pack
	this->RenderWidget = vtkKWRenderWidget::New();
    this->RenderWidget->SetParent(this->ImageFrame);
    this->RenderWidget->Create();
	this->RenderWidget->SetWidth(MAIN_RENDER_WIDGET_WIDTH);
	this->RenderWidget->SetHeight(MAIN_RENDER_WIDGET_HEIGHT);
    this->RenderWidget->CornerAnnotationVisibilityOn();
	this->RenderWidget->GetRenderer()->GetActiveCamera()->Azimuth(RENDER_WIDGET_INITIAL_AZIMUTH_SINGLE);
	this->RenderWidget->GetRenderer()->GetActiveCamera()->Elevation(RENDER_WIDGET_INITIAL_ELEVATION);

	//pack render widget
    this->Script("pack %s -expand n -fill both -anchor n", 
              this->RenderWidget->GetWidgetName());

	//create the image scale
	this->ImageScale = vtkKWScale::New();
	this->ImageScale->SetParent(this->ImageFrame);
	this->ImageScale->Create();
	this->ImageScale->SetResolution(1.0);
	this->ImageScale->SetRange(0, SIMULATION_CASES);
	this->ImageScale->SetStateToDisabled();

	//pack the scale
	this->Script("pack %s -anchor s -fill both", 
              this->ImageScale->GetWidgetName());

	//pack frames
	this->Script("pack %s -side left -expand n -fill both -ipadx 0 -padx 1 -anchor nw", 
              this->ImageFrame->GetWidgetName());
}

/******************************************************************************************************************************************
******************************************************************************************************************************************/

void vtkKWProstateErrorMapRenderingWidget::CreateWidget()
{
	// Check if already created
	if (this->IsCreated())
	{
		vtkErrorMacro("class already created");
		return;
	}

	//Call the superclass to create the whole widget
	this->Superclass::CreateWidget();

	//create the menus
	this->FileMenu = vtkKWMenu::New();
	this->FileMenu->SetParent(this->ParentWindow);
	this->FileMenu->Create();
	this->FileMenu->AddCheckButton(GENERATE_MEAN_SHAPE_MENU_STRING);
	this->FileMenu->AddCheckButton(LOAD_PROSTATE_CONTOUR_MENU_STRING);
	this->FileMenu->AddSeparator();
	this->FileMenu->AddCheckButton(EXIT_MENU_STRING);
	this->ParentWindow->GetMenu()->AddCascade("File", this->FileMenu);

	this->DataMenu = vtkKWMenu::New();
	this->DataMenu->SetParent(this->ParentWindow);
	this->DataMenu->Create();
	this->DataMenu->AddCheckButton(RUN_SIMULATION_MENU_STRING);
	this->ParentWindow->GetMenu()->AddCascade("Data", this->DataMenu);

	this->AddDefaultGUIObjects();

	this->AddObservers();

	//set strings to let user know how much segmentation
	//error is applied to each section of the prostate
	baseSegmentationErrorString = "Base Segmentation Error: ";
	baseSegmentationErrorString.append("Uniform Random on [");
	baseSegmentationErrorString.append(this->GetStringFromDouble(MIN_BASE_SEGMENTATION_ERROR, 1));
	baseSegmentationErrorString.append(" mm, ");
	baseSegmentationErrorString.append(this->GetStringFromDouble(MAX_BASE_SEGMENTATION_ERROR, 1));
	baseSegmentationErrorString.append(" mm]");

	midSegmentationErrorString = "Mid Segmentation Error: ";
	midSegmentationErrorString.append("Uniform Random on [");
	midSegmentationErrorString.append(this->GetStringFromDouble(MIN_MID_SEGMENTATION_ERROR, 1));
	midSegmentationErrorString.append(" mm, ");
	midSegmentationErrorString.append(this->GetStringFromDouble(MAX_MID_SEGMENTATION_ERROR, 1));
	midSegmentationErrorString.append(" mm]");

	apexSegmentationErrorString = "Apex Segmentation Error: ";
	apexSegmentationErrorString.append("Uniform Random on [");
	apexSegmentationErrorString.append(this->GetStringFromDouble(MIN_APEX_SEGMENTATION_ERROR, 1));
	apexSegmentationErrorString.append(" mm, ");
	apexSegmentationErrorString.append(this->GetStringFromDouble(MAX_APEX_SEGMENTATION_ERROR, 1));
	apexSegmentationErrorString.append(" mm]");
}