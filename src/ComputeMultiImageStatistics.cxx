#include "DefRegEval.h"
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkSubtractImageFilter.h"
#include "itkGradientToMagnitudeImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleFunctionErodeImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMaximumImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkAddImageFilter.h"
//#include "itkScalarImageToHistogramGenerator.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkSampleToHistogramFilter.h"
#include "itkHistogram.h"

typedef float  ScalarPixelType;
typedef itk::Image< ScalarPixelType, ImageDimension > ScalarImageType;


void readTextLineToListOfString(const char* textFileName, std::vector< std::string > &listOfStrings)
{
  /* The file MUST end with an empty line, then each line will be
  stored as an element in the returned vector object. */

  std::ifstream f(textFileName);
  std::string thisLine;

  if (f.good())
  {
    while( std::getline(f, thisLine) )
    {
      listOfStrings.push_back(thisLine);
    }
  }
  else
  {
    std::cerr<<"Error: can not open file:"<<textFileName<<std::endl;
    exit(EXIT_FAILURE);
  }

  f.close();

} 

void WriteHistogram(const char *outputFilename, ScalarImageType::Pointer inputImage, int histogramBinCount, double histogramMin, double histogramMax)
{  
  // It would be simpler to use ScalarImageToHistogramGenerator, but unfortunately it auto-scales the histogram bins

  typedef itk::Statistics::ImageToListSampleAdaptor<ScalarImageType>   AdaptorType;
  typedef AdaptorType::Pointer                   AdaptorPointer;
  typedef ScalarImageType::PixelType                   PixelType;
  typedef itk::NumericTraits< PixelType >::RealType   RealPixelType;

  typedef itk::Statistics::Histogram< double > HistogramType;
  typedef itk::Statistics::SampleToHistogramFilter< AdaptorType, HistogramType > GeneratorType;

  typedef GeneratorType::Pointer                   GeneratorPointer;

  typedef HistogramType::Pointer                   HistogramPointer;
  typedef HistogramType::ConstPointer              HistogramConstPointer;

  AdaptorPointer      m_ImageToListAdaptor;
  GeneratorPointer    m_HistogramGenerator;
  


  m_ImageToListAdaptor = AdaptorType::New();
  m_HistogramGenerator = GeneratorType::New();
  m_HistogramGenerator->SetInput( m_ImageToListAdaptor );

  m_HistogramGenerator->SetAutoMinimumMaximum(false);

  HistogramType::Pointer histogramGenerator = HistogramType::New();

  double histogramBinWidth=(histogramMax-histogramMin)/(histogramBinCount-1);
  
  HistogramType::SizeType size;
  size.SetSize(1);
  size.Fill( histogramBinCount );
  m_HistogramGenerator->SetHistogramSize( size );
  
  //m_HistogramGenerator->SetMarginalScale( 10.0 );

  typedef GeneratorType::HistogramMeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType minVector;
  itk::Statistics::MeasurementVectorTraits::SetLength(minVector,1);
  minVector[0] = histogramMin-histogramBinWidth/2.0;
  m_HistogramGenerator->SetHistogramBinMinimum( minVector );

  typedef GeneratorType::HistogramMeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType maxVector;
  itk::Statistics::MeasurementVectorTraits::SetLength(maxVector,1);
  maxVector[0] = histogramMax+histogramBinWidth/2.0;
  m_HistogramGenerator->SetHistogramBinMaximum( maxVector );

  m_ImageToListAdaptor->SetImage( inputImage );

  m_HistogramGenerator->Update();
  const HistogramType * histogram = m_HistogramGenerator->GetOutput();

  const unsigned int histogramSize = histogram->Size();

  std::ofstream outputFile(outputFilename);

  outputFile << "bin center, frequency, histogramSize: " << histogramSize << std::endl;

  HistogramType::ConstIterator itr = histogram->Begin();
  HistogramType::ConstIterator end = histogram->End();

  unsigned int binNumber = 0;
  while( itr != end )
  {
    outputFile << binNumber*histogramBinWidth << ", " << itr.GetFrequency() << std::endl;
    ++itr;
    ++binNumber;
  }

  outputFile.close();
}

int main( int argc, char * argv[] )
{

  std::string fileListFilename;
  std::string meanImageFilename;
  std::string maxImageFilename;
  std::string meanHistogramFilename;
  std::string maxHistogramFilename;
  bool vectorImageInput=false;
  int histogramBinCount=51;
  double histogramMin=0.0;
  double histogramMax=5.0;

  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);

  args.AddArgument("--fileList", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &fileListFilename, "Text file with the list of files to be analysed");
  args.AddArgument("--meanImage", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &meanImageFilename, "Mean image output filename");
  args.AddArgument("--maxImage", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &maxImageFilename, "Max image output filename");
  args.AddArgument("--meanHistogram", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &meanHistogramFilename, "Histogram of mean image values output filename");
  args.AddArgument("--maxHistogram", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &maxHistogramFilename, "Histogram of max image values output filename");
  args.AddArgument("--histogramBinCount", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &histogramBinCount, "Number of histogram bins");
  args.AddArgument("--histogramMin", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &histogramMin, "Histogram minimum value");
  args.AddArgument("--histogramMax", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &histogramMax, "Histogram maximum value");
  args.AddBooleanArgument("--vectorInput", &vectorImageInput, "Set it to true for processing vector images");
  
  if ( !args.Parse() )
  {
    std::cerr << "Problem parsing arguments" << std::endl;
    std::cout << "Help: " << args.GetHelp() << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector< std::string > listOfStrings;
  readTextLineToListOfString(fileListFilename.c_str(), listOfStrings);

  typedef itk::ImageFileReader< ScalarImageType > ScalarImageReaderType;

  typedef itk::Vector< ScalarPixelType, ImageDimension > VectorPixelType;
  typedef itk::Image< VectorPixelType,  ImageDimension > VectorImageType;
  typedef itk::ImageFileReader< VectorImageType > VectorImageReaderType;


  ScalarImageType::Pointer meanAccumulatorImage;
  ScalarImageType::Pointer maxAccumulatorImage;
  
  ScalarImageType::Pointer readImage;

  ScalarImageReaderType::Pointer scalarImageReader = ScalarImageReaderType::New();
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();

  typedef itk::GradientToMagnitudeImageFilter< VectorImageType, ScalarImageType > MagnitudeFilterType;
  MagnitudeFilterType::Pointer magnitudeFilter = MagnitudeFilterType::New();

  int numberOfImages=0;
  for (std::vector< std::string >::iterator it=listOfStrings.begin(); it!=listOfStrings.end(); ++it)
  {
    numberOfImages++;
    if (!vectorImageInput)
    {
      scalarImageReader->SetFileName( *it );
      try
      {
        scalarImageReader->Update();
      }
      catch( itk::ExceptionObject & excp )
      {
        std::cerr << "Exception thrown " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
      }
      scalarImageReader->Update(); 
      readImage=scalarImageReader->GetOutput();
    }
    else
    {
      vectorImageReader->SetFileName( *it );
      try
      {
        vectorImageReader->Update();
      }
      catch( itk::ExceptionObject & excp )
      {
        std::cerr << "Exception thrown " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
      }
      magnitudeFilter->SetInput(vectorImageReader->GetOutput());
      magnitudeFilter->Update(); 
      readImage=magnitudeFilter->GetOutput();
    }

    if (meanAccumulatorImage.IsNull())
    {
      typedef itk::ImageDuplicator< ScalarImageType >   DuplicatorType;
      DuplicatorType::Pointer duplicator  = DuplicatorType::New();
      duplicator->SetInputImage( readImage );      
      duplicator->Update();
      meanAccumulatorImage = duplicator->GetOutput();
    }
    else
    {
      typedef itk::ImageDuplicator< ScalarImageType >   DuplicatorType;
      DuplicatorType::Pointer duplicator  = DuplicatorType::New();
      duplicator->SetInputImage( meanAccumulatorImage );
      duplicator->Update();
      typedef itk::AddImageFilter<ScalarImageType,ScalarImageType,ScalarImageType> AddFilterType;
      AddFilterType::Pointer addFilter = AddFilterType::New();
      addFilter->SetInput1(duplicator->GetOutput());
      addFilter->SetInput2(readImage);
      addFilter->Update();
      meanAccumulatorImage = addFilter->GetOutput();
    }

    if (maxAccumulatorImage.IsNull())
    {
      typedef itk::ImageDuplicator< ScalarImageType >   DuplicatorType;
      DuplicatorType::Pointer duplicator  = DuplicatorType::New();
      duplicator->SetInputImage( readImage );
      duplicator->Update();
      maxAccumulatorImage = duplicator->GetOutput();
    }
    else
    {
      typedef itk::ImageDuplicator< ScalarImageType >   DuplicatorType;
      DuplicatorType::Pointer duplicator  = DuplicatorType::New();
      duplicator->SetInputImage( maxAccumulatorImage );
      duplicator->Update();
      typedef itk::MaximumImageFilter< ScalarImageType, ScalarImageType, ScalarImageType  > MaxFilterType;
      MaxFilterType::Pointer maxFilter = MaxFilterType::New();
      maxFilter->SetInput1(duplicator->GetOutput());
      maxFilter->SetInput2(readImage);
      maxFilter->Update();
      maxAccumulatorImage = maxFilter->GetOutput();
    }

  }

  typedef itk::ShiftScaleImageFilter<ScalarImageType, ScalarImageType> DivideFilterType;
  DivideFilterType::Pointer divider = DivideFilterType::New();
  divider->SetInput(meanAccumulatorImage);
  double scaleFactor=1.0/numberOfImages;
  divider->SetScale(scaleFactor);
  divider->SetShift(0);
  divider->Update();

  if (!meanImageFilename.empty())
  {

    typedef itk::ImageFileWriter< ScalarImageType  > meanWriterType;
    meanWriterType::Pointer meanWriter = meanWriterType::New();
    meanWriter->SetFileName( meanImageFilename );
    meanWriter->SetInput(divider->GetOutput());
    try
    {
      meanWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (!maxImageFilename.empty())
  {

    typedef itk::ImageFileWriter< ScalarImageType  > maxWriterType;
    maxWriterType::Pointer maxWriter = maxWriterType::New();
    maxWriter->SetFileName( maxImageFilename );
    maxWriter->SetInput(maxAccumulatorImage);
    try
    {
      maxWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (!meanHistogramFilename.empty())
  {
    WriteHistogram(meanHistogramFilename.c_str(), divider->GetOutput(), histogramBinCount, histogramMin, histogramMax);
  }

  if (!maxHistogramFilename.empty())
  {
    WriteHistogram(maxHistogramFilename.c_str(), maxAccumulatorImage, histogramBinCount, histogramMin, histogramMax);
  }

  return EXIT_SUCCESS;
}
