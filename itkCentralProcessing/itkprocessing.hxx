#ifndef ITKPROCESSING_HEADER
#define ITKPROCESSING_HEADER

#include "itkhelpers.hxx"
#include "itkImage.h"
#include "itkanglehists.hxx"
#include "itkBinaryThresholdImageFilter.h"

MidImageType::Pointer ThresholdVolume(MidImageType::Pointer origimage,int lowerThreshold, int upperThreshold);

#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

typedef itk::BinaryImageToLabelMapFilter<MidImageType> BinaryImageToLabelMapFilterType;

BinaryImageToLabelMapFilterType::Pointer LabelComponents(MidImageType::Pointer image, int minimumobjectsize);
BinaryImageToLabelMapFilterType::Pointer LabelComponents(MidImageType::Pointer image, int minimumobjectsize, int maximumobjectsize);

#include "itkmathops.hxx"

#include "vectorfileread.hxx"
#include <string>

class ObjectsSamples {
public:
  FRC::VOL_T objectsvectors;
  int Nobjects;
  void convert(const BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter, const double *dims);
  void read(const std::string &filename);
};

class ObjectOrientations {
public:
  ObjectOrientations(void) : Nobjects(-1),fibrethickness(1) {}
  static const unsigned int MeasurementVectorLength=3;
  typedef itk::Vector<double,MeasurementVectorLength> MeasurementVectorType;
  typedef itk::Statistics::ListSample<MeasurementVectorType> SampleType;
  typedef itk::Statistics::CovarianceSampleFilter< SampleType > CovarianceAlgorithmType;

  void  DetermineOrientations(const ObjectsSamples &objectssamples);
  FRC::VOL_T GetEndPoints(void);
  void DetermineOrientationsFromEndPnts(const FRC::VOL_T &endpoints);
  typedef std::vector<PointType> VecPoints;
  void WritePoints(const VecPoints &pntstowrite, const string &filename);
  VecPoints means;
  VecPoints aers;
  VecPoints eigvals;

  EntropyImageType::Pointer entropyimage;
  EntropyDataType DetermineEntropy(const MidImageType::Pointer origimage, const int bgThreshold, const double *dims, const float radialscale, const float distancescale );
  EntropyDataType DetermineMultiScaleEntropy(const MidImageType::Pointer origimage, const int bgThreshold, const double *dims);
  EntropyDataType DetermineMultiScaleEntropy(void); //not spatial
  EntropyDataType DetermineMultiScaleSOP(void); //not spatial, Scalar Order Parameter
  EntropyDataType DetermineSingleScaleSOP(void); //not spatial, Scalar Order Parameter
  EntropyDataType DetermineSOPfromEndPnts(const FRC::VOL_T &endpoints);
  void DetermineEffFactorsfromEndPnts(const FRC::VOL_T &endpoints, std::vector<EntropyDataType> &efffacts);
  int Nobjects;

  VecPoints x1points; //starting points of each fibre
  VecPoints x2points; //ending points of each fibre
  void SetImageSizes(PointType _imagesize) {imagesize=_imagesize;}
  PointType imagesize;

  VecPoints x1_2points; //orthogonal start points (2nd eigenvector / value)
  VecPoints x2_2points; //orthogonal end points
  VecPoints x1_3points; //orthogonal start points (3rd eigenvector / value)
  VecPoints x2_3points; //orthogonal end points

  void SetFibreThickness(float _fibrethickness) {fibrethickness=_fibrethickness;}
  float fibrethickness;

  //determines if a line between points x1 and x2 
  //is well overlapped by the thresholded volume
  float DetermineLineProportion(const MidImageType::Pointer origimage, const PointType _x1, const PointType _x2, const PixelType threshold);
  static const float INOUTPROPORTION_THRESHOLD=0.4;
  void PruneLineSegments(const MidImageType::Pointer origimage, PixelType threshold);

  static const float EIGVALRATIO=5.;
  void PruneLineSegmentsEccentricity(void);
};

//template <class LABELIMGT, class IMG_OUTT, typename PIXEL_OUTT> void LabelMapToImg(typename itk::BinaryImageToLabelMapFilter<LABELIMGT>::Pointer labels, const PIXEL_OUTT &pixval, typename IMG_OUTT::pointer &imgout) {
template <typename labelT, typename pixelT, class imgT> void LabelMapToImg(typename labelT::Pointer labels, pixelT val, typename imgT::Pointer imgout)
 {
    for(unsigned int i = 0; i < labels->GetOutput()->GetNumberOfLabelObjects(); i++) {
      // Get the ith region
      BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labels->GetOutput()->GetNthLabelObject(i);
      // Output the pixels composing the region
      for(unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) 
	imgout->SetPixel(labelObject->GetIndex(pixelId),val);
    }
}

template <typename imgT> MidImageType::Pointer ThresholdVolume(typename imgT::Pointer origimage, int lowerThreshold, int upperThreshold) {
  //  std::cout <<lowerThreshold<<","<<upperThreshold<<std::endl;
  typedef typename itk::BinaryThresholdImageFilter<imgT,MidImageType> BinaryThresholdImageFilterType;
  typename BinaryThresholdImageFilterType::Pointer thresholdFilter=BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(origimage);
  thresholdFilter->SetLowerThreshold(lowerThreshold);
  thresholdFilter->SetUpperThreshold(upperThreshold);
  thresholdFilter->SetInsideValue(65535);
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->Update();
  MidImageType::Pointer image=thresholdFilter->GetOutput();
  image->Update();
  return image;
}

#endif
