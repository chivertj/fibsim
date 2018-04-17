#include <string>

#include "itkImageFileReader.h"

#include "../itkCentralProcessing/itkhelpers.hxx"
#include "../itkCentralProcessing/vectorfileread.hxx"
#include "../itkCentralProcessing/itkprocessing.hxx"

//#define FLOATIMAGET
//possible compile time definitions: FLOATIMAGET, WITH_AERS
#ifdef FLOATIMAGET
typedef float PixelTypeL;
#else
typedef PixelType PixelTypeL;
#endif
typedef itk::Image<PixelTypeL,OutputImageDimension> ImageTypeL;
typedef itk::ImageFileReader< ImageTypeL > ImageReaderType;

int main(int argc, char *argv[] ) {
  if (argc < 5) {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " \n\t[1] input \n\t[2] lower threshold \n\t[3] upper threshold \n\t[4] minimum object size \n\t[5] maximum object size (optional)\n" << std::endl;
    return EXIT_FAILURE;
  }
 
  std::string ipfilename(argv[1]);
  int lowerthreshold=atoi(argv[2]),upperthreshold=atoi(argv[3]);
  int minobjectsize=atoi(argv[4]);

  int maxobjectsize=10000;
  if (argc==6)
    maxobjectsize=atoi(argv[5]);

  ImageReaderType::Pointer reader = ImageReaderType::New();
 
  reader->SetFileName( ipfilename );
  reader->UpdateLargestPossibleRegion();

  ImageTypeL::Pointer origimage=reader->GetOutput();

  //  MidImageType::Pointer thresholdedimage=ThresholdVolume(origimage,lowerthreshold,upperthreshold);
  MidImageType::Pointer thresholdedimage=ThresholdVolume<ImageTypeL>(origimage,lowerthreshold,upperthreshold);
  WriteImage<MidImageType>(thresholdedimage,"thresholdedimage.mhd");

  BinaryImageToLabelMapFilterType::Pointer labelledcomponents=LabelComponents(thresholdedimage,minobjectsize,maxobjectsize);

  //template <typename LABELIMGT, typename IMG_OUTT, typename PIXEL_OUTT> void LabelMapToImg(itk::BinaryImageToLabelMapFilter<LABELIMGT>::Pointer labels, const PIXEL_OUTT &pixval, typename IMG_OUTT::pointer &imgout) {
  //  LabelMapToImg<BinaryImageToLabelMapFilterType,PixelTypeL,ImageTypeL>(labelledcomponents,PixelTypeL(65535),origimage);
  //  WriteImage<ImageTypeL>(origimage,"labelledorigimage.mhd");
  thresholdedimage->FillBuffer( itk::NumericTraits< PixelType >::Zero );
  LabelMapToImg<BinaryImageToLabelMapFilterType,PixelType,MidImageType>(labelledcomponents,PixelTypeL(65535),thresholdedimage);
  WriteImage<MidImageType>(thresholdedimage,"labelledorigimage.mhd");

  double dims[3];
  PointType imagesize;
  for (int i=0;i<3;i++) {
    dims[i]=origimage->GetSpacing()[i];
    imagesize[i]=origimage->GetLargestPossibleRegion().GetSize()[i];
  }

  ObjectsSamples objectssamples;
  objectssamples.convert(labelledcomponents,dims);

  ObjectOrientations objectorientations;
  objectorientations.SetImageSizes(imagesize);

  //  objectorientations.SetFibreThickness(fibrethickness);
  objectorientations.DetermineOrientations(objectssamples);
  //  std::cout <<"Nobjects:"<<objectorientations.Nobjects<<std::endl;
  //  objectorientations.PruneLineSegmentsEccentricity();
  //  std::cout <<"Nobjects:"<<objectorientations.Nobjects<<std::endl;
  //  EntropyDataType wholevolumeentropy=objectorientations.DetermineEntropy(origimage,bgThreshold,dims);

#ifdef WITH_AERS
  //output aers for all accepted objected...
  for (int i=0;i<objectorientations.Nobjects;i++) 
    std::cout <<objectorientations.aers[i][0]<<","<<objectorientations.aers[i][1]<<","<<objectorientations.aers[i][2]<<std::endl;
#else
#ifdef WITH_LINEPROFILES 
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.x1points[i][0]<<","<<objectorientations.x1points[i][1]<<","<<objectorientations.x1points[i][2]<<","<<objectorientations.x2points[i][0]<<","<<objectorientations.x2points[i][1]<<","<<objectorientations.x2points[i][2]<<std::endl;

  std::cout <<std::endl;

  //2nd eigen vector
  for (int i=0;i<objectorientations.Nobjects;i++)  {
    std::cout <<objectorientations.x1_2points[i][0]<<","<<objectorientations.x1_2points[i][1]<<","<<objectorientations.x1_2points[i][2]<<",";
    std::cout <<objectorientations.x2_2points[i][0]<<","<<objectorientations.x2_2points[i][1]<<","<<objectorientations.x2_2points[i][2]<<std::endl;
  }
  std::cout <<std::endl;
  //3rd eigen vector
  for (int i=0;i<objectorientations.Nobjects;i++)   {
    std::cout <<objectorientations.x1_3points[i][0]<<","<<objectorientations.x1_3points[i][1]<<","<<objectorientations.x1_3points[i][2]<<",";
    std::cout <<objectorientations.x2_3points[i][0]<<","<<objectorientations.x2_3points[i][1]<<","<<objectorientations.x2_3points[i][2]<<std::endl;
  }

#else
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.x1points[i][0]<<","<<objectorientations.x1points[i][1]<<","<<objectorientations.x1points[i][2]<<","<<objectorientations.x2points[i][0]<<","<<objectorientations.x2points[i][1]<<","<<objectorientations.x2points[i][2]<<std::endl;

  std::cout <<std::endl;

  //2nd eigen vector: -ve E2 lambda2
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.x1_2points[i][0]<<","<<objectorientations.x1_2points[i][1]<<","<<objectorientations.x1_2points[i][2]<<","
	      <<objectorientations.means[i][0]<<","<<objectorientations.means[i][1]<<","<<objectorientations.means[i][2]<<std::endl;
  std::cout <<std::endl;
  //2nd eigen vector: +ve E2 lambda2
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.means[i][0]<<","<<objectorientations.means[i][1]<<","<<objectorientations.means[i][2]<<","
	      <<objectorientations.x2_2points[i][0]<<","<<objectorientations.x2_2points[i][1]<<","<<objectorientations.x2_2points[i][2]<<std::endl;
  std::cout <<std::endl;

  //3rd eigen vector: -ve E3 lambda3
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.x1_3points[i][0]<<","<<objectorientations.x1_3points[i][1]<<","<<objectorientations.x1_3points[i][2]<<","
	      <<objectorientations.means[i][0]<<","<<objectorientations.means[i][1]<<","<<objectorientations.means[i][2]<<std::endl;
  std::cout <<std::endl;
  //2nd eigen vector: +ve E2 lambda2
  for (int i=0;i<objectorientations.Nobjects;i++)  
    std::cout <<objectorientations.means[i][0]<<","<<objectorientations.means[i][1]<<","<<objectorientations.means[i][2]<<","
	      <<objectorientations.x2_3points[i][0]<<","<<objectorientations.x2_3points[i][1]<<","<<objectorientations.x2_3points[i][2]<<std::endl;


#endif  
#endif
  

#if(0)
  EntropyDataType wholevolumeentropy=objectorientations.DetermineMultiScaleEntropy(origimage,0,dims);  

  WriteImage<EntropyImageType>(objectorientations.entropyimage,argv[2]);
  //  std::cout <<wholevolumeentropy<<std::endl;

  //  std::cout <<std::endl;
  typedef std::vector<EntropyImageType::PixelType> SLICESTATS_T;
  SLICESTATS_T slicemeans=ComputeSliceMeans<EntropyImageType>(objectorientations.entropyimage,vnl_math::eps);
  SLICESTATS_T slicemaxs=ComputeSliceMaxs<EntropyImageType>(objectorientations.entropyimage);
  SLICESTATS_T slicemins=ComputeSliceMins<EntropyImageType>(objectorientations.entropyimage);
  std::cout <<std::accumulate(slicemeans.begin(),slicemeans.end(),0.0)/slicemeans.size()<<","<<*std::min_element(slicemins.begin(),slicemins.end())<<","<<*std::max_element(slicemaxs.begin(),slicemaxs.end())<<std::endl<<std::endl;
  for (int z=0;z<slicemeans.size();z++) 
    std::cout <<slicemeans[z]<<","<<slicemins[z]<<","<<slicemaxs[z]<<std::endl;
#endif  
  return EXIT_SUCCESS;
}
