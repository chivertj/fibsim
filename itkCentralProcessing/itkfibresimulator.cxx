#include "itkfibresimulator.hxx"
#include <itkGaussianDistribution.h>

#include "itkLineIterator.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"


#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkFlatStructuringElement.h"

#include "itkMath.h"

#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkCastImageFilter.h"


//////////////////////////////////////////////////////////////////////////////
void  SetBackgroundValue(MidImageType::Pointer highresimage,unsigned short pixelval) {
  std::cout <<"Setting (mean) background values..."<<std::endl;
  MidImageType::RegionType region=highresimage->GetLargestPossibleRegion();
  itk::ImageRegionIterator<MidImageType> imageIt(highresimage,region);
  PixelType backgroundpixelval=pixelval;
  while(!imageIt.IsAtEnd()) {
    imageIt.Set(backgroundpixelval);
    ++imageIt;
  }
}
//////////////////////////////////////////////////////////////////////////////
FRC::VOL_T AddFibres(MidImageType::Pointer highresimage,MidImageType::RegionType::SizeType superdims,MidImageType::SpacingType super_spacing,unsigned short fibrepixelvalue,int Nfibres,float fibrelengths, int *maxdeltaangles, int *Nanglequantizations) {
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();

  std::cout <<"Adding fibres...#"<<Nfibres<<std::endl;
  MidImageType::IndexType src,dest;
  PixelType fibrepixelval=fibrepixelvalue;
  float initU=generator->GetIntegerVariate(Nanglequantizations[0]);
  float initV=generator->GetIntegerVariate(Nanglequantizations[1]);
  float theta,phi;
  float U,V;

  FRC::PNT_T frc_pnt(3);
  FRC::VOL_T frc_allfibres(Nfibres);
  int fibreidx=0;
  for (int f=0;f<Nfibres;f++) {
    do {
      for (int i=0;i<3;i++)
	src[i]=generator->GetIntegerVariate(superdims[i]-1);
      do {
	U=generator->GetIntegerVariate(maxdeltaangles[0]);
	U+=-maxdeltaangles[0]/2+initU;
      }while (U<0 || U>=Nanglequantizations[0]);
      U/=Nanglequantizations[0];
      do {
	V=generator->GetIntegerVariate(maxdeltaangles[1]);
	V+=-maxdeltaangles[1]/2+initV;
      } while(V<0 || V>=Nanglequantizations[1]);
      V/=Nanglequantizations[1];
      //sphere sampling as described by wolfram, http://mathworld.wolfram.com/SpherePointPicking.html
      theta=2*itk::Math::pi*U;
      phi=acos(2*V-1);
      //transformation according to: http://mathworld.wolfram.com/SphericalCoordinates.html
      dest[0]=int(fibrelengths/super_spacing[0]*cos(theta)*sin(phi)+0.5)+src[0];
      dest[1]=int(fibrelengths/super_spacing[1]*sin(theta)*sin(phi)+0.5)+src[1];
      dest[2]=int(fibrelengths/super_spacing[2]*cos(phi)+0.5)+src[2];
    } while (dest[0]<0 || dest[0]>=superdims[0] || dest[1]<0 || dest[1]>=superdims[1] || dest[2]<0 || dest[2]>=superdims[0]);
    itk::LineIterator<MidImageType> lineIt(highresimage,src,dest);
    FRC::FIBRE_T frc_fibre;
    lineIt.GoToBegin();
    //lineIt.GetIndex()[i] etc
    while (!lineIt.IsAtEnd()) {
      lineIt.Set(fibrepixelval);
      for (int i=0;i<3;i++)
	frc_pnt[i]=lineIt.GetIndex()[i]*super_spacing[i];
      frc_fibre.push_back(frc_pnt);
      ++lineIt;
    }
    frc_allfibres[fibreidx++]=frc_fibre;
  }
  return frc_allfibres;
}
//////////////////////////////////////////////////////////////////////////////
// this one is for ground truth generation as it does not add them to a volume 
FRC::VOL_T AddFibres(MidImageType::RegionType::SizeType superdims,MidImageType::SpacingType super_spacing,int Nfibres,float fibrelengths, int *maxdeltaangles, int *Nanglequantizations) {
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();

  //  std::cout <<"Adding fibres...#"<<Nfibres<<std::endl;
  FRC::PNT_T src(3),dest(3);
  float initU=generator->GetIntegerVariate(Nanglequantizations[0]);
  float initV=generator->GetIntegerVariate(Nanglequantizations[1]);
  float theta,phi;
  float U,V;
  FRC::PNT_T frc_pnt(3);
  FRC::VOL_T frc_allfibres(Nfibres);

  FRC::FIBRE_T frc_angles(Nfibres);
  for (int i=0;i<Nfibres;i++)
    frc_angles[i]=FRC::PNT_T(2);
#if(0)
  for (int i=0;i<3;i++) 
    src[i]=superdims[i]/2;
#endif
  for (int f=0;f<Nfibres;f++) {
    do {
#if(1)
      for (int i=0;i<3;i++) 
	src[i]=generator->GetIntegerVariate(superdims[i]-1);
#endif
      do {
	U=generator->GetIntegerVariate(maxdeltaangles[0]);
	U+=-maxdeltaangles[0]/2+initU;
      }while (U<0 || U>=Nanglequantizations[0]);

      U/=Nanglequantizations[0];
      //      U=fmod(U,1);
      do {
	V=generator->GetIntegerVariate(maxdeltaangles[1]);
	V+=-maxdeltaangles[1]/2+initV;
      } while(V<0 || V>=Nanglequantizations[1]);

      V/=Nanglequantizations[1];
      V=fmod(V,1);
      //sphere sampling as described by wolfram, http://mathworld.wolfram.com/SpherePointPicking.html
      theta=2*itk::Math::pi*U;
      phi=acos(2*V-1);

      frc_angles[f][0]=theta*57.296;
      frc_angles[f][1]=phi*57.296;

      //transformation according to: http://mathworld.wolfram.com/SphericalCoordinates.html
      dest[0]=int(fibrelengths/super_spacing[0]*cos(theta)*sin(phi)+0.5)+src[0];
      dest[1]=int(fibrelengths/super_spacing[1]*sin(theta)*sin(phi)+0.5)+src[1];
      dest[2]=int(fibrelengths/super_spacing[2]*cos(phi)+0.5)+src[2];
    } while (dest[0]<0 || dest[0]>=superdims[0] || dest[1]<0 || dest[1]>=superdims[1] || dest[2]<0 || dest[2]>=superdims[0]);
    //    itk::LineIterator<MidImageType> lineIt(highresimage,src,dest);
    FRC::FIBRE_T frc_fibre;
    frc_fibre.push_back(src);
    frc_fibre.push_back(dest);
    frc_allfibres[f]=frc_fibre;
  }
  FRC::vectorfilewrite("aes.csv",frc_angles);
  return frc_allfibres;
}
//////////////////////////////////////////////////////////////////////////////
MidImageType::Pointer DilateFibres(MidImageType::Pointer highresimage,unsigned short radiusval,unsigned short fibrepixelvalue) {
  typedef itk::BinaryBallStructuringElement<PixelType, 3> StructuringElementType;
  //  typedef itk::BinaryCrossStructuringElement<PixelType, 3> StructuringElementType;
  unsigned int radiusvalue=radiusval;
  StructuringElementType::RadiusType radius;
  radius.Fill(radiusvalue);

  StructuringElementType structuringElement;  
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <MidImageType, MidImageType, StructuringElementType> BinaryDilateImageFilter3DType;  
  BinaryDilateImageFilter3DType::Pointer dilateFilter=BinaryDilateImageFilter3DType::New();
  dilateFilter->SetDilateValue(fibrepixelvalue);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->SetInput(highresimage);
  dilateFilter->Update();
  highresimage=dilateFilter->GetOutput();
  return highresimage;
}
////////////////////////////////////////////////////////////////////////////////
void SetNoise(MidImageType::Pointer highresimage, double noiseparamval) {
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();
  std::cout <<"Setting noise values..."<<std::endl;
  MidImageType::RegionType region2=highresimage->GetLargestPossibleRegion();
  MidImageType::PixelType pixelval;
  itk::ImageRegionIterator<MidImageType> imageIt2(highresimage,region2);
  //  imageIt.GoToBegin();
  itk::Statistics::GaussianDistribution::Pointer gaussian = itk::Statistics::GaussianDistribution::New();
  while(!imageIt2.IsAtEnd()) {
    pixelval=imageIt2.Value();
    //    pixelval+=generator->GetIntegerVariate(backgroundpixelvalue)-noiseparamval;
    pixelval+=gaussian->EvaluateInverseCDF(generator->GetUniformVariate(0,1),0,noiseparamval);
    imageIt2.Set(pixelval);
    ++imageIt2;
  }
}
////////////////////////////////////////////////////////////////////////////////
MidImageType::Pointer DownSampleVolume(MidImageType::Pointer highresimage, float *lowressize, MidImageType::SpacingType highres_spacing, float upsamplerate, float interplanedownsamplerate) {

  std::cout <<"Downsampling..."<<std::endl;

  const     unsigned int    Dimension = 3;
  typedef   float           InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension >   InternalImageType;

  typedef itk::CastImageFilter<MidImageType,InternalImageType> CastFilterType;
  CastFilterType::Pointer castFilter=CastFilterType::New();
  castFilter->SetInput(highresimage);
  InternalImageType::Pointer highresimage_f=castFilter->GetOutput();

  typedef itk::RecursiveGaussianImageFilter<InternalImageType,
					    InternalImageType > GaussianFilterType;
  GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
  GaussianFilterType::Pointer smootherY = GaussianFilterType::New();
  GaussianFilterType::Pointer smootherZ = GaussianFilterType::New();

  smootherX->SetSigma(lowressize[0]);
  smootherY->SetSigma(lowressize[1]);
  smootherZ->SetSigma(lowressize[2]);

  smootherX->SetDirection(0);
  smootherY->SetDirection(1);
  smootherZ->SetDirection(2);

  smootherX->SetOrder(GaussianFilterType::ZeroOrder);
  smootherY->SetOrder(GaussianFilterType::ZeroOrder);
  smootherZ->SetOrder(GaussianFilterType::ZeroOrder);

  smootherX->SetNormalizeAcrossScale(false);
  smootherY->SetNormalizeAcrossScale(false);
  smootherZ->SetNormalizeAcrossScale(false);

  smootherX->SetInput( highresimage_f );
  smootherY->SetInput( smootherX->GetOutput() );
  smootherZ->SetInput( smootherY->GetOutput() );

  std::cout <<"Running smoothers..."<<std::endl;
  smootherZ->Update();

  typedef itk::CastImageFilter<InternalImageType,MidImageType> CastFilterOpType;
  CastFilterOpType::Pointer castOpFilter=CastFilterOpType::New();
  castOpFilter->SetInput(smootherZ->GetOutput());
  castOpFilter->Update();
  //  return castOpFilter->GetOutput();

  MidImageType::Pointer subsampledimage=SubSampler(castOpFilter->GetOutput(),highres_spacing,upsamplerate,interplanedownsamplerate);
  return subsampledimage;
}
////////////////////////////////////////////////////////////////////////////////
int ExportVolume(const MidImageType::Pointer highresimage,const std::string &opfilename) {
  typedef itk::ImageFileWriter< MidImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(highresimage);
  writer->SetFileName( opfilename.c_str() );
 
  try {
    writer->Update();
  }
  catch( itk::ExceptionObject & excp ) {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
////////////////////////////////////////////////////////////////////////////////
MidImageType::Pointer SubSampler(const MidImageType::Pointer filteredimage, MidImageType::SpacingType highres_spacing, float upsamplerate, float interplanedownsamplerate) {
  MidImageType::RegionType filteredregion=filteredimage->GetLargestPossibleRegion();
  MidImageType::SpacingType inputSpacing=filteredimage->GetSpacing();
  MidImageType::SizeType ipSize=filteredregion.GetSize();

  MidImageType::SpacingType lowres_spacing;
  MidImageType::RegionType::SizeType opSize;
  typedef MidImageType::RegionType::SizeType::SizeValueType SizeValueType;
  for (int i=0;i<2;i++) {
    lowres_spacing[i]=highres_spacing[i]*upsamplerate;
    opSize[i] = static_cast<SizeValueType>( ipSize[i]/upsamplerate );
  }
  lowres_spacing[2]=highres_spacing[2]*interplanedownsamplerate;
  opSize[2] = static_cast<SizeValueType>( ipSize[2]/interplanedownsamplerate );

  MidImageType::RegionType subsampledregion;
  subsampledregion.SetSize(opSize);
  MidImageType::IndexType start;
  start[0]=0; start[1]=0; start[2]=0;
  subsampledregion.SetIndex(start);

  MidImageType::Pointer subsampledimage=MidImageType::New();
  subsampledimage->SetRegions(subsampledregion);
  subsampledimage->SetSpacing(lowres_spacing);
  subsampledimage->Allocate();

  MidImageType::IndexType hridx,lridx;

  for (lridx[2]=0;lridx[2]<opSize[2];lridx[2]++) {
    if (lridx[2]>0)
      hridx[2]=lridx[2]*interplanedownsamplerate;
    else
      hridx[2]=5; //smoothers seem to not filter first slice
    for (lridx[1]=0;lridx[1]<opSize[1];lridx[1]++) {
      hridx[1]=lridx[1]*upsamplerate;
      for (lridx[0]=0;lridx[0]<opSize[0];lridx[0]++) {
	hridx[0]=lridx[0]*upsamplerate;
	subsampledimage->SetPixel(lridx,filteredimage->GetPixel(hridx));
      }
    }
  }

  return subsampledimage;
}
////////////////////////////////////////////////////////////////////////////////
