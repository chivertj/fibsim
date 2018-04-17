#ifndef ITKHELPERS_HEADER
#define ITKHELPERS_HEADER

#include "itkImage.h"
#include "itkImageFileWriter.h"

typedef double EntropyDataType;
typedef itk::Image<EntropyDataType,3> EntropyImageType;
typedef EntropyImageType OPImageT;

typedef unsigned short  PixelType;
const unsigned int OutputImageDimension = 3;
typedef itk::Image< PixelType, OutputImageDimension >   MidImageType;

void CreateImage(const MidImageType::SizeType &size,  const MidImageType::SpacingType &spacing, OPImageT::Pointer image);

template <class IMAGETYPE>
void CreateImage(const typename IMAGETYPE::SizeType &size,  const typename IMAGETYPE::SpacingType &spacing, typename IMAGETYPE::Pointer image) {
  typename IMAGETYPE::RegionType region;
  typename IMAGETYPE::IndexType start;
  start[0]=0; start[1]=0; start[2]=0;
  
  region.SetSize(size);
  region.SetIndex(start);
  
  //  image=OPImageT::New();
  image->SetSpacing(spacing);
  image->SetRegions(region);
  image->Allocate();
}


template <class ImageType> void WriteImage(const typename ImageType::Pointer image, const std::string &fileName) {
  typedef  typename itk::ImageFileWriter< ImageType  > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName);
  writer->SetInput(image);
  writer->Update();
}


#include "itkExtractImageFilter.h"
typedef unsigned char ExtractPixelType;
typedef itk::Image<ExtractPixelType, 3> ImageExtract3DType;
typedef itk::Image<ExtractPixelType, 2> ImageExtract2DType;
#include <vector>
typedef std::vector<ImageExtract2DType::Pointer> ImageExtract2DVectorType;

ImageExtract2DVectorType Extract2D(ImageExtract3DType::Pointer ipimage3D);

ImageExtract3DType::Pointer Stack2DInto3D(const ImageExtract2DVectorType &images2D, double zspacing);

template <class T> std::vector<typename T::PixelType> ComputeSliceMeans(typename T::Pointer ipimage3D, typename T::PixelType threshold) {
  typename T::RegionType region=ipimage3D->GetLargestPossibleRegion();
  int X=region.GetSize()[0],Y=region.GetSize()[1],Z=region.GetSize()[2];
  std::vector<typename T::PixelType> slicemeans(Z);
  typename T::IndexType pxlidx;
  int count;
  for (int z=0;z<Z;z++) {
    slicemeans[z]=0;
    count=0;
    pxlidx[2]=z;
    for (int y=0;y<Y;y++) {
      pxlidx[1]=y;
      for (int x=0;x<X;x++) {
	pxlidx[0]=x;
	if (ipimage3D->GetPixel(pxlidx)>threshold) {
	  slicemeans[z]+=ipimage3D->GetPixel(pxlidx);
	  count++;
	}
      }
    }
    if (count>0)
      slicemeans[z]/=count;
    else slicemeans[z]=threshold;
  }
  return slicemeans;
}

template <class T> std::vector<typename T::PixelType> ComputeSliceMaxs(typename T::Pointer ipimage3D) {
  typename T::RegionType region=ipimage3D->GetLargestPossibleRegion();
  int X=region.GetSize()[0],Y=region.GetSize()[1],Z=region.GetSize()[2];
  std::vector<typename T::PixelType> slicemaxs(Z);
  typename T::IndexType pxlidx;
  for (int z=0;z<Z;z++) {
    slicemaxs[z]=0;
    pxlidx[2]=z;
    for (int y=0;y<Y;y++) {
      pxlidx[1]=y;
      for (int x=0;x<X;x++) {
	pxlidx[0]=x;
	if (slicemaxs[z]<ipimage3D->GetPixel(pxlidx))
	  slicemaxs[z]=ipimage3D->GetPixel(pxlidx);
      }
    }
//    slicemaxs[z]/=(X*Y);
  }
  return slicemaxs;
}

template <class T> std::vector<typename T::PixelType> ComputeSliceMins(typename T::Pointer ipimage3D) {
  typename T::RegionType region=ipimage3D->GetLargestPossibleRegion();
  int X=region.GetSize()[0],Y=region.GetSize()[1],Z=region.GetSize()[2];
  std::vector<typename T::PixelType> slicemins(Z);
  typename T::IndexType pxlidx;
  for (int z=0;z<Z;z++) {
    slicemins[z]=1e100;
    pxlidx[2]=z;
    for (int y=0;y<Y;y++) {
      pxlidx[1]=y;
      for (int x=0;x<X;x++) {
	pxlidx[0]=x;
	if (slicemins[z]>ipimage3D->GetPixel(pxlidx) && ipimage3D->GetPixel(pxlidx)>vnl_math::eps)
	  slicemins[z]=ipimage3D->GetPixel(pxlidx);
      }
    }
//    slicemins[z]/=(X*Y);
  }
  return slicemins;
}

#endif
