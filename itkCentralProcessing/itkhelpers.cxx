#include "itkhelpers.hxx"
#include <vector>

#include "itkTileImageFilter.h"
#include "itkChangeInformationImageFilter.h"


//template <class OPImageT, class MidImageType> 
void CreateImage(const MidImageType::SizeType &size,  const MidImageType::SpacingType &spacing, OPImageT::Pointer image) {
  OPImageT::RegionType region;
  OPImageT::IndexType start;
  start[0]=0; start[1]=0; start[2]=0;
  
  region.SetSize(size);
  region.SetIndex(start);
  
  //  image=OPImageT::New();
  image->SetSpacing(spacing);
  image->SetRegions(region);
  image->Allocate();
}
                 

ImageExtract2DVectorType Extract2D(ImageExtract3DType::Pointer ipimage3D) {
  typedef itk::ExtractImageFilter<ImageExtract3DType,ImageExtract2DType> ExtractImageFilterType;


  ImageExtract3DType::RegionType ipregion=ipimage3D->GetLargestPossibleRegion();
  ImageExtract3DType::SizeType size=ipregion.GetSize();
  unsigned int N=size[2];
  size[2]=0;

  ImageExtract3DType::IndexType start=ipregion.GetIndex();
  ImageExtract3DType::RegionType desiredRegion;
  desiredRegion.SetSize(size);

  ImageExtract2DVectorType opvector;
  for (int i=0;i<N;i++) {
    ExtractImageFilterType::Pointer extract2D=ExtractImageFilterType::New();
    extract2D->InPlaceOn();
    extract2D->SetDirectionCollapseToSubmatrix();
    start[2]=i;
    desiredRegion.SetIndex(start);
    extract2D->SetExtractionRegion(desiredRegion);
    extract2D->SetInput(ipimage3D);
    extract2D->Update();
    opvector.push_back(extract2D->GetOutput());
  }
  return opvector;
}

ImageExtract3DType::Pointer Stack2DInto3D(const ImageExtract2DVectorType &images2D, double zspacing) {
  assert(images2D.size()>1);
  const unsigned int OutputImageDimension = 3;
  typedef itk::TileImageFilter<ImageExtract2DType,ImageExtract3DType> TilerType;
  TilerType::Pointer tiler=TilerType::New();
  itk::FixedArray<unsigned int, OutputImageDimension> layout;
  layout[0]=1; layout[1]=1; layout[2]=0;
  tiler->SetLayout(layout);
  for (int i=0;i<images2D.size();i++) {
    images2D[i]->DisconnectPipeline();
    tiler->SetInput(i,images2D[i]);
  }
  PixelType filler=0;
  tiler->SetDefaultPixelValue(filler);
  tiler->Update();
  ImageExtract3DType::Pointer image3D=tiler->GetOutput();
  image3D->Update();

  ImageExtract3DType::SpacingType spacing=image3D->GetSpacing();
  spacing[2]=zspacing;
  typedef itk::ChangeInformationImageFilter<ImageExtract3DType> ChangeInformationType;
  ChangeInformationType::Pointer changeInformation=ChangeInformationType::New();
  changeInformation->SetInput(image3D);
  changeInformation->SetOutputSpacing(spacing);
  changeInformation->ChangeSpacingOn();
  image3D=changeInformation->GetOutput();
  image3D->Update();
  return image3D;
}


