#include "itklineprocessing.hxx"

void Profiles(MidImageType::Pointer image, const ExtIndices &x1points, const ExtIndices &x2points) {
  MidImageType::PixelType intensity;
  MidImageType::RegionType region=image->GetLargestPossibleRegion();

  MidImageType::IndexType x1,x2;
  MidImageType::SpacingType inputSpacing=image->GetSpacing();
  MidImageType::SizeType imagesize=region.GetSize();
  for (int i=0;i<x1points.size();i++) {
    for (int j=0;j<3;j++) {
      x1[j]=int(x1points[i][j]/inputSpacing[j]+0.5);
      x2[j]=int(x2points[i][j]/inputSpacing[j]+0.5);
      if (x1[j]<0)
	x1[j]=0;
      else if (x1[j]>(imagesize[j]-1))
	x1[j]=imagesize[j]-1;
      if (x2[j]<0)
	x2[j]=0;
      else if (x2[j]>(imagesize[j]-1))
	x2[j]=imagesize[j]-1;
    }
    itk::LineIterator<MidImageType> it1(image,x1,x2);
    it1.GoToBegin();
    while (!it1.IsAtEnd()) {
      intensity=it1.Value();
      std::cout <<intensity<<",";
      ++it1;
    }
    std::cout <<std::endl;
  }
}
////////////////////////////////////////////////////////////////
std::string& replaceAll(std::string& context, const std::string& from, const std::string& to) {
  size_t lookHere = 0;
  size_t foundHere;
  while((foundHere = context.find(from, lookHere))
      != std::string::npos) {
    context.replace(foundHere, from.size(), to);
    lookHere = foundHere + to.size();
  }
  return context;
}
////////////////////////////////////////////////////////////////
void DetermineMinMaxMid(MidImageType::Pointer image, MidImageType::IndexType x1, MidImageType::IndexType x2, MidImageType::SpacingType inputSpacing, float &min, float &max, float &mid, int &discretelength) {
  MidImageType::PixelType intensity;
  itk::LineIterator<MidImageType> it1(image,x1,x2);
  it1.GoToBegin();
  min=65535; max=0;
  discretelength=0;
  while (!it1.IsAtEnd()) {
    intensity=it1.Value();
    if (intensity<min)
      min=intensity;
    else if (intensity>max)
      max=intensity;
    ++it1;
    discretelength++;
  }
  mid=(max+min)/2.;
}
////////////////////////////////////////////////////////////////
double TraceLength(MidImageType::Pointer image, const MidImageType::IndexType x1, const MidImageType::IndexType x2, const MidImageType::SpacingType inputSpacing, const float min, const float max, const float mid, const int discretelength) {
  MidImageType::PixelType intensity;
  itk::LineIterator<MidImageType> it1(image,x1,x2);
  it1.GoToBegin();
  int lowcount=0;
  int distance=0;
  while (!it1.IsAtEnd()) {
    intensity=it1.Value();
    if (intensity<mid) {
      if (lowcount==0) 
	lowcount++;
      else if (lowcount==1)
	break;
    }
    else if (lowcount==1) //then reset as false alarm
      lowcount=0;
    ++it1;
    distance++;
  }
  float linefraction=distance/float(discretelength);

  itk::Point<double,3> p1,p2,p_end;
  for (int i=0;i<3;i++) {
    p1[i]=x1[i];
    p1[i]*=inputSpacing[i];
    p2[i]=x2[i];
    p2[i]*=inputSpacing[i];
    p_end[i]=(p2[i]-p1[i])*linefraction+p1[i];
  }
  return p1.EuclideanDistanceTo(p_end);
}
////////////////////////////////////////////////////////////////
void CalcRadiuses(MidImageType::Pointer image, const ExtIndices &x1points, const ExtIndices &x2points) {
  MidImageType::PixelType intensity;
  MidImageType::RegionType region=image->GetLargestPossibleRegion();

  MidImageType::IndexType x1,x2;
  MidImageType::SpacingType inputSpacing=image->GetSpacing();
  MidImageType::SizeType imagesize=region.GetSize();
  float min,max,mid;
  int discretelength;
  double continuouslength;
  for (int i=0;i<x1points.size();i++) {
    for (int j=0;j<3;j++) {
      x1[j]=int(x1points[i][j]/inputSpacing[j]+0.5);
      x2[j]=int(x2points[i][j]/inputSpacing[j]+0.5);
      if (x1[j]<0)
	x1[j]=0;
      else if (x1[j]>(imagesize[j]-1))
	x1[j]=imagesize[j]-1;
      if (x2[j]<0)
	x2[j]=0;
      else if (x2[j]>(imagesize[j]-1))
	x2[j]=imagesize[j]-1;
    }
    DetermineMinMaxMid(image,x1,x2,inputSpacing,min,max,mid,discretelength);
    continuouslength=TraceLength(image,x1,x2,inputSpacing,min,max,mid,discretelength);
    std::cout <<discretelength<<","<<continuouslength<<std::endl;
  }
}
////////////////////////////////////////////////////////////////
