#ifndef ITKLINEPROCESSING_HEADER
#define ITKLINEPROCESSING_HEADER

#include "itkhelpers.hxx"
#include "itkLineIterator.h"


typedef itk::Point<double,3> PointType3D;
typedef itk::Point<double,2> PointType2D;
//////////////////////////////////////////////////////////////////////////
//from https://github.com/ComplexSystemsModeling/DCMTK-ITK/blob/master/Modules/ThirdParty/VNL/src/vxl/core/vnl/vnl_cross.h
template<class T>
inline vnl_vector_fixed<T,3> vnl_cross_3d( const vnl_vector_fixed<T,3>& v1, const vnl_vector_fixed<T,3>& v2 )
{
  vnl_vector_fixed<T,3> result;
  result[0] = v1[1] * v2[2] - v1[2] * v2[1]; // work for both col/row
  result[1] = v1[2] * v2[0] - v1[0] * v2[2]; // representation
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}
//////////////////////////////////////////////////////////////////////////
template<class T> inline T vnl_cross_2d(vnl_vector_fixed<T,2> const& v1, vnl_vector_fixed<T,2> const& v2) {
  return v1[0] * v2[1] - v1[1] * v2[0];
}
//////////////////////////////////////////////////////////////////////////
//from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
template <class T> T DetermineDistance3D(const itk::Point<T,3> &x0, const itk::Point<T,3> &x1, const itk::Point<T,3> &x2) {
  itk::Point<T,3> A=x0-x1,B=x0-x2;
  T denominator=x2.EuclideanDistanceTo(x1);
  vnl_vector_fixed<T,3> vnl_A, vnl_B;
  for (int i=0;i<3;i++) {
    vnl_A[i]=A[i];
    vnl_B[i]=B[i];
  }
  vnl_vector_fixed<T,3> vnl_numvec=vnl_cross_3d<T>(vnl_A,vnl_B);
  return vnl_numvec.two_norm()/denominator;
}
//////////////////////////////////////////////////////////////////////////
//for line segments
//////////////////////////////////////////////////////////////////////////
template <class T> T DetermineDistance3Dseg(const itk::Point<T,3> &x0, const itk::Point<T,3> &x1, const itk::Point<T,3> &x2) {
  T d=DetermineDistance3D(x0,x1,x2);
  T x1x2d=x2.EuclideanDistanceTo(x1);
  T maxh=sqrt(pow(x1x2d,2)+pow(d,2));

  T x0x1d=x0.EuclideanDistanceTo(x1);
  T x0x2d=x0.EuclideanDistanceTo(x2);

  if (x0x1d>maxh || x0x2d>maxh)  //must be beyond line segment extent
    return vnl_math_min(x0x1d,x0x2d);    
  return d;
}
//for 2D
template <class T> T DetermineDistance2D(const itk::Point<T,2> &x0, const itk::Point<T,2> &x1, const itk::Point<T,2> &x2) {
  //  itk::Point<T,2> A=x0-x1,B=x0-x2;
  itk::Point<T,2> A=x0-x1,B=x0-x2;
  T denominator=x2.EuclideanDistanceTo(x1);
  vnl_vector_fixed<T,2> vnl_A, vnl_B;
  for (int i=0;i<2;i++) {
    vnl_A[i]=A[i];
    vnl_B[i]=B[i];
  }
  T numerator=abs(vnl_cross_2d<T>(vnl_A,vnl_B));
  return numerator/denominator;
}
////////////////////////////////////////////////////////////////
typedef float ExtIndexType;
typedef vnl_vector_fixed<ExtIndexType,3> ExtPointType;
typedef std::vector<ExtPointType> ExtIndices;
////////////////////////////////////////////////////////////////
void Profiles(MidImageType::Pointer image, const ExtIndices &x1points, const ExtIndices &x2points);
////////////////////////////////////////////////////////////////
std::string& replaceAll(std::string& context, const std::string& from, const std::string& to);
////////////////////////////////////////////////////////////////
template <class T> void CSVSize(const std::string &filename, const std::string &seperator, int &ncols, int &nrows) {
  std::ifstream csvfile(filename.c_str());
  std::string line;
  std::getline(csvfile,line);
  if (seperator!=" ")
    replaceAll(line,seperator," ");
  std::stringstream ss(line);
  T dataelement;
  ncols=0;
  while (ss>>dataelement)
    ncols++;
  nrows=1;
  while (std::getline(csvfile,line)) {
    if (line!="")
      nrows++;
  }
}
////////////////////////////////////////////////////////////////
template <class T> void ReadCSVtoLinesSegments(const std::string &filename, const std::string &seperator, ExtIndices &x1points, ExtIndices &x2points ) {
  int ncols,nrows;
  CSVSize<T>(filename,seperator,ncols,nrows);
  if (ncols!=6) {std::cout <<"CSV File not reading correctly."<<std::endl; exit(0);}
  x1points.resize(nrows);
  x2points.resize(nrows);

  std::ifstream csvfile(filename.c_str());
  std::string line;

  T vals[6];
  for (int j=0;j<nrows;j++) {
    std::getline(csvfile,line);
    if (seperator!=" ")
      replaceAll(line,seperator," ");
    std::stringstream ss(line);
    for (int i=0;i<ncols; i++)
      ss >> vals[i];
    for (int i=0;i<3;i++) {
      x1points[j][i]=vals[i];
      x2points[j][i]=vals[i+3];
    }
  }
}
////////////////////////////////////////////////////////////////
void DetermineMinMaxMid(MidImageType::Pointer image, MidImageType::IndexType x1, MidImageType::IndexType x2, MidImageType::SpacingType inputSpacing, float &min, float &max, float &mid, int &discretelength);
////////////////////////////////////////////////////////////////
double TraceLength(MidImageType::Pointer image, const MidImageType::IndexType x1, const MidImageType::IndexType x2, const MidImageType::SpacingType inputSpacing, const float min, const float max, const float mid, const int discretelength);
////////////////////////////////////////////////////////////////
void CalcRadiuses(MidImageType::Pointer image, const ExtIndices &x1points, const ExtIndices &x2points);
////////////////////////////////////////////////////////////////
#endif
