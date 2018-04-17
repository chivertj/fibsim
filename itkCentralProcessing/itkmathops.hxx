#ifndef ITKMATHOPS_HEADER
#define ITKMATHOPS_HEADER

#include "itkVector.h"
#include "itkListSample.h"
#include "itkMeanSampleFilter.h"
#include "itkCovarianceSampleFilter.h"

#include "itkSymmetricEigenAnalysis.h"
#include "vnl/vnl_matrix.h"
#include "itkFixedArray.h"
#include "itkMatrix.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkhelpers.hxx"
#include "itk1Dhist.hxx"

typedef vnl_matrix< double > InputMatrixType;
typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
typedef itk::SymmetricEigenAnalysis< InputMatrixType,  
				     EigenValuesArrayType, EigenVectorMatrixType > SymmetricEigenAnalysisType;

class EigenMaker {
public:
  EigenValuesArrayType eigenvalues;
  EigenVectorMatrixType eigenvectors;
  SymmetricEigenAnalysisType symmetricEigenSystem;

  EigenMaker(void) : symmetricEigenSystem(3) {}
  void Make(InputMatrixType m) {
    symmetricEigenSystem.ComputeEigenValuesAndVectors(m, eigenvalues, eigenvectors);
  }
};

#include "vectorfileread.hxx"
  typedef  vnl_vector_fixed<double,3> VEC_T;
  typedef  vnl_matrix_fixed<double,3,3> MAT_T;
typedef std::vector<VEC_T> StdVnlVecs;
const StdVnlVecs DetermineNorms(const FRC::VOL_T &_endpoints);
  //////////////////////
//const vnl_vector_fixed<double,3>& frcVecToVNL(const FRC::PNT_T &v) {
const VEC_T frcVecToVNL(const FRC::PNT_T &v);

////////////////////////////////////////////////////////////////////
//calculation from Suuronen et al. J. Mater Sci 2013 (48):1358-1367
class SOPfromAlignMatrix {
public:
  VEC_T result;
  FRC::VOL_T endpoints;
  StdVnlVecs normalvecs;
  /////////////////////
  //calculation from Suuronen et al. J. Mater Sci 2013 (48):1358-1367
  MAT_T alignm;
  const MAT_T& CalcAlignMatrix(void) {
    //    vnl_matrix_fixed<double,3,3> alignm;
    alignm=0;
    for (int i=0;i<normalvecs.size();i++) 
      alignm+=outer_product(normalvecs[i],normalvecs[i]);
    alignm/=normalvecs.size();
    //    std::cout <<alignm[0][0]<<","<<alignm[1][1]<<","<<alignm[2][2]<<std::endl;
    //this makes the SOP values more noisy - is that correct?! traceless-symmetric matrix should be
    //    alignm[0][0]=0; alignm[1][1]=0; alignm[2][2]=0;
    return alignm;
  }
  ///////////////////////
  //
  double sopval;
  double CalcSOPVal(void) {
    EigenMaker em;
    em.Make(alignm);
    sopval=em.eigenvalues[2];
    //return sopval*3/2;
    return sopval;
  }
  ////////////////////////////
  double CalcAll(const FRC::VOL_T &_endpoints) {
    endpoints=_endpoints;
    normalvecs=DetermineNorms(_endpoints);    
#if(1)
    CalcAlignMatrix();
    return CalcSOPVal();
#else
    return CalcFiberAppSOPVals();
#endif
  }
  double CalcFiberAppSOPVals(void) {
    double A=0,B=0,N=normalvecs.size();
    for (int i=0;i<N;i++) {
      A+=pow(normalvecs[i][0],2);
      B+=normalvecs[i][0]*normalvecs[i][1]*normalvecs[i][2];
    }
    return sqrt(pow(2*A-N,2)+4*pow(B,2))/N*1.5-0.5;
  }
};

////////////////////////////////////////////////////////////////////
//calculation from Vincente et al. NDTE 2014.
class OrientEfficiencyFactors {
public:
  VEC_T result;
  FRC::VOL_T endpoints;
  StdVnlVecs normalvecs;
  StdVnlVecs angles;
  static const int nbins=30;
  histQx1D<double,3,double,nbins> anghists;
  std::vector<double> efffactors;
  //////////////////////
  OrientEfficiencyFactors(void) {
    anghists.initialize(0,itk::Math::pi);
  }
  ///////////////////////////
  const StdVnlVecs CalcAngles(const StdVnlVecs &_normalvecs) {
    int Nobjects=_normalvecs.size();
    assert(Nobjects>0);
    StdVnlVecs _angles=StdVnlVecs(Nobjects);
    for (int i=0;i<Nobjects;i++) {
      for (int j=0;j<3;j++) {
	_angles[i][j]=acos(_normalvecs[i][j]);
	if (_angles[i][j]<0)
	  _angles[i][j]+=itk::Math::pi; //so that they lie between 0 and pi. Later (in eff fact calcs) will need to adjust to between -pi/2 to pi/2
      }
    }
    return _angles;
  }
  ////////////////////////////
  void CalcAngHists(const StdVnlVecs &_angles) {
    for (int i=0;i<_angles.size();i++) 
      anghists.add(_angles[i]);
    anghists.normalize();
  }
  ////////////////////////////
  void CalcEffFactors(void) {
    efffactors=std::vector<double>(3);
    efffactors[0]=0.; efffactors[1]=0.; efffactors[2]=0.; 
    for (int i=0;i<nbins;i++) {
      for (int j=0;j<3;j++) 
	efffactors[j]+=2*anghists[j][i]*cos(anghists[j].midbinvalue(i)-itk::Math::pi/2.);
    }
  }
  ////////////////////////////
#if(0)
  void CalcAll(const FRC::VOL_T &_endpoints, std::vector<EntropyDataType> &_efffactors) {
    endpoints=_endpoints;
    normalvecs=DetermineNorms(endpoints);
    angles=CalcAngles(normalvecs);
    CalcAngHists(angles);
    CalcEffFactors();
    _efffactors=efffactors;
  }
#else
  void CalcAll(const FRC::VOL_T &_endpoints, std::vector<EntropyDataType> &_efffactors) {
    endpoints=_endpoints;
    normalvecs=DetermineNorms(endpoints);
    efffactors=std::vector<double>(3);
    efffactors[0]=0.; efffactors[1]=0.; efffactors[2]=0.; 
    for (int i=0;i<normalvecs.size();i++) {
      for (int j=0;j<3;j++)
	efffactors[j]+=2*abs(normalvecs[i][j]);
    }
    for (int j=0;j<3;j++)
      efffactors[j]/=normalvecs.size();
    _efffactors=efffactors;
  }
#endif
};


#endif //ITKMATHOPS_HEADER
