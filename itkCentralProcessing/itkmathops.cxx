#include "itkmathops.hxx"
///////////////////////////////////////////////////////////////////
  const StdVnlVecs DetermineNorms(const FRC::VOL_T &_endpoints) {
    int Nobjects=_endpoints.size();
    assert(Nobjects>0);
    StdVnlVecs _normalvecs=StdVnlVecs(_endpoints.size());
    double norm;
    for(unsigned int i = 0; i < Nobjects; i++) {
      norm=0;
      for (int j=0;j<3;j++) {
	//	means[i][j]=(endpoints[i][0][j]+endpoints[i][1][j])/2.;
	_normalvecs[i][j]=_endpoints[i][0][j]-_endpoints[i][1][j];
	norm+=pow(_normalvecs[i][j],2);
      }
      for (int j=0;j<3;j++)
	_normalvecs[i][j]/=sqrt(norm);
    }
    return _normalvecs;
  }
///////////////////////////////////////////////////////////////////
const VEC_T frcVecToVNL(const FRC::PNT_T &v) {
  assert(v.size()==3);
  VEC_T result;
  for (int i=0;i<3;i++)
    result[i]=v[i];
  return result;
}
///////////////////////////////////////////////////////////////////
