#include "itkanglehists.hxx"

unsigned int DiscretizeAER(const PointType &aer) {
  double factA=aer[0]-fmod(aer[0],angleDelta);
  double factE=aer[1]-fmod(aer[1],angleDelta);
  double opbinA=factA/angleDelta+abs(minA/angleDelta);
  double opbinE=factE/angleDelta+abs(minE/angleDelta);
  return opbinA+opbinE*(maxA-minA)/angleDelta;
}

unsigned int DiscretizeAER(const PointType &aer, const double _angleDelta) {
  double factA=aer[0]-fmod(aer[0],_angleDelta);
  double factE=aer[1]-fmod(aer[1],_angleDelta);
  double opbinA=factA/_angleDelta+abs(minA/_angleDelta);
  double opbinE=factE/_angleDelta+abs(minE/_angleDelta);
  return opbinA+opbinE*(maxA-minA)/_angleDelta;
}

unsigned int ComputeHistSize(const double _angleDelta) {
  return int((maxE-minE)/_angleDelta*(maxA-minA)/_angleDelta);
}

//returns center of range of angles covered by index bin
#if(0)
const PointType UnDiscretizeAER(unsigned int index, const double _angleDelta) {
  PointType aer;
  //index=opbinA+opbinE*(maxA-minA)/_angleDelta

#if(0)
  double Asize=(maxA-minA)*_angleDelta;
  double opbinE=fmod(index,Asize);
#else
  double Asize=(maxA-minA);
  double opbinE=(maxA-minA)
fmod(index,_angleDelta);
#endif
  double opbinA=(index-opbinE)/Asize;
  aer[0]=opbinA; //*_angleDelta+_angleDelta/2;
  aer[1]=opbinE; //*_angleDelta+_angleDelta/2;
  aer[2]=0;

  return aer;
}
#else
const PointType UnDiscretizeAER(unsigned int index, const double _angleDelta) {
  double totalangle=index*_angleDelta;
  double Eangle=fmod(totalangle,maxA-minA);
  double Atotalangle=totalangle-Eangle;
  double Aangle=Atotalangle/(maxA-minA)*_angleDelta;
  PointType aer;
  aer[0]=Aangle+_angleDelta/2; aer[1]=Eangle+_angleDelta/2; aer[2]=0;
  return aer;
}
#endif
