#ifndef ITKANGLEHISTS_HEADER
#define ITKANGLEHISTS_HEADER

#include <vector>

#include "itkPoint.h"
typedef itk::Point<double,3> PointType;
const double minA=0;
const double maxA=180;
const double minE=0;
const double maxE=180;
const double angleDelta=15;//30 //45; //60
const unsigned int AERHistSize = int((maxE-minE)/angleDelta*(maxA-minA)/angleDelta);

unsigned int DiscretizeAER(const PointType &aer);
unsigned int DiscretizeAER(const PointType &aer, const double _angleDelta);
const PointType UnDiscretizeAER(unsigned int index, const double _angleDelta);

unsigned int ComputeHistSize(const double _angleDelta);


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
template <typename T> class MSHists {
  static const int N=4; //the number of elements in the following
  static const double msAngleDeltas[N];
  static const int AERmsHistSizes[N]; 
  ///
  //hists is an array of 1 dimensional histograms
  //a 1 dimensional histogram has bins accessed through DiscretizeAER
  typedef std::vector<T> histT;
  std::vector< histT > hists;
  std::vector< histT > probs;
  ///
  histT sums;
  histT entropies;
  histT factors;
  T ttlentropy;
  T meanfactor;
public:
  ///
  MSHists (void):hists(N),sums(N),probs(N),entropies(N),factors(N) { 
    for (int i=0;i<N;i++) {
      //      std::cout <<"Hist Sizes:"<<AERmsHistSizes[i]<<std::endl;
      hists[i]=histT(AERmsHistSizes[i]); 
      probs[i]=histT(AERmsHistSizes[i]); 
    }
    Reset();
  }
  ///
  void AddToBin(const PointType &aer) {
    int histbin;
    for (int i=0;i<N;i++) {
      histbin=DiscretizeAER(aer,msAngleDeltas[i]);
      //      std::cout <<histbin<<"/"<<hists[i].size()<<std::endl;
      //      if (histbin>=hists[i].size()) std::cout <<std::string("*************************\n histbin larger than hist[i].size()\n******************************")<<std::endl;
      //      std::cout <<hists.size()<<" "<<histbin<<"/"<<hists[i].size()<<std::endl;
#if(1)
      hists[i][histbin]+=1;
#endif
    }
  }
  ///
  histT& operator [](int i) { return hists[i]; }
  ///
  const T GetFactor(int i) const { return factors[i];}
  ///
  const T GetEntropy(int i) const { return entropies[i];}
  ///
  T GetSum(int i) {return sums[i];}
  ///
  T GetMeanSum(void) {T mean=0;for (int i=0;i<N;i++) mean+=sums[i]; return mean/N; }
  ///
  void Reset(void) {
    ttlentropy=0;
    meanfactor=0;
    for (int i=0;i<N;i++) {
      for (int j=0;j<hists[i].size();j++) {
	hists[i][j]=0;    
	probs[i][j]=0;
      }
      sums[i]=0;
    }
  }
  ///
  static const double GetDelta(unsigned int i) { return msAngleDeltas[i]; }
  ///
  static const int Size(void) { return N; }
  ///
  void StoreSumAndNormalize(void) {
    if (sums[0]>vnl_math::eps) //already normalised
      return;
    for (int i=0;i<N;i++) 
      sums[i]=0;
    for (int i=0;i<N;i++) {
      for (int j=0;j<hists[i].size();j++) 
	sums[i]+=hists[i][j];
    }
    for (int i=0;i<N;i++) {
      if (sums[i]>vnl_math::eps) {
	for (int j=0;j<hists[i].size();j++) 
	  probs[i][j]=hists[i][j]/sums[i];
      }
      else {
	for (int j=0;j<hists[i].size();j++) 
	  probs[i][j]=0;
      }
    }
  }
  ///
  T ComputeAllEntropyAllSteps(void) {
    StoreSumAndNormalize();
    ttlentropy=0;
    for (int i=0;i<N;i++) {
      entropies[i]=ComputeEntropy(i);
      ttlentropy+=entropies[i];
    }
    return ttlentropy;
  }
  ///
  T ComputeAllFactorAllSteps(void) {
    StoreSumAndNormalize();
    meanfactor=0;
    for (int i=0;i<N;i++) {
      factors[i]=ComputeFactor(i);
      meanfactor+=factors[i];
    }
    meanfactor/=N;
    return meanfactor;
  }
  ///
private:
  T ComputeEntropy(int i) {
    T entropyval=0;
    for (int j=0;j<probs[i].size();j++)  {
      if (probs[i][j]>vnl_math::eps)
	entropyval+=probs[i][j]*log2(1./probs[i][j]);
    }
    return entropyval;
  }
  T ComputeFactor(int i) {
#if(1)
    //compute mean: theta_hat = sum theta*p(theta); phi_hat = sum phi*p(phi)
    PointType aer_hat,aer;
    aer_hat[0]=0; aer_hat[1]=0; aer_hat[2]=0;
    for (int j=0;j<probs[i].size();j++) {
      //      std::cout <<probs[i][j]<<",";
      //      probs[i][j]=1/probs[i].size();
      aer=UnDiscretizeAER(j,msAngleDeltas[i]);	
      if (probs[i][j]>vnl_math::eps) {
	aer_hat[0]+=aer[0]*probs[i][j];
	aer_hat[1]+=aer[1]*probs[i][j];
	//	std::cout <<probs[i][j]<<",";
      }
      //      std::cout <<aer[0]<<","<<aer[1]<<" ";
    }
    //compute differences: theta_alpha = theta - theta_hat; phi_alpha = phi - phi_hat
    PointType aer_alpha,aer_moment;
    aer_moment[0]=0; aer_moment[1]=0; aer_moment[2]=0;
    //compute weighted moment: sum cos(theta_alpha)*p(theta); sum cos(phi_alpha)*p(alpha)
    for (int j=0;j<probs[i].size();j++) {
      if (probs[i][j]>vnl_math::eps) {
	aer=UnDiscretizeAER(j,msAngleDeltas[i]);
	aer_alpha[0]=aer[0]-aer_hat[0];
	aer_alpha[1]=aer[1]-aer_hat[1];
	aer_moment[0]+=(1.5*pow(cos(aer_alpha[0]*0.017453293),2)-0.5)*probs[i][j];
	aer_moment[1]+=(1.5*pow(cos(aer_alpha[1]*0.017453293),2)-0.5)*probs[i][j];
      }
    }
    //    return aer_moment[1];
    //    return aer_hat[0];
    //    std::cout <<"hat:"<<aer_hat[0]<<","<<aer_hat[1]<<"||m:"<<aer_moment[0]<<","<<aer_moment[1]<<" ";
    //std::cout <<aer_moment[0]<<" ";
    return aer_moment[1];
#else
    T factorval=0;
    for (int j=0;j<probs[i].size();j++) {
      //      std::cout <<msAngleDeltas[i]/2<<","<<msAngleDeltas[i]*j<<"->"<<msAngleDeltas[i]/2+msAngleDeltas[i]*j<<std::endl;

      if (probs[i][j]>vnl_math::eps) {
	PointType aer=UnDiscretizeAER(j,msAngleDeltas[i]);
	factorval+=probs[i][j]*(3*pow(cos(msAngleDeltas[i]/2+msAngleDeltas[i]*j),2)-1)/2;
      }
    }
    return factorval;
#endif
  }
};
//////////
//template<typename T> const double MSHists<T>::msAngleDeltas[N]={1,5,15,30,45,60};
//template<typename T> const double MSHists<T>::msAngleDeltas[N]={30,36,45,60};
#if(1)
template<typename T> const double MSHists<T>::msAngleDeltas[N]={15,30,45,60};
#else
template<typename T> const double MSHists<T>::msAngleDeltas[N]={15};
#endif
//template<typename T> const double MSHists<T>::msAngleDeltas[N]={1};
//template<typename T> const double MSHists<T>::msAngleDeltas[N]={5};
//template<typename T> const int MSHists<T>::AERmsHistSizes[N]= {
//  int((maxE-minE)/msAngleDeltas[0]*(maxA-minA)/msAngleDeltas[0]) };
//  int((maxE-minE)/msAngleDeltas[0]*(maxA-minA)/msAngleDeltas[0]),
  //  int((maxE-minE)/msAngleDeltas[1]*(maxA-minA)/msAngleDeltas[1]),
//  int((maxE-minE)/msAngleDeltas[2]*(maxA-minA)/msAngleDeltas[2]),
//  int((maxE-minE)/msAngleDeltas[3]*(maxA-minA)/msAngleDeltas[3]),
//  int((maxE-minE)/msAngleDeltas[3]*(maxA-minA)/msAngleDeltas[4]),
//  int((maxE-minE)/msAngleDeltas[3]*(maxA-minA)/msAngleDeltas[5])
//};
#if(1)
template<typename T> const int MSHists<T>::AERmsHistSizes[N]= {
  int((maxE-minE)/msAngleDeltas[0]*(maxA-minA)/msAngleDeltas[0]),
  int((maxE-minE)/msAngleDeltas[1]*(maxA-minA)/msAngleDeltas[1]),
  int((maxE-minE)/msAngleDeltas[2]*(maxA-minA)/msAngleDeltas[2]),
  int((maxE-minE)/msAngleDeltas[3]*(maxA-minA)/msAngleDeltas[3])
};
#else
template<typename T> const int MSHists<T>::AERmsHistSizes[N]= {
  int((maxE-minE)/msAngleDeltas[0]*(maxA-minA)/msAngleDeltas[0])};
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//a container for multiple MSHists. Each MSHists is also a multidimensional
//histogram for multiscale in the radial dimension.
//This is multiscale in terms of variance (distance).
template <typename T> class MS2DHists {
  static const int M=4;
  static const double invvars[M];
  typedef MSHists<T>  MSHistT;
  typedef std::vector<MSHistT> MSHistsT;
  MSHistsT mshists;
  T ttlentropy;
public:
  MS2DHists(void) : mshists(M) {
  }  
  void Reset(void) {
    for (int j=0;j<M;j++)
      mshists[j].Reset();
  }
  MSHistT& operator[] (int j) { return mshists[j]; }
  static const int Size(void) { return M; }
  void StoreSumAndNormalize(void) {
    for (int j=0;j<M;j++)
      mshists[j].StoreSumAndNormalize();
  }
  T GetMeanSum(void) {
    T meansum=0;
    for (int j=0;j<M;j++)
      meansum+=mshists[j].GetMeanSum();
    return meansum/M;
  }
  T ComputeAllEntropyAllSteps(void) {
    ttlentropy=0;
    for (int j=0;j<M;j++)
      ttlentropy+=mshists[j].ComputeAllEntropyAllSteps();
    return ttlentropy;
  }
  void AddToBin(const PointType &aer, const T distance, const T varscale) {
    std::vector<int> histbins(MSHistT::Size());
    //    std::cout <<"aer:"<<aer[0]<<","<<aer[1]<<","<<aer[2]<<" ";
    //    std::cout <<"histbin values:";
    for (int i=0;i<MSHistT::Size();i++) {
      histbins[i]=DiscretizeAER(aer,MSHistT::GetDelta(i));
      //      std::cout <<"<"<<MSHistT::GetDelta(i)<<" "<<histbins[i]<<",";
    }
    //    std::cout <<" hist sizes:";
    T distsq=pow(distance,2);
    for (int j=0;j<M;j++) {
      T binval=exp(-distsq*invvars[j]/varscale);
      for (int i=0;i<MSHistT::Size();i++) {
#if(1)
	mshists[j][i][histbins[i]]+=binval;
#else
      std::cout <<mshists[j][i].size()<<",";
	mshists[j][i][0]+=binval;
#endif
      }
    }
	//    std::cout <<std::endl;
  }
};
/////////
  template<typename T> const double MS2DHists<T>::invvars[M]={1./2,1./(2*100),1./(2*10000),1./(2*1000000)};
//  template<typename T> const double MS2DHists<T>::invvars[M]={1./(2)};
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

#endif //ITKANGLEHISTS_HEADER
