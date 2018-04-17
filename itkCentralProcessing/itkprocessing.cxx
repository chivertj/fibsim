#include "itkprocessing.hxx"

#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"

#include "vectorfileread.hxx"
#include "itklineprocessing.hxx"
#include "itkLineIterator.h"
#include "itkanglehists.hxx"
#include <fstream>

////////////////////////////////////////////////////////////////////////////
MidImageType::Pointer ThresholdVolume(MidImageType::Pointer origimage,int lowerThreshold, int upperThreshold) {
  typedef itk::BinaryThresholdImageFilter<MidImageType,MidImageType> BinaryThresholdImageFilterType;
  BinaryThresholdImageFilterType::Pointer thresholdFilter=BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(origimage);
  thresholdFilter->SetLowerThreshold(lowerThreshold);
  thresholdFilter->SetUpperThreshold(upperThreshold);
  thresholdFilter->SetInsideValue(65535);
  thresholdFilter->SetOutsideValue(0);
  MidImageType::Pointer image=thresholdFilter->GetOutput();
  image->Update();
  return image;
}

////////////////////////////////////////////////////////////////////////////
BinaryImageToLabelMapFilterType::Pointer LabelComponents(MidImageType::Pointer image, int minimumobjectsize) {
  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter=BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput(image);
  binaryImageToLabelMapFilter->Update();
  //  std::cout <<"no objects:"<<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()<<std::endl;
  //  std::cout <<"(minimum object size:"<<minimumobjectsize<<")"<<std::endl;
  std::vector<unsigned long> labelsToRemove;
  for (unsigned int i=0;i<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();i++) {
    BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType *labelObject=binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    if (labelObject->Size()<minimumobjectsize) 
      labelsToRemove.push_back(labelObject->GetLabel());
  }

  for (unsigned int i=0;i<labelsToRemove.size();i++) 
    binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
  int Nobjects=binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();

  binaryImageToLabelMapFilter->Update();
  //  std::cout <<"revised no objects:"<<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()<<std::endl;
  return binaryImageToLabelMapFilter;
}

////////////////////////////////////////////////////////////////////////////
BinaryImageToLabelMapFilterType::Pointer LabelComponents(MidImageType::Pointer image, int minimumobjectsize, int maximumobjectsize) {
  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter=BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput(image);
  binaryImageToLabelMapFilter->Update();
  //  std::cout <<"no objects:"<<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()<<std::endl;
  //  std::cout <<"(minimum object size:"<<minimumobjectsize<<")"<<std::endl;
  std::vector<unsigned long> labelsToRemove;
  for (unsigned int i=0;i<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();i++) {
    BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType *labelObject=binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    if (labelObject->Size()<minimumobjectsize || labelObject->Size()>maximumobjectsize)
      labelsToRemove.push_back(labelObject->GetLabel());
  }

  for (unsigned int i=0;i<labelsToRemove.size();i++) 
    binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
  int Nobjects=binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();

  binaryImageToLabelMapFilter->Update();
  //  std::cout <<"revised no objects:"<<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()<<std::endl;
  return binaryImageToLabelMapFilter;
}

///////////////////////////////////////////////////////
// version that uses FRC::VOL_T from (initially) ground
// truth end points, i.e. a pair of points - so no 
// estimation required.
///////////////////////////////////////////////////////
void ObjectOrientations::DetermineOrientationsFromEndPnts(const FRC::VOL_T &endpoints){
  assert(endpoints.size()>0);
  Nobjects=endpoints.size();

  means=VecPoints(Nobjects);
  aers=VecPoints(Nobjects);
  eigvals=VecPoints(Nobjects);

  x1points=VecPoints(Nobjects);
  x2points=VecPoints(Nobjects);

  PointType vec;
  double norm;
  for(unsigned int i = 0; i < Nobjects; i++) {
    norm=0;
    for (int j=0;j<MeasurementVectorLength;j++) {
      means[i][j]=(endpoints[i][0][j]+endpoints[i][1][j])/2.;
      vec[j]=endpoints[i][0][j]-endpoints[i][1][j];
      norm+=pow(vec[j],2);
    }
    for (int j=0;j<MeasurementVectorLength;j++)
      vec[j]/=sqrt(norm);

    aers[i][0]=atan2(vec[1],vec[0])/(2*vnl_math::pi)*360;
    aers[i][1]=acos(vec[2])/(2*vnl_math::pi)*360;

#if(0)
    aers[i][0]+=180;
    aers[i][1]+=180;
    aers[i][0]=fmod(aers[i][0],180);
    aers[i][1]=fmod(aers[i][1],180);
#endif
    aers[i][2]=sqrt(pow(vec[2],2)+pow(vec[1],2)+pow(vec[0],2));

    for (int el=0;el<3;el++) {
      x1points[i][el]=endpoints[i][0][el];
      x2points[i][el]=endpoints[i][1][el];
    }
  }
}

///////////////////////////////////////////////////////
// version that uses ObjectsSamples
///////////////////////////////////////////////////////
void  ObjectOrientations::DetermineOrientations(const ObjectsSamples &objectssamples) {
  assert(objectssamples.Nobjects>0);
  Nobjects=objectssamples.Nobjects;
  SampleType::Pointer sample=SampleType::New();
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = CovarianceAlgorithmType::New();

  MeasurementVectorType mv;

  means=VecPoints(Nobjects);
  aers=VecPoints(Nobjects);
  eigvals=VecPoints(Nobjects);

  x1points=VecPoints(Nobjects);
  x2points=VecPoints(Nobjects);

  x1_2points=VecPoints(Nobjects);
  x2_2points=VecPoints(Nobjects);
  x1_3points=VecPoints(Nobjects);
  x2_3points=VecPoints(Nobjects);

  PointType mins;
  PointType maxs;
  for(unsigned int i = 0; i < Nobjects; i++) {
    SampleType::Pointer sample =SampleType::New();
    for (int l=0;l<MeasurementVectorLength;l++) {
      mins[l]=9999999;
      maxs[l]=0;
      means[i][l]=0;
    }
    //////////////////////
    // Get the ith region and extract out the samples that will be used for determining the sample vector
    int objectsize=objectssamples.objectsvectors[i].size();
    for(unsigned int pixelId = 0; pixelId < objectsize; pixelId++) {
      for (int l=0;l<MeasurementVectorLength;l++) {
	mv[l]=objectssamples.objectsvectors[i][pixelId][l];
	if (mv[l]<mins[l])
	  mins[l]=mv[l];
	if (mv[l]>maxs[l])
	  maxs[l]=mv[l];
      }
      sample->PushBack(mv);
    }
    //////////////////////
    covarianceAlgorithm->SetInput(sample);
    covarianceAlgorithm->Update();

    CovarianceAlgorithmType::MatrixType covarmat=covarianceAlgorithm->GetCovarianceMatrix();
    means[i]=covarianceAlgorithm->GetMean();
    EigenMaker eigenmaker;
    eigenmaker.Make(covarmat.GetVnlMatrix());

    PointType vec;
    for (int el=0;el<3;el++) {
      vec[el]=eigenmaker.eigenvectors(2,el);
      eigvals[i][el]=eigenmaker.eigenvalues[el];
    }
    aers[i][0]=atan2(vec[1],vec[0])/(2*vnl_math::pi)*360;
    aers[i][1]=acos(vec[2])/(2*vnl_math::pi)*360;

#if(0)
    aers[i][0]+=180;
    aers[i][1]+=180;
#endif
    //    std::cout <<"aers:"<<aers[i][0]<<","<<aers[i][1]<<"->";
    aers[i][0]=fmod(aers[i][0],180);
    aers[i][1]=fmod(aers[i][1],180);
    if (aers[i][0]<0)
      aers[i][0]+=180;
    if (aers[i][1]<0)
      aers[i][1]+=180;
    //    std::cout <<aers[i][0]<<","<<aers[i][1]<<std::endl;
    double normvals[3]={sqrt(eigenmaker.eigenvalues[2])*2,sqrt(eigenmaker.eigenvalues[1])*5,sqrt(eigenmaker.eigenvalues[0])*5};
    aers[i][2]=normvals[0];

    for (int el=0;el<3;el++) {
      x1points[i][el]=-vec[el]*normvals[0]+means[i][el];
      x2points[i][el]=vec[el]*normvals[0]+means[i][el];
    }
    for (int el=0;el<3;el++) {
      x1_2points[i][el]=-eigenmaker.eigenvectors(1,el)*normvals[1]+means[i][el];
      x2_2points[i][el]=eigenmaker.eigenvectors(1,el)*normvals[1]+means[i][el];
      x1_3points[i][el]=-eigenmaker.eigenvectors(0,el)*normvals[2]+means[i][el];
      x2_3points[i][el]=eigenmaker.eigenvectors(0,el)*normvals[2]+means[i][el];
    }
  }
}
////////////////////////////////////////////////////////////////////////////
FRC::VOL_T ObjectOrientations::GetEndPoints(void) {
  assert(Nobjects>0);
  FRC::VOL_T endpoints(Nobjects);
  for (int i=0;i<Nobjects;i++) {
    endpoints[i]=FRC::FIBRE_T(2);
    endpoints[i][0]=FRC::PNT_T(3);
    for (int j=0;j<3;j++) 
      endpoints[i][0][j]=x1points[i][j];
    endpoints[i][1]=FRC::PNT_T(3);
    for (int j=0;j<3;j++) 
      endpoints[i][1][j]=x2points[i][j];
  }
  return endpoints;
}
////////////////////////////////////////////////////////////////////////////
EntropyDataType ObjectOrientations::DetermineEntropy(const MidImageType::Pointer origimage, const int bgThreshold, const double *dims, const float radialscale, const float distancescale) {
  assert(Nobjects>0);

  //now make a histogram of the aers
  int X=origimage->GetLargestPossibleRegion().GetSize()[0];
  int Y=origimage->GetLargestPossibleRegion().GetSize()[1];
  int Z=origimage->GetLargestPossibleRegion().GetSize()[2];

  typedef itk::Point<double,AERHistSize> AERHistType;
  AERHistType aerhist;
  const unsigned int Npixels = X*Y*Z;

  entropyimage=EntropyImageType::New();
  CreateImage<EntropyImageType>( origimage->GetLargestPossibleRegion().GetSize(), origimage->GetSpacing(),entropyimage );

  int pixidx=0,histbin;
  double objdist;
  PointType loc;
  EntropyImageType::IndexType pxlidx;
  EntropyDataType entropyval;
  EntropyDataType wholeVolumeEntropy=0.,maxWholeVolumeEntropy=-1.,minWholeVolumeEntropy=9999999999999;
  int normval=0;
  double sigma=sqrt(distancescale*fibrethickness);
  double lengthpriorvar=3; 
  double lengthpriormean=7.5; //If sqrt(lengthpriormean) multiplied by +/- 3 then good fit lines approximate well but others don't
  double sigmasq=sigma*sigma;
  typedef MSHists<double> MSHistType;
  MSHistType mshists;
  for (int z=0;z<Z;z++) {
    //    std::cout <<z<<"/"<<Z-1<<":"<<wholeVolumeEntropy<<","<<minWholeVolumeEntropy<<","<<maxWholeVolumeEntropy<<std::endl;
    loc[2]=z*dims[2];
    pxlidx[2]=z;
    for (int y=0;y<Y;y++) {
      loc[1]=y*dims[1];
      pxlidx[1]=y;
      for (int x=0;x<X;x++) {
	loc[0]=x*dims[0];
	pxlidx[0]=x;
	entropyval=0.;
	if (origimage->GetPixel(pxlidx)>=bgThreshold) {
	  normval++;
	  //initialise histogram
	  mshists.Reset();
	  //compute histogram
	  for (int i=0;i<Nobjects;i++) {
	    objdist=DetermineDistance3Dseg(loc,x1points[i],x2points[i]);
	    if (objdist<vnl_math::eps)
	      objdist=vnl_math::eps;
	    if (objdist<sigma*3) {
	      //histbin=DiscretizeAER(aers[i],MSHistType::GetDelta(0));
	      histbin=DiscretizeAER(aers[i],radialscale);
	      mshists[0][histbin]+=exp(-pow(objdist,2)/(2*sigmasq));
	    }
	  }
	  entropyval=mshists.ComputeAllEntropyAllSteps();
	  //previously normalized so that the entropy calculation makes sense
	  //however the normalization removes the weighting that would come
	  //from the relative distances to the fibres encountered.
	  entropyval*=mshists.GetSum(0);
	}
	entropyimage->SetPixel(pxlidx,entropyval);
	if (entropyval>maxWholeVolumeEntropy)
	  maxWholeVolumeEntropy=entropyval;
	else if (entropyval<minWholeVolumeEntropy)
	  minWholeVolumeEntropy=entropyval;
	wholeVolumeEntropy+=entropyval;
      }
    }
  }
  wholeVolumeEntropy/=normval;
  //  return wholeVolumeEntropy;
  std::cout <<wholeVolumeEntropy<<","<<minWholeVolumeEntropy<<","<<maxWholeVolumeEntropy<<std::endl;
  return maxWholeVolumeEntropy;
}
////////////////////////////////////////////////////////////////////////////
EntropyDataType ObjectOrientations::DetermineMultiScaleEntropy(const MidImageType::Pointer origimage, const int bgThreshold, const double *dims) {
  assert(Nobjects>0);

  //now make a histogram of the aers
  int X=origimage->GetLargestPossibleRegion().GetSize()[0];
  int Y=origimage->GetLargestPossibleRegion().GetSize()[1];
  int Z=origimage->GetLargestPossibleRegion().GetSize()[2];

  typedef MS2DHists<double> MS2DHistType;
  MS2DHistType mshists;

  const unsigned int Npixels = X*Y*Z;

  entropyimage=EntropyImageType::New();
  CreateImage<EntropyImageType>( origimage->GetLargestPossibleRegion().GetSize(), origimage->GetSpacing(),entropyimage );

  int pixidx=0,histbin;
  double objdist;
  PointType loc;
  EntropyImageType::IndexType pxlidx;
  EntropyDataType entropyval;
  //  EntropyDataType wholeVolumeEntropy=0.,maxWholeVolumeEntropy=-1.;
  EntropyDataType wholeVolumeEntropy=0.,maxWholeVolumeEntropy=-1.,minWholeVolumeEntropy=9999999999999;
  int normval=0;
  double lengthpriorvar=3; 
  double lengthpriormean=7.5; //If sqrt(lengthpriormean) multiplied by +/- 3 then good fit lines approximate well but others don't
  for (int z=0;z<Z;z++) {
    //    std::cout <<z<<"/"<<Z<<std::endl;
    loc[2]=z*dims[2];
    pxlidx[2]=z;
    for (int y=0;y<Y;y++) {
      loc[1]=y*dims[1];
      pxlidx[1]=y;
      for (int x=0;x<X;x++) {
	loc[0]=x*dims[0];
	pxlidx[0]=x;
	entropyval=0.;

	if (origimage->GetPixel(pxlidx)>=bgThreshold) {
	  normval++;
	  //initialise histograms
	  mshists.Reset();
	  //compute histogram
	  for (int i=0;i<Nobjects;i++) {
	    objdist=DetermineDistance3Dseg(loc,x1points[i],x2points[i]);
	    if (objdist<vnl_math::eps)
	      objdist=vnl_math::eps;
	    mshists.AddToBin(aers[i],objdist,fibrethickness);
	  }
	  entropyval=mshists.ComputeAllEntropyAllSteps();
	  //previously normalized so that the entropy calculation makes sense
	  //however the normalization removes the weighting that would come
	  //from the relative distances to the fibres encountered.
	  entropyval*=mshists.GetMeanSum();
	}
	entropyimage->SetPixel(pxlidx,entropyval);
	if (entropyval>maxWholeVolumeEntropy)
	  maxWholeVolumeEntropy=entropyval;
	else if (entropyval<minWholeVolumeEntropy)
	  minWholeVolumeEntropy=entropyval;
	wholeVolumeEntropy+=entropyval;
      }
    }
  }
  //  std::cout <<"normval:"<<normval<<std::endl;
  wholeVolumeEntropy/=normval;
  //  std::cout <<wholeVolumeEntropy<<","<<minWholeVolumeEntropy<<","<<maxWholeVolumeEntropy<<std::endl;
  //  std::cout <<"wholeVolumeEntropy:"<<wholeVolumeEntropy<<std::endl;
  //  std::cout <<"maxWholeVolumeEntropy:"<<maxWholeVolumeEntropy<<std::endl;
  //  return wholeVolumeEntropy;
  entropyimage->Register();
  return maxWholeVolumeEntropy;
}
////////////////////////////////////////////////////////////////////////////
EntropyDataType ObjectOrientations::DetermineMultiScaleEntropy(void) {
  assert(Nobjects>0);

#if(0)
  typedef MS2DHists<double> MS2DHistType;
  MS2DHistType mshists;

  EntropyDataType entropyval;
  //initialise histograms
  mshists.Reset();
  //compute histogram
  for (int i=0;i<Nobjects;i++) 
    mshists.AddToBin(aers[i],1,1);
  entropyval=mshists.ComputeAllEntropyAllSteps();
  return entropyval;
#else
  typedef MSHists<double> MSHistsType;
  MSHistsType mshists;
  mshists.Reset();
  for (int i=0;i<Nobjects;i++)
    mshists.AddToBin(aers[i]);
  EntropyDataType entropyval=mshists.ComputeAllEntropyAllSteps();
  //  EntropyDataType entropyval=mshists.ComputeAllFactorAllSteps();
  return entropyval;
#endif
}
////////////////////////////////////////////////////////////////////////////
EntropyDataType ObjectOrientations::DetermineSOPfromEndPnts(const FRC::VOL_T &endpoints) {
  SOPfromAlignMatrix sopfromalignmatrix;
  return sopfromalignmatrix.CalcAll(endpoints);
}
////////////////////////////////////////////////////////////////////////////
void ObjectOrientations::DetermineEffFactorsfromEndPnts(const FRC::VOL_T &endpoints, std::vector<EntropyDataType> &efffacts) {
  OrientEfficiencyFactors efffactors;
  efffactors.CalcAll(endpoints,efffacts);
}
////////////////////////////////////////////////////////////////////////////
EntropyDataType ObjectOrientations::DetermineSingleScaleSOP(void) {
  assert(Nobjects>0);

  typedef MSHists<double> MSHistsType;
  MSHistsType mshists;
  mshists.Reset();
  for (int i=0;i<Nobjects;i++)
    mshists.AddToBin(aers[i]);
  //  EntropyDataType entropyval=mshists.ComputeAllEntropyAllSteps();
  EntropyDataType entropyval=mshists.ComputeAllFactorAllSteps();
  return entropyval;
}
////////////////////////////////////////////////////////////////////////////
void ObjectsSamples::convert(const BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter, const double *dims) {
  Nobjects=binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
  objectsvectors=FRC::VOL_T(Nobjects);

  ObjectOrientations::MeasurementVectorType mv;
  
  for(unsigned int i = 0; i < Nobjects; i++) {
    ObjectOrientations::SampleType::Pointer sample=ObjectOrientations::SampleType::New();
    // Get the ith region and extract out the samples that will be used for determining the sample vector
    BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    int objectsize=labelObject->Size();
    FRC::FIBRE_T objectpnts(objectsize);
    FRC::PNT_T frc_mv(ObjectOrientations::MeasurementVectorLength);;
    for(unsigned int pixelId = 0; pixelId < objectsize; pixelId++) {
      for (int l=0;l<ObjectOrientations::MeasurementVectorLength;l++) {
	mv[l]=labelObject->GetIndex(pixelId)[l]*dims[l];
	frc_mv[l]=mv[l];
      }
      sample->PushBack(mv);
      objectpnts[pixelId]=frc_mv;
    }
    objectsvectors[i]=objectpnts;
    //////////////////////
  }  
}
////////////////////////////////////////////////////////////////////////////
void ObjectsSamples::read(const std::string &filename) {
  objectsvectors=FRC::vectorfileread(filename);
  Nobjects=objectsvectors.size();
}
////////////////////////////////////////////////////////////////////////////
float ObjectOrientations::DetermineLineProportion(const MidImageType::Pointer origimage, const PointType _x1, const PointType _x2, const PixelType threshold) {
  MidImageType::SpacingType spacing=origimage->GetSpacing();
  MidImageType::SizeType size=origimage->GetLargestPossibleRegion().GetSize();
  MidImageType::IndexType x1,x2;
  for (int i=0;i<3;i++) {
    x1[i]=int(_x1[i]/spacing[i]+0.5);
    x2[i]=int(_x2[i]/spacing[i]+0.5);
  }
  itk::LineIterator<MidImageType> lineIt(origimage,x1,x2);
  MidImageType::PixelType imgval;
  lineIt.GoToBegin();
  float in=0.f,total=0.f;
  while (lineIt.GetIndex()[0]>=0 && lineIt.GetIndex()[0]<size[0] && lineIt.GetIndex()[1]>=0 && lineIt.GetIndex()[1]<size[1] && lineIt.GetIndex()[2]>=0 && lineIt.GetIndex()[2]<size[2] && !lineIt.IsAtEnd()) {
    imgval=lineIt.Value();
    if (imgval>threshold) {
      in++;
      //      lineIt.Set(200);
    }
    //    else
    //      lineIt.Set(200);
    total++;
    ++lineIt;
  }
  return in/total;
}
////////////////////////////////////////////////////////////////////////////
void ObjectOrientations::PruneLineSegments(const MidImageType::Pointer origimage, PixelType threshold) {
  int incount=0;
  VecPoints _x1points,_x2points,_means,_aers;
  float lineinoutproportion;
  PointType x1,x2;
  for (int i=0;i<Nobjects;i++) {
    for (int j=0;j<3;j++) {
      x1[j]=x1points[i][j];
      x2[j]=x2points[i][j];
    }
    lineinoutproportion=DetermineLineProportion(origimage,x1,x2,threshold);
    if (lineinoutproportion>=INOUTPROPORTION_THRESHOLD) {
      _x1points.push_back(x1points[i]);
      _x2points.push_back(x2points[i]);
      _aers.push_back(aers[i]);
      _means.push_back(means[i]);
      incount++;
    }
  }
  x1points=_x1points;
  x2points=_x2points;
  aers=_aers;
  means=_means;
  Nobjects=incount;
}
////////////////////////////////////////////////////////////////////////////
void ObjectOrientations::PruneLineSegmentsEccentricity(void) {
  int incount=0;
  VecPoints _x1points,_x2points,_means,_aers,_eigvals;
  float eigvalratio[2];
  PointType x1,x2;
  for (int i=0;i<Nobjects;i++) {
    for (int j=0;j<3;j++) {
      x1[j]=x1points[i][j];
      x2[j]=x2points[i][j];
    }
    eigvalratio[0]=eigvals[i][2]/eigvals[i][0];
    eigvalratio[1]=eigvals[i][2]/eigvals[i][1];
    if (eigvalratio[0]>EIGVALRATIO && eigvalratio[1]>EIGVALRATIO) {
      _x1points.push_back(x1points[i]);
      _x2points.push_back(x2points[i]);
      _aers.push_back(aers[i]);
      _means.push_back(means[i]);
      _eigvals.push_back(eigvals[i]);
      incount++;
    }
  }
  x1points=_x1points;
  x2points=_x2points;
  aers=_aers;
  means=_means;
  eigvals=_eigvals;
  Nobjects=incount;
}
////////////////////////////////////////////////////////////////////////////
void ObjectOrientations::WritePoints(const VecPoints &pntstowrite, const string &filename) {
  int N=pntstowrite.size();
  if (N>0) {
    std::ofstream fileout(filename.c_str());
    for (int i=0;i<pntstowrite.size();i++) {
      fileout <<pntstowrite[i][0]<<","<<pntstowrite[i][1]<<","<<pntstowrite[i][2]<<std::endl;
    }
    fileout.close();
  }
}
////////////////////////////////////////////////////////////////////////////
#if(0)
void LabelMapToImg(BinaryImageToLabelMapFilterType::Pointer labels, MidImageType::Pointer imgout) {
    for(unsigned int i = 0; i < labels->GetOutput()->GetNumberOfLabelObjects(); i++) {
      // Get the ith region
      BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labels->GetOutput()->GetNthLabelObject(i);
      // Output the pixels composing the region
      for(unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) 
	imgout->SetPixel(labelObject->GetIndex(pixelId),65535);
    }
}
#endif
