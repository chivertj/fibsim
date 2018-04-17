#ifndef ITKFIBRESIMULATOR_HEADER
#define ITKFIBRESIMULATOR_HEADER

#include "itkprocessing.hxx"
#include "itkhelpers.hxx"
#include "vectorfileread.hxx"

void  SetBackgroundValue(MidImageType::Pointer highresimage,unsigned short pixelval);
FRC::VOL_T AddFibres(MidImageType::Pointer highresimage,MidImageType::RegionType::SizeType superdims,MidImageType::SpacingType super_spacing,unsigned short fibrepixelvalue,int Nfibres,float fibrelengths, int *maxdeltaangles, int *Nanglequantizations);
FRC::VOL_T AddFibres(MidImageType::RegionType::SizeType superdims,MidImageType::SpacingType super_spacing,int Nfibres,float fibrelengths, int *maxdeltaangles, int *Nanglequantizations);
MidImageType::Pointer DilateFibres(MidImageType::Pointer highresimage,unsigned short radiusval,unsigned short fibrepixelvalue);
void SetNoise(MidImageType::Pointer highresimage, double noiseparamval);
MidImageType::Pointer DownSampleVolume(MidImageType::Pointer highresimage, float *lowressize, MidImageType::SpacingType highres_spacing, float upsamplerate, float interplanedownsamplerate);
int ExportVolume(const MidImageType::Pointer highresimage,const std::string &opfilename);

const float varepsilon=1e-5;

MidImageType::Pointer SubSampler(const MidImageType::Pointer filteredimage, MidImageType::SpacingType highres_spacing, float upsamplerate, float interplanedownsamplerate);

#endif //ITKFIBRESIMULATOR_HEADER
