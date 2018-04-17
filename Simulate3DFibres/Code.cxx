#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "../itkCentralProcessing/itkhelpers.hxx"

#include "../itkCentralProcessing/vectorfileread.hxx"
#include "../itkCentralProcessing/itkprocessing.hxx"

#include "../itkCentralProcessing/itkfibresimulator.hxx"
#include "../itkCentralProcessing/vectorfileread.hxx"

int main(int argc, char *argv[] ) {
  typedef itk::ImageFileReader< MidImageType > ImageReaderType;

  if (argc < 19) {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "\n\t[1] dim0\n\t[2] dim1\n\t[3] dim2\n\t[4] spacing0\n\t[5] spacing1\n\t[6] spacing2\n\t[7] size0\n\t[8] size1\n\t[9] size2\n\t[10] upsample rate(>1)\n\t[11] interplane downsample rate(>1)\n\t[12] #fibres\n\t[13] max delta theta\n\t[14] max delta phi\n\t[15] theta quantization\n\t[16] phi quantization\n\t[17] fibrelengths\n\t[18] fibre thickness (radius in upsample pixels)\n\t[19] background pixel mean value\n\t[20] fibre pixel mean value\n\t[21] noise standard deviation \n\t[22]op filename" << std::endl;
    return EXIT_FAILURE;
  }

  //createimage from itkhelpers (at rel high res)
  //add fibres 
  //add noise
  //downsample
  //export

  //X,Y,Z
  int dims[3]={atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
  float spacing[3]={atof(argv[4]), atof(argv[5]), atof(argv[6])}; //inclusive of size
  float size[3]={atof(argv[7]), atof(argv[8]), atof(argv[9])};
  float upsamplerate=atof(argv[10]),interplanedownsamplerate=atof(argv[11]);
  int Nfibres=atoi(argv[12]);
  int maxdeltaangles[2]={atoi(argv[13]),atoi(argv[14])};
  int Nquantizations[2]={atoi(argv[15]),atoi(argv[16])};
  bool usefixedlength=true;
  float fibrelengths=atof(argv[17]);
  if (fibrelengths<varepsilon)
    usefixedlength=false;
  int fibrethickness=atoi(argv[18]);
  unsigned short backgroundpixelvalue=atoi(argv[19]);
  unsigned short fibrepixelvalue=atoi(argv[20]);
  double noiseparamval=pow(atof(argv[21]),2);
  string opfilename(argv[22]);
  MidImageType::RegionType::SizeType highres;
  MidImageType::SpacingType highres_spacing,lowres_spacing;
  for (int i=0;i<2;i++) {
    lowres_spacing[i]=spacing[i];
    highres[i]=int(dims[i]*upsamplerate);
    highres_spacing[i]=spacing[i]/upsamplerate;
  }
  lowres_spacing[2]=spacing[2];
  highres_spacing[2]=spacing[2]/interplanedownsamplerate;
  highres[2]=int(dims[2]*interplanedownsamplerate);

  MidImageType::Pointer highresimage=MidImageType::New();
  std::cout <<"Creating first image..."<<highres[0]<<","<<highres[1]<<","<<highres[2]<<","<<std::endl;
  CreateImage<MidImageType>(highres,highres_spacing,highresimage);

  std::cout <<"Created first image..."<<std::endl;
  ///////////////////////////////////////////////////
  //setting up random for both points and noise
  //////////////////////////////////////////////////  
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();
  ///////////////////////////////////////////////////
  //setting background values
  //////////////////////////////////////////////////  
  SetBackgroundValue(highresimage,backgroundpixelvalue);
  ///////////////////////////////////////////////////
  //adding fibres
  //////////////////////////////////////////////////  
  FRC::VOL_T fibres=AddFibres(highresimage,highres,highres_spacing,fibrepixelvalue,Nfibres,fibrelengths,maxdeltaangles,Nquantizations);
  FRC::vectorfilewrite("fibrepoints.csv",fibres);
  ///////////////////////////////////////////////////
  int RETURNVAL=0;
#if(1)
  //  MidImageType::Pointer linevol=SubSampler(highresimage,highres_spacing,upsamplerate,interplanedownsamplerate);
  MidImageType::Pointer linevol=DownSampleVolume(highresimage,size,highres_spacing,upsamplerate,interplanedownsamplerate);
  RETURNVAL=ExportVolume(linevol,"linevolume.mhd");
  if (RETURNVAL!=0)
    return RETURNVAL;
#endif
  ///////////////////////////////////////////////////
  //dilating fibres
  //////////////////////////////////////////////////  
  highresimage=DilateFibres(highresimage,fibrethickness,fibrepixelvalue);
  ///////////////////////////////////////////////////
  //setting noise values
  //////////////////////////////////////////////////  
  SetNoise(highresimage,noiseparamval);
  ///////////////////////////////////////////////////
  //downsampling
  //////////////////////////////////////////////////  
  MidImageType::Pointer lowresimage=DownSampleVolume(highresimage,size,highres_spacing,upsamplerate,interplanedownsamplerate);
  ///////////////////////////////////////////////////
  //saving volume
  //////////////////////////////////////////////////  
  RETURNVAL=ExportVolume(lowresimage,opfilename);

  return RETURNVAL;
}
